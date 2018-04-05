#
# Program to put the ct-slices of porous structure together to get a 3d array.
# It can easily be implemented in the LB code,
# since material-domain-voxels and fluid-domain-voxels can be recognised in a binary array
#
# The array is written to a .vtk file to plot the array with paraview <- mavi can only plot small arrays
#

"""
.vtk file is also restricted to a certain number of slices: 250 slices is barely working. I do not know, if it is my computer, that is so bad, or in general, it is
quite costly to visualise voxel pictures

"""

import numpy
from pylab import *
import matplotlib.pyplot as plt
from PIL import Image                       # image processing
#from mpl_toolkits.mplot3d import Axes3D    # 3d plotting stuff
import mayavi.mlab as mavi                  # better 3d plotting stuff
from tvtk.api import tvtk, write_data       # export to .vtk file


# total number of slices
numberOfPics = 1180#191
picList = numpy.arange(numberOfPics)


# get array "obstacle" with dimension (numberOfPics x length x width)
# length and width are obtained from a test picture (dummyPic) since all slices have the same shape
dummyPic = Image.open("0002_00_bin_01071.tif")
dummyArray = numpy.array( dummyPic, dtype = float)

obstacle = numpy.zeros((numberOfPics, dummyArray.shape[0], dummyArray.shape[1]))

#
# convert each ct slice into binary array and put them together in the big 3d obstacle array
#
for ii in range(numberOfPics):

    # read image
    pic = Image.open("0002_00_bin_01071.tif")
    pic = Image.open("0002_00_bin_"+str(ii+71).zfill(5)+".tif")    
    # convert into grey scale picture (actually not necessary, since it is already a black and white pic)
    pic = pic.convert('L')


    # get the grey values and asign black pixel to 0 and white pixel to 1
    # 1 and 0 interchangend 
    values = [1 if v > 128 else 2 for v in range(256)]
    obstacletmp = pic.point(values)
    obstacle[ii] = numpy.array(obstacletmp)     # at this point, the big obstacle array is filled with 1 and 0



centerX = 450
centerY = 450
radius = 370

halbachseA = 345.
halbachesB = 325.



for kk in range(numberOfPics):
    print "kk = ", kk
    for ii in range(obstacle.shape[1]):
        for jj in range(obstacle.shape[2]):
            if ((ii-centerX)/halbachseA)**2 + ((jj-centerY)/halbachesB)**2 > 1.:
                obstacle[kk,ii,jj] = 0
            if numpy.abs(((ii-centerX)/halbachseA)**2 + ((jj-centerY)/halbachesB)**2 - 1.) < 1e-2 :
                obstacle[kk,ii,jj] = 5

oo = obstacle[0,centerX-halbachseA:centerX+halbachseA, centerY-halbachesB:centerY+halbachesB].copy()
#oo = obstacle[0]
plt.contourf(oo)
plt.axis('off')
plt.grid()
plt.colorbar()
plt.show()

domain = numpy.zeros((numberOfPics, oo.shape[0], oo.shape[1]))
for kk in range(numberOfPics):
    domain[kk,:,:] = obstacle[kk,centerX-halbachseA:centerX+halbachseA, centerY-halbachesB:centerY+halbachesB]








# inlet slice
inletSliceTmp = numpy.ones((obstacle.shape[1], obstacle.shape[2]))
for ii in range(inletSliceTmp.shape[0]):
    for jj in range(inletSliceTmp.shape[1]):
        if ((ii-centerX)/halbachseA)**2 + ((jj-centerY)/halbachesB)**2 > 1.:
            inletSliceTmp[ii,jj] = 0
        if numpy.abs(((ii-centerX)/halbachseA)**2 + ((jj-centerY)/halbachesB)**2 - 1.) < 1e-2 :
            inletSliceTmp[ii,jj] = 5
            
inletSlice = inletSliceTmp[centerX-halbachseA:centerX+halbachseA, centerY-halbachesB:centerY+halbachesB]
            
    



def writeToVTK(filename, *array1):
    print "writing vtk ..."
    def namestr(obj, namespace):
        return [name for name in namespace if namespace[name] is obj]       
    name = namestr(array1[0], globals())       
    #write data to .vtk file, in order to plot with paraview
    vtkFile = tvtk.ImageData(spacing = (1, 1, 1), origin = (0, 0, 0)) # get an empty "vtk array"
    vtkFile.point_data.scalars = numpy.ravel(array1[0], order='F') # flattens the obstalce array in C-style and writes the values to the vtk array
    vtkFile.point_data.scalars.name = str(name) # shown in paraview
    vtkFile.dimensions = array1[0].shape
    write_data(vtkFile,str(filename)+str(array1[0].shape[0])+'_'+str(array1[0].shape[1])+'_'+str(array1[0].shape[2])+'.vtk') 
    print "writing vtk ... ok\n"
    
def writeLBGeometry(array1):
    print "writing vtk for LB ..."  
    #write data to .vtk file, in order to plot with paraview
    vtkFile = tvtk.ImageData(spacing = (1, 1, 1), origin = (0, 0, 0)) # get an empty "vtk array"
    vtkFile.point_data.scalars = numpy.ravel(array1, order='A') 
    vtkFile.point_data.scalars.name = 'LBgeometry' # shown in paraview
    vtkFile.dimensions = array1.shape
    write_data(vtkFile,'LBgeometry_'+str(array1.shape[0])+'_'+str(array1.shape[1])+'_'+str(array1.shape[2])+'.vtk') # write to a .vtk file
    print "writing vtk for LB ... ok"
    
    
    
    
def makeInletOutlet(array, firstSlice):
    inletLength = 5
    outletLength = 5
    longArray= numpy.zeros((array.shape[0]+ inletLength+outletLength, array.shape[1], array.shape[2] ))
    longArray[inletLength:longArray.shape[0]-outletLength,:,:] = array
    
    
    #firstSlice = array[:,:,0].copy()
    lastSlice = array[array.shape[0]-1,:,:].copy()
            
    for ii in range(inletLength):
        longArray[ii,:,:] = firstSlice
    for ii in range(outletLength):            
        longArray[-ii-1,:,:] = lastSlice
    
    # inlet outlet 
    for ii in range(array.shape[1]):
        for jj in range(array.shape[2]):
            if longArray[0,ii,jj] ==1:
                longArray[0,ii,jj] = 3
            if longArray[-1,ii,jj] ==1:
                longArray[-1,ii,jj] = 4
    
    return longArray
    


newDomain = makeInletOutlet(domain, inletSlice)


""" newDomain is our target array """


def shrinkArray(array):
    # cut the array, such that we can divide the number by 4
    print array.shape
    print array.shape[0] / 2.
    print array.shape[1] / 4.
    print array.shape[2]%4.
    for ii in range(4):
        if (array.shape[2] - ii)%4 == 0:
            print array.shape[2] - ii, ii
            if ii != 0:
                cut = ii
            else:
                cut = 0
                
    if cut != 0:
        cutArray = array[:,0:-cut,0:-cut]
    else:
        cutArray = array.copy()
       
    
    print cutArray.shape
    
    # get empty array with the proper size
    shrunkenArray = numpy.zeros((cutArray.shape[0]/2., cutArray.shape[1]/2., cutArray.shape[2]/2.))
    
    print "shrunkenArray shape", shrunkenArray.shape
    
    # the voxel at 0,0,0 consists of 0,0,0, 0+1,0,0, 0,0+1,0,  0+1,0+1,0, 0,0+1,0+1, 0+1,0+1,0+1, 0+1,0,0+1, 0,0,0+1
    # for the next voxel at 1,0,0 we have to shift the original array by 2
    for ii in range(shrunkenArray.shape[1]):
        for jj in range(shrunkenArray.shape[2]):
            for kk in range(shrunkenArray.shape[0]):
                
                neighbours = []
                neighbours.append( [cutArray[2*kk,2*ii,2*jj], cutArray[2*kk,2*ii,2*jj+1], 
                                    cutArray[2*kk,2*ii+1, 2*jj], cutArray[2*kk,2*ii+1,2*jj+1],
                                    cutArray[2*kk+1,2*ii,2*jj], cutArray[2*kk+1,2*ii,2*jj+1], 
                                    cutArray[2*kk+1,2*ii+1, 2*jj], cutArray[2*kk+1,2*ii+1,2*jj+1]] )
                
                pixelValue =  numpy.array(neighbours).max()
                shrunkenArray[kk,ii,jj] = pixelValue
    # outlet for badArry, since it got cut 
    for ii in range(shrunkenArray.shape[1]):
        for jj in range(shrunkenArray.shape[2]):
            if shrunkenArray[-1,ii,jj] ==1:
                shrunkenArray[-1,ii,jj] = 4
    
    oo = shrunkenArray[-1]
    plt.contourf(oo)
    plt.axis('off')
    plt.grid()
    plt.colorbar()
    plt.show()  
    
    return shrunkenArray


print newDomain.shape
badResolution = shrinkArray(newDomain)
print "neues array", badResolution.shape

badBadResolution = shrinkArray(badResolution)
print badBadResolution.shape

writeToVTK("withInlet_", newDomain)
writeLBGeometry(newDomain)
writeToVTK("badResolutionWithInlet_", badResolution)
writeLBGeometry(badResolution)
writeToVTK("badBadResolutionWithInlet_", badBadResolution)
writeLBGeometry(badBadResolution)




print "done"