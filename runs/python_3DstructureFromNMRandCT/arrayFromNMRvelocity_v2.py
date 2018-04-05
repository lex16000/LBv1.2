import numpy
from pylab import *
import matplotlib.pyplot as plt
from PIL import Image  
#import vtk                     # image processing
from tvtk.api import tvtk, write_data       # export to .vtk file
import os, sys
import scipy.io as sio

class NMRdata:
   

    def __init__(self, name):

        self.name = name
   
    @staticmethod
    def initialiseFromDatFiles(listOfFiles, direction, ending):
        print "loading slices from txt arrays ..."
        listOfSlices = []
        for files in listOfFiles:
            if files.endswith(ending):
                if files.startswith(direction):
                    listOfSlices.append(files)               
        print "loading slices from txt arrays ... ok"
        numberOfPics = len(listOfSlices) 
        print "number of slices: ", numberOfPics
        
        
        dummyArray = numpy.loadtxt(listOfSlices[0],  dtype=float, delimiter=',')
        print "shape of one slice: ", dummyArray.shape        
        # bigArray: directions. x y z
        velocityArray = numpy.zeros((dummyArray.shape[0],dummyArray.shape[1],numberOfPics))
        # go through v_ii and store in bigArray
        
        for ii in range(len(listOfSlices)):
            velocityArray[:,:,ii] = numpy.loadtxt(sorted(listOfSlices)[ii], dtype=float, delimiter=',')
        
        return velocityArray
   
        
    @staticmethod
    def writeToVTK(filename, *array1):
        print "writing vtk ...", str(filename)
        def namestr(obj, namespace):
            return [name for name in namespace if namespace[name] is obj]       
        name = namestr(array1[0], globals())       
        #write data to .vtk file, in order to plot with paraview
        vtkFile = tvtk.ImageData(spacing = (1, 1, 1), origin = (0, 0, 0)) # get an empty "vtk array"
        vtkFile.point_data.scalars = numpy.ravel(array1[0], order='F') # flattens the obstalce array in C-style and writes the values to the vtk array
        vtkFile.point_data.scalars.name = str(name) # shown in paraview
        vtkFile.dimensions = array1[0].shape
        # add second point data field
        numberOfArrays = len(array1)
        if len(array1) > 1:
            for ii in range(1,numberOfArrays):
                name = namestr(array1[ii], globals())
                vtkFile.point_data.add_array(numpy.ravel(array1[ii], order='F'))
                vtkFile.point_data.get_array(ii).name = str(name)
                vtkFile.point_data.update()
    
        write_data(vtkFile,str(filename)+'.vtk') 
        print "writing vtk ...", str(filename), "... ok"
  
    @staticmethod
    def getGeometry(array):
        print "constructing geometry from velocity values ..."
        newarray = numpy.zeros((array.shape[0],array.shape[1],array.shape[2]))
        for ii in range(array.shape[0]):
            for jj in range(array.shape[1]):
                for kk in range(array.shape[2]):
                    if numpy.abs(array[ii,jj,kk])  < 1e-15:
                        newarray[ii,jj,kk] = 0
                    else:
                        newarray[ii,jj,kk] = 1
        print "constructing geometry from velocity values ... ok"
        return newarray
        
        
    @staticmethod
    def writeLBGeometry(array1):
        print "writing LB vtk ..."
        #write data to .vtk file, in order to plot with paraview
        vtkFile = tvtk.ImageData(spacing = (1, 1, 1), origin = (0, 0, 0)) # get an empty "vtk array"
        vtkFile.point_data.scalars = numpy.ravel(array1, order='A') # flattens the obstalce array in C-style and writes the values to the vtk array
        vtkFile.point_data.scalars.name = 'LBgeometry' # shown in paraview
        vtkFile.dimensions = array1.shape
        write_data(vtkFile,'LBgeometry.vtk') # write to a .vtk file
        print "writing LB vtk ... ok"
        #print "writing to txt file ..."
        #raveled = numpy.ravel(array1, order = "C")
        #numpy.savetxt("LBgeometry.txt", raveled)
        #print "writing to txt file ... ok"
        
        
        
    @staticmethod
    def makeBoundary(domain):
        print "be aware of array: loop begins at node 1, not 0!"
        print "standard materialNumber: 2 for noSlip, 3 for inlet, 4 for outlet"
        for ii in range(1,domain.shape[0]-1):
            for jj in range(1, domain.shape[1]-1):
                #for kk in range(1,domain.shape[2]-1):
                for kk in range(domain.shape[2]):    
                    neighbourList = []
                    neighbourList.append(domain[ii+1,jj,kk])
                    neighbourList.append(domain[ii-1,jj,kk])
                    
                    neighbourList.append(domain[ii,jj+1,kk])
                    neighbourList.append(domain[ii,jj-1,kk])
                    
                    #neighbourList.append(domain[ii,jj,kk+1])
                    #neighbourList.append(domain[ii,jj,kk-1])
                    neighbourList.append(domain[ii,jj,kk])
                    neighbourList.append(domain[ii,jj,kk])
                    
                    neighbourList.append(domain[ii+1,jj+1,kk])
                    neighbourList.append(domain[ii-1,jj+1,kk])
                    neighbourList.append(domain[ii+1,jj-1,kk])
                    neighbourList.append(domain[ii-1,jj-1,kk])

                    #neighbourList.append(domain[ii,jj+1,kk+1])
                    #neighbourList.append(domain[ii,jj-1,kk+1])
                    #neighbourList.append(domain[ii,jj+1,kk-1])
                    #neighbourList.append(domain[ii,jj-1,kk-1])

                    #neighbourList.append(domain[ii+1,jj,kk+1])
                    #neighbourList.append(domain[ii-1,jj,kk+1])
                    #neighbourList.append(domain[ii+1,jj,kk-1])
                    #neighbourList.append(domain[ii-1,jj,kk-1])        
                    neighbourList.append(domain[ii,jj+1,kk])
                    neighbourList.append(domain[ii,jj-1,kk])
                    neighbourList.append(domain[ii,jj+1,kk])
                    neighbourList.append(domain[ii,jj-1,kk])

                    neighbourList.append(domain[ii+1,jj,kk])
                    neighbourList.append(domain[ii-1,jj,kk])
                    neighbourList.append(domain[ii+1,jj,kk])
                    neighbourList.append(domain[ii-1,jj,kk])        
                    
                    
                    xPoints = numpy.unique(neighbourList)
                    
                    
                    if (1 in xPoints) and (0 in xPoints):
                        domain[ii,jj,kk] = 2
                        
        
        for ii in range(domain.shape[0]):
            for jj in range(domain.shape[1]):  
                if numpy.abs(domain[ii,jj,0]) == 1:
                    domain[ii,jj,0] = 3
                if numpy.abs(domain[ii,jj,domain.shape[2]-1]) == 1:
                    domain[ii,jj,domain.shape[2]-1] = 4
        
        
                    
    @staticmethod
    def makeInletOutlet(array):
        inletLength = 3
        outletLength = 5
        longArray= numpy.zeros((array.shape[0], array.shape[1], array.shape[2] + inletLength+outletLength))
        longArray[:,:,inletLength:longArray.shape[2]-outletLength] = array
        
        
        firstSlice = array[:,:,0].copy()
        lastSlice = array[:,:,array.shape[2]-1].copy()
                
        for ii in range(inletLength):
            longArray[:,:,ii] = firstSlice
        for ii in range(outletLength):            
            longArray[:,:,-ii-1] = lastSlice
        
        return longArray
        
 
        



"""
path = os.getcwd()
listofFiles = os.listdir(path)

data = NMRdata("simple_tube")
vz = data.initialiseFromDatFiles(listofFiles,"vz", "dat")
vx = data.initialiseFromDatFiles(listofFiles,"vx", "dat")
vy = data.initialiseFromDatFiles(listofFiles,"vy", "dat")
vMag = numpy.sqrt(vz ** 2 + vx **2 + vy **2)
"""


datah = sio.loadmat("cmpflowresults20170907_VENC15_R0.2_alex.mat")
vx =  datah["vx"]
vy =  datah["vy"]
vz =  datah["vz"]

vMagn = numpy.sqrt(vz ** 2 + vy **2 + vx **2)

#geometry = materialnumber(vMagn)
data = NMRdata("porousTube")


geometry = data.getGeometry(vMagn)

cutX = [20,82]
cutY = [20,80]
cutZ = [20,160]

cutGeometry = geometry[cutX[0]:cutX[1],cutY[0]:cutY[1],cutZ[0]:cutZ[1]]


cutvx = vx[cutX[0]:cutX[1],cutY[0]:cutY[1],cutZ[0]:cutZ[1]]
cutvy = vy[cutX[0]:cutX[1],cutY[0]:cutY[1],cutZ[0]:cutZ[1]]
cutvz = vz[cutX[0]:cutX[1],cutY[0]:cutY[1],cutZ[0]:cutZ[1]]
cutvMagn = vMagn[cutX[0]:cutX[1],cutY[0]:cutY[1],cutZ[0]:cutZ[1]]

superNMRgeometry = numpy.zeros((4,cutvx.shape[0],cutvx.shape[1],cutvx.shape[2]))
superNMRgeometry[0] = cutvx
superNMRgeometry[1] = cutvy
superNMRgeometry[2] = cutvz
superNMRgeometry[3] = cutvMagn


LBgeometry = data.makeInletOutlet(cutGeometry)
compareVx = data.makeInletOutlet(superNMRgeometry[0])
compareVy = data.makeInletOutlet(superNMRgeometry[1])
compareVz = data.makeInletOutlet(superNMRgeometry[2])
compareVmagn = data.makeInletOutlet(superNMRgeometry[3])



data.writeToVTK("nmrData", geometry,vx,vy,vz,vMagn)
data.writeToVTK("NMRdata_LBsetup",compareVx,compareVy,compareVz,compareVmagn )

#data.makeBoundary(LBgeometry)
data.writeLBGeometry(LBgeometry)

print LBgeometry.shape

def countArea(array):
    countlist = []
    for kk in range(1,array.shape[2]-1):
        count = 0
        for ii in range(array.shape[0]):
            for jj in range(array.shape[1]):
                if array[ii,jj,kk] == 1:
                    count += 1
        countlist.append([kk, count])
    return countlist


#area = numpy.array(countArea(LBgeometry))
#plt.plot(area[:,0], area[:,1])
#plt.show()
#numpy.savetxt("crossSectionAreaAlongTube.txt", area)




print "done"