#!/bin/bash

# Definiere Wertebereich für den zu untersuchenden Parameter
# laplace: r = for i in $(seq 8 2 20)
for i in $(seq 8 2 20)
do echo $i

#Erstelle Rechnungsordner (Hier: Gamma0.1, Gamma0.2, ...)
#mkdir Gamma$i
#cd Gamma$i
#Kopiere Input Daten rein
#cp ../Alternating_Ladder.py .
#cp ../SL1D_in.hsd .
#cp ../SL1D_N.py .
#cp ../Stacked_Ladder.py .
#Erstelle 2 leere Ordner, die DFTB braucht ... kannst du vermutlich löschen
#mkdir GS
#mkdir contacts

cp droplet3D_movingBoundaries.cpp.org_laPlace droplet3D_movingBoundaries.cpp

#Sucht in der Datei SL1D_Npy nach dem string 'GammaVal' und ersetzt ihn durch den Wert von i
#sed -i "s/shearVel/$i/g" droplet3D_movingBoundaries.cpp
sed -i "s/radiusFactor/$i/g" droplet3D_movingBoundaries.cpp

# Starte Rechnung
make && ./droplet3D_movingBoundaries
sleep 2



done
