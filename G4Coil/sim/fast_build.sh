# . G4Setup.sh

rm -rf build
mkdir build

cp /home/cristmor/Desktop/Research/Muon/G4Coil/sim/models/G4CoilSingleTurn.stl build
cp /home/cristmor/Desktop/Research/Muon/G4Coil/sim/project_scripts/batch.mac build
# cp /home/cristmor/Desktop/Research/Muon/G4Coil/sim/project_scripts/electron.mac build
cp /home/cristmor/Desktop/Research/Muon/G4Coil/sim/project_scripts/exampleB4a.out build
cp /home/cristmor/Desktop/Research/Muon/G4Coil/sim/project_scripts/exampleB4.in build
cp /home/cristmor/Desktop/Research/Muon/G4Coil/sim/project_scripts/gui.mac build
cp /home/cristmor/Desktop/Research/Muon/G4Coil/sim/project_scripts/init_vis.mac build
cp /home/cristmor/Desktop/Research/Muon/G4Coil/sim/project_scripts/plotHisto.C build
cp /home/cristmor/Desktop/Research/Muon/G4Coil/sim/project_scripts/plotNtuple.C build
cp /home/cristmor/Desktop/Research/Muon/G4Coil/sim/project_scripts/run1.mac build
cp /home/cristmor/Desktop/Research/Muon/G4Coil/sim/project_scripts/run2.mac build
cp /home/cristmor/Desktop/Research/Muon/G4Coil/sim/project_scripts/vis.mac build

cd build
cmake ..
make -j16

./G4Coil -t 1 -m batch.mac > cout.txt
# ./G4Coil -t 1 -m electron.mac > cout.txt

./G4Coil > cout_1.txt
cd ..
