# G4Coil
G4Coil is a project to create a Coil like structure and testing it's scintillation properties in Geant4. 

# How to Run: (Setup) "I have to fix an issue with this process because the CMake file does not copy some of the required files to the build directory. The file are on the /project_scripts"
1.) Without the fast_build.sh

mkdir build
cd build
cmake ..
make -jN               N="Number of Threads on your computer"
./G4Coil

/run/beamOn N          N="Number of Event" "This step is also done when the program is running"

2.) With the fast_build.sh
. fast_build.sh

This will setup everything and run the simulation automatically. You change the run pararmeters in the file "/project_scripts/batch.mac" .
