# soil_compression
Simulation of soil compression mechanical test using discrete element method (soft particles interaction).

1) The purpose of the program.

Simulation of mechanical tests of soil samples. The implemented method is the descrete element method with soft particles interactions (the method is similar to the molecular dynamics simulation).
In the current version, the only implemented test is a soil compression test (square piston loading on a soil sample of a square cuboid shape).

2) Installation.

No installation is required, the program (executable file dem.exe) can be placed in any convinient location. 


4) Usage.

4.1) General usage.

The program is a project oriented. There are two types of projects: the 'sample' project and the 'simulation' project.
The purpose of the 'sample' project is to fabricate a soil sample possesed a specific soil properties. The fabricated soil sample is ready for the mechanical testing where it can be reused a multiple times in different simulation studies.
The purpose of 'simulation' project is, by using the fabricated sample, perform the mechanical tests.

4.2) Command line.

To run the 'sample'project: dem.exe parameter_file.dat
To run the 'simulation' project: dem.exe parameter_file.dat fabricated_sample.smpl

4.3) Program output.

The 'sample' progect generates a 'sample' folder with a following content:
a) 'vtk_files' - the folder which consists of *.vtk files (if a 'record output' option was set to TRUE).
b) 'Equilibrium history' - the file (time stamps of convergence to equilibrium of a soil sample).
c) 'sample.smpl' - the file is a fabricated sample. Should be used in the 'simulation' project.

The 'simulation' project generates a 'simulation' folder with a following content:
a) 'vtk_files' - the folder which consists of *.vtk files (if a 'record output' option set to TRUE).

In addition to the 'sample' and 'simulation' folders, each program run generates the 'DEM_RUNTIME' file.
This file dynamically collects the 'DEM' program runtime activity, therefore, can be traced by user continuously during the program run.

4.4) Input data.

See the example of parameters file, it is self descriptive. Note, no extra lines allowed in this file.

5) Some practical suggestions.

a) It is not always necessary to record an output in a 'sample' project. If a 'record output' option will be set to FALSE,
computations will complete faster.
b) The input parameters values (in the parameter's file) in two consecutive (sample->simulation) projects allowed to be slightely different, but not for those parameters which determine grains sizes. Note, a sample equilibrium in a 'simulation' project is guarantied only if the parameters files in both types of projects are equal. 
