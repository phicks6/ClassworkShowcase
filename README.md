# Classwork Showcase
The propose of the repository is to showcase a number of the more interesting labs/projects that I've done over my undergrad and grad career.

## Parallel Matrix Multiplication
In this project I programmed Cannon's algorithm and the Parallel DNS algorithm to multiple matrixs. Cannon's algorithm works by splitting the computations of a matrix multiplication evenly across a 2D square mesh of processors. The Parallel DNS algorithm works by splitting the computation across a 3D cube of stacked meshes of processors. There is a much higher comunication cost for Parallel DNS as distributing the neccessary parts of the two matrixes takes longer when it needs to travel in 3 dimensions. However, if communications are quick enough then the Parallel DNS scales better.

### Running
MPI is required to run these programs. To compile them on linux: ```mpicc "cannon.c/DNS.c" -lm```  
To run: ```mpirun -np "NumProcessors" ./a.out "MatrixSize"```  
With "NumProcessors" being a perfect square for cannon.c or a perfect cube for DNS.c and "MatrixSize" being some multiple of "NumProcessors"  

## CPU Cache Simulator
In this project I made a CPU cache simulator that also simulates other aspects of the memory hierarchy. It simulates a 2 level cache with the second level being toggleable. Each level of cache has two writing policies: Write-Through/No Write Allocate and Write-Back/Write Allocate. Input can be physical addresses or virtual addresses that go through a simulated page table and option TLB. The simulator reads a config file named ```trace.config``` to set sizes, policies and toggles.  

### Running
The simulator can be made with the included makefile, which produces ```memhier```. It can be ran with some example input with ```./memhier < trace.dat```
