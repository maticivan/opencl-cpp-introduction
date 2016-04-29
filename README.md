# Introduction to parallel programming on graphic cards with opencl and c++
This repository contains the class GraphicCard that allows beginners to start programming with opencl/c++ without worrying about proper management of platforms, contexts, devices, kernels, workgroups, and memory objects.

It also contains a simple parallel random number generator.

## Compilation and usage of the class

###Step 1: Place the files GraphicCardOpenCL.cpp and gc_parallel_rnum.cl in the folder that will contain your source code and place the line 

  #include "GraphicCardOpenCL.cpp"

in the header of your file. 

In order to run examples given in the folder you will also need generatingSequenceFromFile.cpp. 

In order to compile and run example01.cpp, make sure input01.txt is also in the folder and type the following in Mac OSX terminal window:

  c++ -framework opencl -o myprogram01 example01.cpp

After the compilation, you can run the program by typing

  ./myprogram01





