# Introduction to parallel programming on graphic cards with opencl and c++
This repository contains the class GraphicCard that allows beginners to start programming with opencl/c++ without worrying about proper management of platforms, contexts, devices, kernels, workgroups, and memory objects.

It also contains a simple parallel random number generator.

## Compilation and usage of the class

###Step 1: Placing files in correct folders
Place the files GraphicCardOpenCL.cpp and gc_parallel_rnum.cl in the folder that will contain your source code and place the line 

  #include "GraphicCardOpenCL.cpp"

in the header of your file. 

In order to run examples given in the folder you will also need generatingSequenceFromFile.cpp. 

In order to compile and run example01.cpp, make sure input01.txt is also in the folder.

###Step 2: Running the compilers
#### Mac OSX 10.11.**
As long as you have xcode installed on your system, the following command in terminal would do the job: 

  c++ -framework opencl -o myprogram01 example01.cpp

After the compilation, you can run the program by typing (in terminal window)

  ./myprogram01

#### Ubuntu 16.04
If you have proper drivers for grahpic card installed, you need to type the following in the terminal
  c++ -o myprogram01 example01.cpp -lOpenCL -std=c++11

The compiler may return warnings if your version of opencl is 2.0. However, it should still compile well and the program can be run by typing
  ./myprogram01
