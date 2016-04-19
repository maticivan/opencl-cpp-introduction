# Introduction to parallel programming on graphic cards with opencl and c++
This repository contains the class GraphicCard that allows beginners to start programming with opencl/c++ without worrying about proper management of platforms, contexts, devices, kernels, workgroups, and memory objects.

It also contains a simple parallel random number generator.

#GPGPU 
The abbreviation GPGPU refers to general purpose programming on graphics processing units. The graphic cards of modern computers can be used to do serious scientific calculations. They contain thousands of computing cores (that we will call processing elements), and as such they are ideal for parallel programming. The processing elements can be thought of as small CPUs that are numerous, but not as advanced as the central processors. High end CPUs can have speeds in range of 4GHz, but a quad-core system will have only four of these units. On the other hand, the graphic card can have 4000 processing elements each of which runs at 1GHz. In addition, not all operations are permitted and not all data structures are available to the processing elements of a typical graphic card. 

# Class GraphicCard
The interaction between GPU and the host platform is fairly complex for the beginners. This tutorial introduces the class GraphicCard written in C++ which does most of the work associated with memory management on GPU. The usage of this class removes the necessity of understanding completely the mechanisms behind the host/device interaction. 
