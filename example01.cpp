// Compile by typing "c++ -framework opencl -o myprogram01 example01.cpp" in the terminal window
// After that you can run your code by typing "./myprogram01" in the terminal window

// This program that generates a sequence based on "input01.txt" and doubles each of its elements.



#include <iostream>
#include <fstream>


#include "GraphicCardOpenCL.cpp"
#include "generatingSequenceFromFile.cpp"


int main(){
    
    int theNumberOfElementsToOutput=30;
    
    int *aRAM;
    int nRAM;
    
    readSequenceFromFile("input01.txt", &aRAM, &nRAM );
    
    // Allocate an object myCard that is of the type GraphicCard
    // myCard object will contain the kernel written in parallel01.cl
    GraphicCard myCard("parallel01.cl");
    
    // We will place the sequence aRAM to the device memory.
    // The kernel will look for a sequence "a" on the device memory.
    // We will create this sequence by copying the sequence aRAM from the host memory.
    myCard.writeDeviceMemory("a",aRAM,nRAM);
    // Similarly, we need to place the variable "n" on the graphic card.
    // Since the only allowed objects on the graphic cards are sequences, and n is
    // a single element, it must be treated as an one-element sequence.
    // Luckily, in C++ every pointer is an one-element sequence, and so is &nRAM.
    myCard.writeDeviceMemory("n",&nRAM,1);
    
    
    // The GPU is given instruction to execute the kernel on at least nRAM processing elements.
    // Caution: GPU may actually invoke more than nRAM elements, and you can't control this.
    //          Hardware design of GPU organizes the processing elements in groups and only the entire
    //          group can be called.
    myCard.executeKernel("eachTermGetsDoubled",nRAM);
    
    
    
    // After the kernel execution we need to retrieve the new modified sequence "a".
    // The sequence has to be copied from the device memory back to the RAM memory.
    // We can accomplish this by first creating the space in RAM sufficient to store the new sequence.
    int * resultRAM;
    resultRAM=new int[nRAM];
    
    // Now since the space is created, we call the function that reads the content of the
    // device memory and stores it in the sequence resultRAM.
    myCard.readDeviceMemory("a",resultRAM,nRAM);
    
    //The final step is to output the result in the file output01.txt
    std::ofstream fOutput;
    fOutput.open("output01.txt");

    int mLen=nRAM;
    if(mLen>theNumberOfElementsToOutput){mLen=theNumberOfElementsToOutput;}
    for(int i=0;i<mLen;i++){
        fOutput<<resultRAM[i]<<" ";
    }
    fOutput<<std::endl;
    
    
    fOutput.close();
    delete [] aRAM;
    delete [] resultRAM;
    
    
    
    return 1;
}
