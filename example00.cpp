// Compile by typing "c++ -framework opencl -o myprogram00 example00.cpp" in the terminal window
// After that you can run your code by typing "./myprogram00" in the terminal window

// The following program is an illustration on how to use readSequenceFromFile

// The file input00.txt contains the number 2000 in the first line. This number indicates that
// we wish to form a sequence that contains 2000 elements.
// After the number 2000, the file input00.txt contains several integers, whose number is clearly less than 2000.
// These will be the first few elements of our sequence.
// Once the function reads the number -9, it will stop reading from the file and start generating random numbers
// until the total number of elements in the sequence is 2000.


// After generating the sequence, the program immediately outputs it to the file output00.txt
// The program will output only the first 30 elements. Of course, if you really want to see these random numbers
// you may change 30 to some other number of your choosing.


#include <iostream>
#include <fstream>

#include "GraphicCardOpenCL.cpp"
#include "generatingSequenceFromFile.cpp"



int main(){
    
    int theNumberOfElementsToOutput=30;
    
    int *aRAM;
    int nRAM;
    
    // The following command does all the magic.
    readSequenceFromFile("input00.txt", &aRAM, &nRAM );
    // The sequence is now stored in the variable aRam. Its length is stored in nRAM.
    
    // The final step is to output the result in the file output00.txt
    std::ofstream fOutput;
    fOutput.open("output00.txt");
    
    int mLen=nRAM;
    if(mLen>theNumberOfElementsToOutput){mLen=theNumberOfElementsToOutput;}
    for(int i=0;i<mLen;i++){
        fOutput<<aRAM[i]<<" ";
    }
    fOutput<<std::endl;
    
    
    fOutput.close();
    
    // It remains to free the memory allocated to the sequence aRAM.
    delete [] aRAM;

    return 1;
}
