#include <iostream>
#include <fstream>



#include "GraphicCardOpenCL.cpp"
#include "generatingSequenceFromFile.cpp"




int main(){
    std::ofstream fOutput;
    fOutput.open("output02.txt");
    int *aR;
    int *wR;
    myint nRAM;
    myint m;
    
    std::cout << "How long is the array w?" << std::endl;
    std::cin >> m;
    std::cout << "Enter array of "<<m <<" elements - 0 or 1 please" << std::endl;
    wR = new int[m];
    for (int j=0; j<m; j++) {
        std::cout<<"Enter the "<<j<<"th element. ";
        std::cin >> wR[j];
    }
    
    
    readSequenceFromFile("input02.txt", &aR, &nRAM );
    

    
    
    
    
    GraphicCard myCard("parallel02.cl");

    
    
    myCard.writeDeviceMemory("a",aR,nRAM);
    myCard.writeDeviceMemory("b",aR,nRAM);
    myCard.writeDeviceMemory("n",&nRAM,1);
    
    for (int k=0; k < m; k++) {
        if (wR[k] == 0) {
            std::cout<<"Executing kernel X"<<std::endl;
            myCard.executeKernel("X",nRAM);
            myCard.executeKernel("copyBtoA",nRAM);
        } else {
            std::cout<<"Executing kernel Y"<<std::endl;
            myCard.executeKernel("Y",nRAM);
            myCard.executeKernel("copyBtoA",nRAM);
        }
    }
    
    
    int * resultRAM;
    resultRAM=new int[nRAM];
    

    
    myCard.readDeviceMemory("a",resultRAM,nRAM);
    
    myint mLen=nRAM;
    if(mLen>30){mLen=30;}
    for(int i=0;i<mLen;i++){
        fOutput<<resultRAM[i]<<" ";
    }
    fOutput<<std::endl;
    
    
    fOutput.close();
    delete [] aR;
    delete [] resultRAM;

    

    return 1;
}
