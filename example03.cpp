#include <iostream>
#include <fstream>



#include "GraphicCardOpenCL.cpp"
#include "generatingSequenceFromFile.cpp"




int main(){
    std::ofstream fOutput;
    fOutput.open("output03.txt");
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
    
    
    readSequenceFromFile("input03.txt", &aR, &nRAM );
    

    
    
    
    
    GraphicCard myCard("parallel03.cl");

    
    
    myCard.writeDeviceMemory("aDM",aR,nRAM);
    myCard.writeDeviceMemory("bDM",aR,nRAM);
    myCard.writeDeviceMemory("nDM",&nRAM,1);
                             
    std::string* parameters;
    parameters=new std::string[3];
    parameters[0]="aDM";
    parameters[1]="bDM";
    parameters[2]="nDM";
    std::string helpS;
    
    for (int k=0; k < m; k++) {
        if (wR[k] == 0) {
            myCard.setKernelArguments("X",parameters,3);
            myCard.executeKernel("X",nRAM);

        } else {
            myCard.setKernelArguments("Y",parameters,3);
            myCard.executeKernel("Y",nRAM);
        }
        helpS=parameters[0];
        parameters[0]=parameters[1];
        parameters[1]=helpS;
    }
    
    
    int * resultRAM;
    resultRAM=new int[nRAM];
    

    
    myCard.readDeviceMemory(parameters[0],resultRAM,nRAM);
    
    myint mLen=nRAM;
    if(mLen>30){mLen=30;}
    for(int i=0;i<mLen;i++){
        fOutput<<resultRAM[i]<<" ";
    }
    fOutput<<std::endl;
    
    
    fOutput.close();
    delete [] aR;
    delete [] resultRAM;
    delete [] parameters;
    

    return 1;
}
