



#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif


#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <iterator>
#include <regex>
#include <string>
#include <streambuf>
#include <random>
#include <chrono>







typedef int mysint;
typedef  int myint;




class Strin_gList{
public:
    std::string content;
    Strin_gList* next;
};

void deleteStrin_gListInGC(Strin_gList* first){
    if(first->next!=0){
        deleteStrin_gListInGC(first->next);
    }
    delete first;
}



class GraphicCard{
public:
    GraphicCard(std::string="parallel.cl",mysint =0);
    ~GraphicCard();
    mysint findAddKernel(std::string,mysint=0);
    //Arguments: kernelName, preventAddition
    //If the kernel with given name exists, the
    //function will return the kernel's number.
    //Otherwise, the function will return -1.
    //If the second parameter is set to 0, the function will add the kernel with the given name.
    mysint deviceMemoryAccess(std::string,myint *, myint, mysint=0, myint=0);
    //The arguments are: std::string memBlockName, myint * memorySequence,
    //                      myint sLength, mysint action, myint writingShift=0
    //First the function will determine the index of the element of memObjectsGC that
    //corresponds to the memBlockName.
    //If it does not exist then the index will be announced to be -1.
    //action==0 - Locating the memory block named memBlockName, or creating a new one, if it
    //            does not exist. The content will be copied or overwritten from memorySequence.
    //            The first writingShift members of memoryBlock will be skipped.
    //action==1 - Reading the memory block of length sLength (starting from writingShift)
    //              and storing it into the sequence
    //action==2 - Updating only the length of the memory block. The new length is provided in sLength.
    //action==3 - Just determine the id number of the memory block with name memBlockName.
    
    
    mysint deviceMemoryAccess(std::string,cl_double *, myint, mysint=0, myint=0);
    // Same as before but accepts sequence of real numbers of type double
    
    mysint writeDeviceMemory(std::string, myint*, myint);
    mysint readDeviceMemory(std::string, myint*, myint);
    mysint writeDeviceMemory(std::string, double*, myint);
    mysint readDeviceMemory(std::string, double*, myint);
    
    mysint setKernelArguments(std::string,std::string*,mysint);
    // The arguments are: kernelName, listOfParameters, numberOfParametersInTheList
    mysint executeKernel(std::string,myint);
    

    mysint generateRandomNumbers(myint,myint,mysint=0);
    //Runs random number generator (in the first run, it will initialize them as well).
    //Arguments:    N - Size of each axis in the seed
    //              r - length of the numbers in binary
    //              indicator
    // - If indicator is 1, then the permutations will be created again.
    // - Note that permutations will be created again if N and/or r is changed from what it were before
    
    
protected:
    cl_context contextGC;
    cl_command_queue commandQueueGC;
    cl_program programGC;
    cl_device_id deviceGC;
    cl_kernel* kernflsGC;
    cl_mem *memObjectsGC;
    mysint numMemObjectsGC;
    mysint numberOfKernelsGC;
    std::string * kernelNamesGC;
    std::string * memObjNamesGC;
    size_t* preferred_workgroup_sizeGC;
    
    mysint *variablesCorrectlySetInKernelGC;
    
    //For random number generators
    myint *axisSizeGC;
    myint *xAxisGC;
    myint *yAxisGC;
    myint *lengthInBinaryGC;
    myint *balancedNumbersGC;
    myint *numBalancedNumbersGC;
    myint *shufflingPrimeGC;
    myint *xPermutationsGC;
    myint *yPermutationsGC;
    myint *pascalTriangleGC;
    myint *sizePascalTrGC;
    
    myint indFirstRNInitGC;
    myint currentSeedGC;
    myint uniformLimitGC;
    
    
    
    myint *randNumbersGC;
    myint *powersOfTwoGC;
    
    std::string kernelFileGC;
    std::string kernelPreLoadedGC;
    
    
    
    
    
    cl_context CreateContext(mysint =0);
    cl_command_queue CreateCommandQueue(cl_context, cl_device_id*);
    cl_program CreateProgram(cl_context, cl_device_id , const char* );
    myint treatError(cl_int, cl_context, cl_command_queue,
                     cl_program, cl_kernel *, myint,cl_mem *, myint);
    
    
    void Cleanup();
    mysint generatePermutationsGC(myint*,myint,myint);
    //myint binomGC(myint**, myint, myint);
    myint createBalancedNumbersGC(myint *,  myint *, myint* , myint *, myint*);
    myint createPascalTriangleGC(myint);
};



GraphicCard::GraphicCard(std::string programName,mysint forceCPU){
    contextGC=nullptr;
    commandQueueGC=nullptr;
    programGC=nullptr;
    deviceGC=0;
    kernflsGC=nullptr;
    numMemObjectsGC=0;
    memObjectsGC=nullptr;
    memObjNamesGC=nullptr;
    numberOfKernelsGC=0;
    kernelNamesGC=nullptr;
    preferred_workgroup_sizeGC=new size_t;
    variablesCorrectlySetInKernelGC=nullptr;
    
    kernelPreLoadedGC="gc_parallel_rnum.cl";
    axisSizeGC=new myint;
    *axisSizeGC=0;
    xAxisGC=nullptr;
    yAxisGC=nullptr;
    randNumbersGC=nullptr;
    lengthInBinaryGC=new myint;
    *lengthInBinaryGC=0;
    shufflingPrimeGC=new myint;
    *shufflingPrimeGC=3;// figure out some better algorithm for determining this.
    
    numBalancedNumbersGC=new myint;
    *numBalancedNumbersGC=0;
    
    xPermutationsGC=nullptr;
    yPermutationsGC=nullptr;
    indFirstRNInitGC=0;
    
    
    currentSeedGC=std::chrono::high_resolution_clock::now().time_since_epoch().count();
    
    powersOfTwoGC=nullptr;
    balancedNumbersGC=nullptr;
    contextGC=CreateContext(forceCPU);
    kernelFileGC=programName;
    
    pascalTriangleGC=nullptr;
    sizePascalTrGC=new myint;
    *sizePascalTrGC=0;
    
    
    
    mysint shouldQuit=0;
    
    if (contextGC == nullptr)
    {
        std::cerr << "Failed to create OpenCL context." << std::endl;
        shouldQuit=1;
    }
    if(shouldQuit==0){
    commandQueueGC = CreateCommandQueue(contextGC, &deviceGC);
    }
    if((shouldQuit==0)&& (commandQueueGC == nullptr))
    {
        Cleanup();
        shouldQuit= 1;
    }
    if(shouldQuit==0){
        const char * cChar = kernelFileGC.c_str();
        
        
        programGC = CreateProgram(contextGC, deviceGC, cChar);
    }
    if((shouldQuit==0)&&(programGC == nullptr))
    {
        Cleanup();
        shouldQuit=1;
    }

    
    
}
mysint GraphicCard::findAddKernel(std::string nameToAdd, mysint preventAddition){
    mysint i=0;
    mysint foundName=-1;
    while((i<numberOfKernelsGC)&&(foundName==-1)){
        if(kernelNamesGC[i]==nameToAdd){
            foundName=i;
        }
        i++;
    }
    if((foundName==-1)&&(preventAddition==0)){
        std::string *newKernelNames;
        cl_kernel *newKernfls;
        mysint *newVarCS;
        newKernelNames=new std::string[numberOfKernelsGC+1];
        newKernfls=new cl_kernel[numberOfKernelsGC+1];
        newVarCS=new mysint[numberOfKernelsGC+1];
        for(i=0;i<numberOfKernelsGC;i++){
            newKernelNames[i]=kernelNamesGC[i];
            newKernfls[i]=kernflsGC[i];
            newVarCS[i]=variablesCorrectlySetInKernelGC[i];
            
        }
        if(numberOfKernelsGC>0){
            delete[] kernelNamesGC;
            delete[] kernflsGC;
            delete[] variablesCorrectlySetInKernelGC;
        }
        kernelNamesGC=newKernelNames;
        kernflsGC=newKernfls;
        variablesCorrectlySetInKernelGC=newVarCS;
        
        kernelNamesGC[numberOfKernelsGC]=nameToAdd;
        const char * cChar = nameToAdd.c_str();
        cl_int error_ret;
        //std::cout<<"Kernel addition started: "<<std::endl;
        kernflsGC[numberOfKernelsGC]=clCreateKernel(programGC,cChar,&error_ret);
        myint jErr=0;
        while((jErr<1000)&&(error_ret!=0)){
            //std::cout<<"Creating kernel attempt: "<<jErr+1<<" name: "<< cChar<<std::endl;
            kernflsGC[numberOfKernelsGC]=clCreateKernel(programGC,cChar,&error_ret);
            
            jErr++;
        }
        if(numberOfKernelsGC<3){
            clGetKernelWorkGroupInfo (kernflsGC[numberOfKernelsGC],
                                      deviceGC,
                                      CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
                                      sizeof(size_t),
                                      preferred_workgroup_sizeGC,
                                      nullptr);
        }
        variablesCorrectlySetInKernelGC[numberOfKernelsGC]=0;
        numberOfKernelsGC++;
    }
    return foundName;
}









cl_context GraphicCard::CreateContext(mysint forceCPU)
{
    cl_int errNum;
    cl_uint numPlatforms;
    cl_platform_id firstPlatformId;
    cl_context context = nullptr;
    
    
    errNum = clGetPlatformIDs(1, &firstPlatformId, &numPlatforms);
    if (errNum != CL_SUCCESS || numPlatforms <= 0)
    {
        std::cerr << "No OpenCL platforms." << std::endl;
        return nullptr;
    }
    
    cl_context_properties contextProperties[] =
    {
        CL_CONTEXT_PLATFORM,
        (cl_context_properties)firstPlatformId,
        0
    };
    
    //Select GPU:
    
    
    if(forceCPU==0){
        context = clCreateContextFromType(contextProperties, CL_DEVICE_TYPE_GPU,nullptr, nullptr, &errNum);
        if (errNum != CL_SUCCESS)
        {
            std::cout << "No GPU context. Trying CPU." << std::endl;
            context = clCreateContextFromType(contextProperties, CL_DEVICE_TYPE_CPU,
                                              nullptr, nullptr, &errNum);
            if (errNum != CL_SUCCESS)
            {
                std::cerr << "Failed to create an OpenCL GPU or CPU context." << std::endl;
                return nullptr;
            }
        }
    }
    //Comparison: Select CPU:
    
    if(forceCPU==1){
        context = clCreateContextFromType(contextProperties, CL_DEVICE_TYPE_CPU,
                                          nullptr, nullptr, &errNum);
        if (errNum != CL_SUCCESS)
        {
            std::cerr << "Failed to create an OpenCL GPU or CPU context." << std::endl;
            return nullptr;
        }
    }
    
    
    
    return context;
}


cl_command_queue GraphicCard::CreateCommandQueue(cl_context context, cl_device_id *device)
{
    cl_int errNum;
    cl_device_id *devices;
    cl_command_queue commandQueue = nullptr;
    size_t deviceBufferSize = -1;
    
    // First get the size of the devices buffer
    errNum = clGetContextInfo(context, CL_CONTEXT_DEVICES, 0, nullptr, &deviceBufferSize);
    
    if (errNum != CL_SUCCESS)
    {
        std::cerr << "Failed call to clGetContextInfo(...,GL_CONTEXT_DEVICES,...)";
        return nullptr;
    }
    
    if (deviceBufferSize <= 0)
    {
        std::cerr << "No devices available.";
        return nullptr;
    }
    
    // Allocate memory for the devices buffer
    myint numDevices=(myint)(deviceBufferSize / sizeof(cl_device_id));
    devices = new cl_device_id[deviceBufferSize / sizeof(cl_device_id)];
    
    //std::cout<<"Number of ids "<<deviceBufferSize / sizeof(cl_device_id)<<std::endl;
    errNum = clGetContextInfo(context, CL_CONTEXT_DEVICES, deviceBufferSize, devices, nullptr);
    if (errNum != CL_SUCCESS)
    {
        delete [] devices;
        std::cerr << "Failed to get device IDs";
        return nullptr;
    }
    
    // If I am at computer at work then deviceToChoose should be 0; if I am at home it should be 1;
    
    myint deviceToChoose=numDevices-1;
    // deviceToChoose=0;
    //std::cout<<"Selecting device "<<deviceToChoose<<std::endl;
    
    commandQueue = clCreateCommandQueue(context, devices[deviceToChoose], 0, nullptr);
    //std::cout<<"Device ID "<<devices[deviceToChoose]<<std::endl;
    if (commandQueue == nullptr)
    {
        delete [] devices;
        std::cerr << "Failed to create commandQueue for device 0";
        return nullptr;
    }
    
    *device = devices[deviceToChoose];
    delete [] devices;
    return commandQueue;
}

cl_program GraphicCard::CreateProgram(cl_context context, cl_device_id device, const char* fileName)
{
    cl_int errNum;
    cl_program program;
    const char * fileName0 = kernelPreLoadedGC.c_str();
    std::ifstream kernelFile0(fileName0, std::ios::in);
    if (!kernelFile0.is_open())
    {
        std::cerr << "Failed to open file for reading: " << fileName0 << std::endl;
        return nullptr;
    }
    
    std::ostringstream oss0;
    oss0 << kernelFile0.rdbuf();
    
    std::string srcStdStr0 = oss0.str();
    
    std::ifstream kernelFile(fileName, std::ios::in);
    if (!kernelFile.is_open())
    {
        std::cerr << "Failed to open file for reading: " << fileName << std::endl;
        return nullptr;
    }
    
    std::ostringstream oss;
    oss << kernelFile.rdbuf();
    
    std::string srcStdStr = srcStdStr0;
    srcStdStr+=oss.str();
    
    const char *srcStr = srcStdStr.c_str();
   // std::cout<<srcStr<<std::endl;
    program = clCreateProgramWithSource(context, 1,
                                        (const char**)&srcStr,
                                        nullptr, nullptr);
    if (program == nullptr)
    {
        std::cerr << "Failed to create CL program from source." << std::endl;
        return nullptr;
    }
    
    errNum = clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
    if (errNum != CL_SUCCESS)
    {
        // Determine the reason for the error
        char buildLog[16384];
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG,
                              sizeof(buildLog), buildLog, nullptr);
        
        std::cerr << "Error in kernel: " << std::endl;
        std::cerr << buildLog;
        clReleaseProgram(program);
        return nullptr;
    }
    
    return program;
}



void GraphicCard::Cleanup()
{
    for (myint i = 0; i < numMemObjectsGC; i++)
    {
        if (memObjectsGC[i] != nullptr) 
            clReleaseMemObject(memObjectsGC[i]);
    }
    if (commandQueueGC != nullptr) 
        clReleaseCommandQueue(commandQueueGC);
    
    for(myint k=0;k<numberOfKernelsGC;k++){
        if(kernflsGC[k]!=nullptr){ 
            clReleaseKernel(kernflsGC[k]);
        }
    }
    
    
    
    
    if (programGC != nullptr) 
        clReleaseProgram(programGC);
    
    if (contextGC != nullptr) 
        clReleaseContext(contextGC);
    
    delete[] kernelNamesGC;
    delete[] kernflsGC;
    delete[] memObjNamesGC;
    delete[] memObjectsGC;
    delete preferred_workgroup_sizeGC;
    delete axisSizeGC;
    delete lengthInBinaryGC;
    delete shufflingPrimeGC;
    delete[] variablesCorrectlySetInKernelGC;
    delete[] powersOfTwoGC;
    
 
    delete[] xAxisGC;

    delete[] yAxisGC;

    delete[] xPermutationsGC;

    delete[] yPermutationsGC;

    delete[] balancedNumbersGC;

    delete numBalancedNumbersGC;

    delete[] pascalTriangleGC;

    delete sizePascalTrGC;




    delete randNumbersGC;
     
}



mysint GraphicCard::deviceMemoryAccess(std::string memBlockName,myint *memorySequence, myint sLength, mysint action, myint writingShift){
    mysint i=0;
    mysint foundName=-1;
    cl_int errNum;
    cl_event eventE;
    myint *pSLength;
    pSLength=new myint;
    *pSLength=sLength;
    while((i<numMemObjectsGC)&&(foundName==-1)){
        if(memObjNamesGC[i]==memBlockName){
            foundName=i;
        }
        i++;
    }
    if((foundName==-1)&&(action==0)){
        std::string *newMemObjNames;
        cl_mem *newMemObjects;
        newMemObjNames=new std::string[numMemObjectsGC+1];
        newMemObjects=new cl_mem[2*(numMemObjectsGC+1)];
        for(i=0;i<numMemObjectsGC;i++){
            newMemObjNames[i]=memObjNamesGC[i];
            newMemObjects[2*i]=memObjectsGC[2*i];
            newMemObjects[2*i+1]=memObjectsGC[2*i+1];
        }
        if(numMemObjectsGC>0){
            delete[] memObjNamesGC;
            delete[] memObjectsGC;
        }
        

        memObjNamesGC=newMemObjNames;
        memObjectsGC=newMemObjects;
        memObjNamesGC[numMemObjectsGC]=memBlockName;
        const char * cChar = memBlockName.c_str();
        
        memObjectsGC[2*numMemObjectsGC] = clCreateBuffer(contextGC, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                       sizeof(myint) , pSLength, nullptr);
        memObjectsGC[2*numMemObjectsGC+1] = clCreateBuffer(contextGC, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                       sizeof(myint) * sLength , memorySequence, nullptr);
        
        numMemObjectsGC++;

    }
    if((foundName!=-1)&&(action==1)){
        errNum = clEnqueueReadBuffer(commandQueueGC, memObjectsGC[2*foundName+1], CL_TRUE,
                                     sizeof(myint) * writingShift,   sizeof(myint)* sLength, memorySequence,
                                     0, nullptr, &eventE);
        clWaitForEvents(1, &eventE);
    }
    
    if((foundName!=-1)&&(action==0)){
        clEnqueueWriteBuffer(commandQueueGC, memObjectsGC[2*foundName+1], CL_TRUE, sizeof(myint) * writingShift,   sizeof(myint)* sLength,memorySequence,0,nullptr, &eventE);
        clWaitForEvents(1, &eventE);
    }

    if((foundName!=-1)&&(action==2)){
        clEnqueueWriteBuffer(commandQueueGC, memObjectsGC[2*foundName], CL_TRUE, 0,     sizeof(myint),pSLength,0,nullptr, &eventE);
        clWaitForEvents(1, &eventE);
    }
    delete pSLength;
    return foundName;
}





mysint GraphicCard::deviceMemoryAccess(std::string memBlockName,cl_double *memorySequence, myint sLength, mysint action, myint writingShift){
    mysint i=0;
    mysint foundName=-1;
    cl_int errNum;
    cl_event eventE;
    myint *pSLength;
    pSLength=new myint;
    *pSLength=sLength;
    while((i<numMemObjectsGC)&&(foundName==-1)){
        if(memObjNamesGC[i]==memBlockName){
            foundName=i;
        }
        i++;
    }
    if((foundName==-1)&&(action==0)){
        std::string *newMemObjNames;
        cl_mem *newMemObjects;
        newMemObjNames=new std::string[numMemObjectsGC+1];
        newMemObjects=new cl_mem[2*(numMemObjectsGC+1)];
        for(i=0;i<numMemObjectsGC;i++){
            newMemObjNames[i]=memObjNamesGC[i];
            newMemObjects[2*i]=memObjectsGC[2*i];
            newMemObjects[2*i+1]=memObjectsGC[2*i+1];
        }
        if(numMemObjectsGC>0){
            delete[] memObjNamesGC;
            delete[] memObjectsGC;
        }
        
        
        memObjNamesGC=newMemObjNames;
        memObjectsGC=newMemObjects;
        memObjNamesGC[numMemObjectsGC]=memBlockName;
        const char * cChar = memBlockName.c_str();
        
        memObjectsGC[2*numMemObjectsGC] = clCreateBuffer(contextGC, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                                         sizeof(myint) , pSLength, nullptr);
        memObjectsGC[2*numMemObjectsGC+1] = clCreateBuffer(contextGC, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                                           sizeof(cl_double) * sLength , memorySequence, nullptr);
        
        numMemObjectsGC++;
        
    }
    if((foundName!=-1)&&(action==1)){
        errNum = clEnqueueReadBuffer(commandQueueGC, memObjectsGC[2*foundName+1], CL_TRUE,
                                     sizeof(cl_double) * writingShift,   sizeof(cl_double)* sLength, memorySequence,
                                     0, nullptr, &eventE);
        clWaitForEvents(1, &eventE);
    }
    
    if((foundName!=-1)&&(action==0)){
        clEnqueueWriteBuffer(commandQueueGC, memObjectsGC[2*foundName+1], CL_TRUE, sizeof(cl_double) * writingShift,   sizeof(cl_double)* sLength,memorySequence,0,nullptr, &eventE);
        clWaitForEvents(1, &eventE);
    }
    
    if((foundName!=-1)&&(action==2)){
        clEnqueueWriteBuffer(commandQueueGC, memObjectsGC[2*foundName], CL_TRUE, 0,     sizeof(myint),pSLength,0,nullptr, &eventE);
        clWaitForEvents(1, &eventE);
    }
    delete pSLength;
    return foundName;
}




mysint GraphicCard::setKernelArguments(std::string kernelName,std::string* parNames,mysint numberOfParameters){
    mysint kernelNumber=findAddKernel(kernelName);
    mysint forReturn=0;
    myint *memSequence=0;
    mysint helpNum;
    cl_int errNum;

    if(kernelNumber==-1){
        kernelNumber=findAddKernel(kernelName);
    }
    
    mysint *parNumbers;
    parNumbers=new mysint[numberOfParameters];
    for(mysint i=0;i<numberOfParameters;i++){
        parNumbers[i]=deviceMemoryAccess(parNames[i],memSequence , 0, 3, 0);
        if(parNumbers[i]==-1){
            forReturn=-1;
        }
    }
    if(forReturn!=-1){
        for(mysint i=0;i<numberOfParameters;i++){
            helpNum=2*parNumbers[i];
            helpNum++;
            errNum=clSetKernelArg(kernflsGC[kernelNumber], i, sizeof(cl_mem), &memObjectsGC[helpNum]);
            if (errNum != CL_SUCCESS)
            {
                std::cerr << "Error setting kernel argument " << i<<" "<<errNum<<std::endl;
                std::cerr << "Kernel number: " << kernelNumber<<", kernel code:"<<kernflsGC[kernelNumber]<<std::endl;
                
                if(errNum==CL_INVALID_KERNEL){
                    std::cout<<"CL_INVALID_KERNEL"<<std::endl;
                }
                if(errNum==CL_INVALID_ARG_INDEX){
                    std::cout<<"CL_INVALID_ARG_INDEX"<<std::endl;
                }
                if(errNum==CL_INVALID_ARG_VALUE){
                    std::cout<<"CL_INVALID_ARG_VALUE"<<std::endl;
                }
                if(errNum==CL_INVALID_MEM_OBJECT){
                    std::cout<<"CL_INVALID_MEM_OBJECT"<<std::endl;
                }
                if(errNum==CL_INVALID_SAMPLER){
                    std::cout<<"CL_INVALID_SAMPLER"<<std::endl;
                }
                if(errNum==CL_INVALID_ARG_SIZE){
                    std::cout<<"CL_INVALID_ARG_SIZE"<<std::endl;
                }
                if(errNum==CL_OUT_OF_RESOURCES){
                    std::cout<<"CL_OUT_OF_RESOURCES"<<std::endl;
                }
                if(errNum==CL_OUT_OF_HOST_MEMORY){
                    std::cout<<"CL_OUT_OF_HOST_MEMORY"<<std::endl;
                }
                forReturn=-2;
            }
            
            if(forReturn!=-2){
                variablesCorrectlySetInKernelGC[kernelNumber]=1;
            }
            
        }
    }
    
    delete[] parNumbers;
    return forReturn;
    
}






mysint GraphicCard::executeKernel(std::string kernelName, myint numberOfProcessingElements){
    mysint kernelNumber=findAddKernel(kernelName, 1);
    cl_int errNum;
    cl_event eventE;
    size_t *globalWorkSizeP;
    size_t *localWorkSizeP;
    globalWorkSizeP=new size_t;
    localWorkSizeP=new size_t;
    size_t pws;
    
    if(kernelNumber==-1){
        kernelNumber=findAddKernel(kernelName);
        if(kernelNumber==-1){
            kernelNumber=findAddKernel(kernelName,1);
        }
    }
    
    if(kernelNumber!=-1){
        
        if(variablesCorrectlySetInKernelGC[kernelNumber]==0){
          //  The variables were not set correctly
          //  We need to set them up.

            
            std::ifstream t(kernelFileGC);
            std::string str;
            
            t.seekg(0, std::ios::end);
            str.reserve(t.tellg());
            t.seekg(0, std::ios::beg);
            
            str.assign((std::istreambuf_iterator<char>(t)),
                       std::istreambuf_iterator<char>());
 
            replace (str.begin(), str.end(), '*' , ' ');
            replace (str.begin(), str.end(), ',' , ' ');
            replace (str.begin(), str.end(), ';' , ' ');
            replace (str.begin(), str.end(), ')' , ' ');
            replace (str.begin(), str.end(), '&' , ' ');
            replace (str.begin(), str.end(), '(' , ' ');
            
            size_t index;
            while ((index = str.find("const")) != std::string::npos)
                str.replace(index, 5, " ");
            
            while ((index = str.find("void")) != std::string::npos)
                str.replace(index, 4, " ");
            while ((index = str.find("__local")) != std::string::npos)
                str.replace(index, 7, " __global ");
            
            
            std::regex word_regex("(\\S+)");
            std::sregex_iterator words_begin =
            std::sregex_iterator(str.begin(), str.end(), word_regex);
            std::sregex_iterator  words_end = std::sregex_iterator();

            std::sregex_iterator cWord;
            std::string wordNext;
            mysint ik=0;
            mysint jk=0;
            cWord=words_begin;
            
            Strin_gList* firstListEl, *currentListEl;
            firstListEl=new Strin_gList;
            currentListEl=firstListEl;
            (*firstListEl).next=0;
            while(cWord!=words_end){
                wordNext=(*cWord).str();
                if(wordNext=="__kernel"){
                    ++cWord;
                    if(cWord!=words_end){
                        wordNext=(*cWord).str();
                        if(wordNext==kernelName){
                            ++cWord;
                            while(cWord!=words_end){
                                wordNext=(*cWord).str();
                                if(wordNext=="__global"){
                                    ++cWord;
                                    if(cWord!=words_end){
                                        ++cWord;
                                        if(cWord!=words_end){
                                            ik++;
                                            (*currentListEl).content=(*cWord).str();
                                            (*currentListEl).next=new Strin_gList;
                                            currentListEl=(*currentListEl).next;
                                            (*currentListEl).next=0;
                                        }
                                    }
                                }
                                if(cWord!=words_end){++cWord;}
                                if(cWord!=words_end){
                                    if(((*cWord).str())=="__kernel"){
                                        cWord=words_end;
                                    }
                                }

                            }
                        }
                    }
                }
                if(cWord!=words_end){++cWord;}
            }
            std::string* myStList;
            myStList=new std::string[ik];
            currentListEl=firstListEl;
            for(myint i=0;i<ik;i++){
                myStList[i]=currentListEl->content;
                currentListEl=currentListEl->next;
            }
 
            
            deleteStrin_gListInGC(firstListEl);
            
            setKernelArguments(kernelName,myStList,ik);

            delete[] myStList;
 
        }

        pws= (size_t) numberOfProcessingElements;
        pws= ( (pws/ (*preferred_workgroup_sizeGC))+1) * (*preferred_workgroup_sizeGC);
        *globalWorkSizeP=pws;
        *localWorkSizeP=*preferred_workgroup_sizeGC;
        
        
        errNum=clEnqueueNDRangeKernel(commandQueueGC, kernflsGC[kernelNumber], 1, nullptr,
                                      globalWorkSizeP, localWorkSizeP,
                                      0, nullptr, &eventE);
        
        if (errNum != CL_SUCCESS){
            treatError(errNum,contextGC, commandQueueGC, programGC, kernflsGC, numberOfKernelsGC, memObjectsGC, numMemObjectsGC);
            return 1;
        }
        clWaitForEvents(1, &eventE);
       
    }
    delete globalWorkSizeP;
    delete localWorkSizeP;
    
    
    return kernelNumber;
    
}


myint GraphicCard::treatError(cl_int errNum, cl_context context, cl_command_queue commandQueue,
                              cl_program program, cl_kernel *kernfls, myint numKernels,cl_mem *memObjects, myint numMemObjects){
    std::cerr << "Error queuing kernel for execution. Kernel 1! " << " "<<errNum<< std::endl;
    if(errNum==CL_INVALID_PROGRAM_EXECUTABLE){std::cout<<"CL_INVALID_PROGRAM_EXECUTABLE"<<std::endl;}
    if(errNum==CL_INVALID_COMMAND_QUEUE){std::cout<<"CL_INVALID_COMMAND_QUEUE"<<std::endl;}
    if(errNum==CL_INVALID_KERNEL){std::cout<<"CL_INVALID_KERNEL"<<std::endl;}
    if(errNum==CL_INVALID_CONTEXT){std::cout<<"CL_INVALID_CONTEXT"<<std::endl;}
    if(errNum==CL_INVALID_KERNEL_ARGS){std::cout<<"CL_INVALID_KERNEL_ARGS"<<std::endl;}
    if(errNum==CL_INVALID_WORK_DIMENSION){std::cout<<"CL_INVALID_WORK_DIMENSION"<<std::endl;}
    if(errNum==CL_INVALID_GLOBAL_WORK_SIZE){std::cout<<"CL_INVALID_GLOBAL_WORK_SIZE"<<std::endl;}
    if(errNum==CL_INVALID_GLOBAL_OFFSET){std::cout<<"CL_INVALID_GLOBAL_OFFSET"<<std::endl;}
    if(errNum==CL_INVALID_WORK_GROUP_SIZE){std::cout<<"CL_INVALID_WORK_GROUP_SIZE"<<std::endl;}
    if(errNum==CL_INVALID_WORK_ITEM_SIZE){std::cout<<"CL_INVALID_WORK_ITEM_SIZE"<<std::endl;}
    if(errNum==CL_MISALIGNED_SUB_BUFFER_OFFSET){std::cout<<"CL_MISALIGNED_SUB_BUFFER_OFFSET"<<std::endl;}
    if(errNum==CL_OUT_OF_RESOURCES){std::cout<<"CL_OUT_OF_RESOURCES"<<std::endl;}
    if(errNum==CL_MEM_OBJECT_ALLOCATION_FAILURE){std::cout<<"CL_MEM_OBJECT_ALLOCATION_FAILURE"<<std::endl;}
    if(errNum==CL_INVALID_EVENT_WAIT_LIST){std::cout<<"CL_INVALID_EVENT_WAIT_LIST"<<std::endl;}
    if(errNum==CL_OUT_OF_HOST_MEMORY){std::cout<<"CL_OUT_OF_HOST_MEMORY"<<std::endl;}
    
    
    Cleanup();
    return 1;
}

mysint GraphicCard::writeDeviceMemory(std::string memBlockName,myint *memorySequence, myint sLength){
    return deviceMemoryAccess(memBlockName,memorySequence,sLength);
}
mysint GraphicCard::readDeviceMemory(std::string memBlockName,myint *memorySequence, myint sLength){
    return deviceMemoryAccess(memBlockName,memorySequence,sLength,1);
}
mysint GraphicCard::writeDeviceMemory(std::string memBlockName,cl_double *memorySequence, myint sLength){
    return deviceMemoryAccess(memBlockName,memorySequence,sLength);
}
mysint GraphicCard::readDeviceMemory(std::string memBlockName,cl_double *memorySequence, myint sLength){
    return deviceMemoryAccess(memBlockName,memorySequence,sLength,1);
}
GraphicCard::~GraphicCard(){
    Cleanup();
    
}











