



#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <map>
#include <set>
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







typedef long int mysint;
typedef long int myint;




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

    
    //    mysint deviceMemoryAccess(std::string,myint *, myint, mysint=0, myint=0);
    
    
    template <typename int_doub> mysint deviceMemoryAccess(std::string,int_doub *, myint, mysint=0, myint=0);
    
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
    //action==4 - Delete memory block;
    
   // mysint deviceMemoryAccess(std::string,cl_double *, myint, mysint=0, myint=0);
    // Same as before but accepts sequence of real numbers of type double
    
    
    
    
    
    template <typename int_doub> mysint writeDeviceMemory(std::string, int_doub *, myint);
    template <typename int_doub> mysint readDeviceMemory(std::string, int_doub *, myint);
    //mysint writeDeviceMemory(std::string, double*, myint);
   // mysint readDeviceMemory(std::string, double*, myint);
    mysint freeDeviceMemory(std::string);
    
    mysint setKernelArguments(std::string,std::string*,mysint);
    // The arguments are: kernelName, listOfParameters, numberOfParametersInTheList
    mysint executeKernel(std::string,myint);
    
    
    mysint generateRandomNumbers(myint,myint,mysint=0);
    //Runs random number generator (in the first run, it will initialize them as well).
    //Arguments:    N - Size of each axis in the seed
    //              r - length of the numbers in binary
    //              indicator
    // - If indicator is 1, then the permutations will be created again.
    // - Note that permutations will be created the first time the generator is executed.
   
  
    mysint generateNormalBeasleySpringerMoro(myint, myint, myint, myint*,
                                             myint*,myint,//these two are for debugging
                                             cl_double =0.0, cl_double =1.0, mysint=0);
    
    //Generates normal distribution
    //Arguments:    N - Size of each axis in the seed
    //              r - length of the numbers in binary
    //              precisionRequest - the number of integers from the discrete distribution that will be used
    //                                  to generate uniform [0,1] distribution.
    //              sampleLength - (myint*) contains the size of the generated sample
    //                              sampleLength = N^2 /  precisionRequest
    //              parameter1 - used for normal and denotes the mean
    //              parameter2 - for normal it denotes sigma, for exponential, it denotes lambda.
    //              indicator - If the indicator is 1, then the permutations will be created again.
    
    
    mysint generateExponential(myint, myint, myint, myint*,
                                             myint*,myint,//these two are for debugging
                                             cl_double =1.0, mysint=0);
    
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
    
    std::set<myint>* kernelResponsibilitiesGC;
    
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
    myint *highestPrecisionGC;
    myint *precisionRequestGC;
    
    myint *exponentKGC;
    myint *sizeForRejectionSamplingGC;
    myint *sampleLengthGC;
    myint *inspectorRejSampGC;
    
    
    cl_double *parameter1GC;
    cl_double *parameter2GC;
    
    cl_double *randomSampleGC;

    
    
    myint indFirstRNInitGC;
    
    myint inspectorExecutedGC;
    myint sampleGenKernExecutedGC;
    
    myint currentSeedGC;
    myint uniformLimitGC;
    
    myint lastRandAlgUsedGC; //each random distribution has a code. Normal is 7; exponential is 8
                             //when the same algorithm with same sample sizes is used repeatedly, some
                             //savings in time are possible.
                             //This variable is updated to contain the label of the last algorithm used.
    
    
    myint *randNumbersGC;
    myint *powersOfTwoGC;
    
    myint numberOfDistributionsGC;
    std::string* distributionKernelNamesGC;
    std::string normalBSMNameGC;
    std::string exponentialDistNameGC;
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
    *preferred_workgroup_sizeGC=128;
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
    inspectorExecutedGC=0;
    sampleGenKernExecutedGC=0;
    uniformLimitGC=10000;
    currentSeedGC=std::chrono::high_resolution_clock::now().time_since_epoch().count();
    
    powersOfTwoGC=nullptr;
    balancedNumbersGC=nullptr;
    contextGC=CreateContext(forceCPU);
    kernelFileGC=programName;
    
    pascalTriangleGC=nullptr;
    sizePascalTrGC=new myint;
    *sizePascalTrGC=0;
    
    highestPrecisionGC=new myint;
    *highestPrecisionGC=100000;
    
    precisionRequestGC=new myint;
    *precisionRequestGC=0;
    exponentKGC=new myint;
    *exponentKGC=3;
    sizeForRejectionSamplingGC=new myint;
    *sizeForRejectionSamplingGC=5;
    
    sampleLengthGC=new myint;
    *sampleLengthGC=0;
    
    parameter1GC=new cl_double;
    *parameter1GC=0;
    parameter2GC=new cl_double;
    *parameter2GC=1;
    randomSampleGC=nullptr;
    
    inspectorRejSampGC=nullptr;
    kernelResponsibilitiesGC=nullptr;

    
    normalBSMNameGC="normalBSMGC";
    exponentialDistNameGC="exponentialDistGC";
    lastRandAlgUsedGC=-1;
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
    delete[] kernelResponsibilitiesGC;
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
    
    delete highestPrecisionGC;
    delete precisionRequestGC;
    
    delete exponentKGC;
    delete sizeForRejectionSamplingGC;
    delete sampleLengthGC;
    
    delete parameter1GC;
    delete parameter2GC;
    
    delete[] randomSampleGC;
    
    delete[] inspectorRejSampGC;


    
    delete randNumbersGC;
    
}


template <typename int_doub>
mysint GraphicCard::deviceMemoryAccess(std::string memBlockName,int_doub *memorySequence,
                                       myint sLength, mysint action, myint writingShift){
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
        std::set<myint>*newKernelResp;
        newMemObjNames=new std::string[numMemObjectsGC+1];
        newKernelResp=new std::set<myint>[numMemObjectsGC+1];
        newMemObjects=new cl_mem[2*(numMemObjectsGC+1)];
        for(i=0;i<numMemObjectsGC;i++){
            newMemObjNames[i]=memObjNamesGC[i];
            newKernelResp[i]=kernelResponsibilitiesGC[i];
            newMemObjects[2*i]=memObjectsGC[2*i];
            newMemObjects[2*i+1]=memObjectsGC[2*i+1];
        }
        if(numMemObjectsGC>0){
            delete[] memObjNamesGC;
            delete[] memObjectsGC;
            delete[] kernelResponsibilitiesGC;
        }
        
        
        memObjNamesGC=newMemObjNames;
        memObjectsGC=newMemObjects;
        kernelResponsibilitiesGC=newKernelResp;
        memObjNamesGC[numMemObjectsGC]=memBlockName;
        //kernelResponsibilitiesGC[numMemObjectsGC] is empty set without us having to do anything
        const char * cChar = memBlockName.c_str();
        
        memObjectsGC[2*numMemObjectsGC] = clCreateBuffer(contextGC, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                                         sizeof(myint) , pSLength, nullptr);
        memObjectsGC[2*numMemObjectsGC+1] = clCreateBuffer(contextGC, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                                           sizeof(int_doub) * sLength , memorySequence, nullptr);
        
        numMemObjectsGC++;
        
    }
    if((foundName!=-1)&&(action==1)){
        errNum = clEnqueueReadBuffer(commandQueueGC, memObjectsGC[2*foundName+1], CL_TRUE,
                                     sizeof(int_doub) * writingShift,   sizeof(int_doub)* sLength, memorySequence,
                                     0, nullptr, &eventE);
        clWaitForEvents(1, &eventE);
    }
    
    if((foundName!=-1)&&(action==0)){
        clEnqueueWriteBuffer(commandQueueGC, memObjectsGC[2*foundName+1], CL_TRUE, sizeof(int_doub) * writingShift,   sizeof(int_doub)* sLength,memorySequence,0,nullptr, &eventE);
        clWaitForEvents(1, &eventE);
    }
    
    if((foundName!=-1)&&(action==2)){
        clEnqueueWriteBuffer(commandQueueGC, memObjectsGC[2*foundName], CL_TRUE, 0,     sizeof(myint),pSLength,0,nullptr, &eventE);
        clWaitForEvents(1, &eventE);
    }
    if((foundName!=-1)&&(action==4)){
        if (memObjectsGC[2*foundName] != nullptr){
            clReleaseMemObject(memObjectsGC[2*foundName]);
        }
        if (memObjectsGC[2*foundName+1] != nullptr){
            clReleaseMemObject(memObjectsGC[2*foundName+1]);
        }
        numMemObjectsGC--;
        for(myint i=foundName;i<numMemObjectsGC;i++){
            memObjNamesGC[i]=memObjNamesGC[i+1];
            memObjectsGC[2*i]=memObjectsGC[2*(i+1)];
            memObjectsGC[2*i+1]=memObjectsGC[2*(i+1)+1];
            kernelResponsibilitiesGC[i]=kernelResponsibilitiesGC[i+1];
        }
    }
    
    
    delete pSLength;
    return foundName;
}




/*
mysint GraphicCard::deviceMemoryAccess(std::string memBlockName,myint *memorySequence, myint sLength, mysint action, myint writingShift){
    return  deviceMemoryAccTemp(memBlockName, memorySequence,  sLength, action, writingShift);
}


mysint GraphicCard::deviceMemoryAccess(std::string memBlockName,cl_double *memorySequence, myint sLength, mysint action, myint writingShift){
    return   deviceMemoryAccTemp(memBlockName, memorySequence,  sLength, action, writingShift);
}


*/

mysint GraphicCard::setKernelArguments(std::string kernelName,std::string* parNames,mysint numberOfParameters){
    mysint kernelNumber=findAddKernel(kernelName);
    mysint forReturn=0;
    myint *memSequence=nullptr;
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
            kernelResponsibilitiesGC[parNumbers[i]].insert(kernelNumber);
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
            std::cout<<"Error in kernel "<<kernelNumber<<std::endl;
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
    std::cerr << "Error queuing kernel for execution! " << " "<<errNum<< std::endl;
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

template <typename int_doub> mysint GraphicCard::writeDeviceMemory(std::string memBlockName,int_doub *memorySequence, myint sLength){
    
    mysint memId=deviceMemoryAccess(memBlockName,memorySequence,sLength,3);
    if(memId!=-1){
        // We will check whether the new length is bigger than the allocated length.
        // If that is the case, we need to reallocate more space.
        myint currentLength;
        cl_int errNum;
        cl_event eventE;
        errNum = clEnqueueReadBuffer(commandQueueGC, memObjectsGC[2*memId], CL_TRUE,
                                     //sizeof(myint) * writingShift,   sizeof(myint)* sLength, memorySequence,
                                           0,   sizeof(myint) , &currentLength,
                                     0, nullptr, &eventE);
        clWaitForEvents(1, &eventE);
        if(currentLength<sLength){
            std::set<myint>::iterator it;
            for(it=(kernelResponsibilitiesGC[memId]).begin();it!=(kernelResponsibilitiesGC[memId]).end();++it){
                variablesCorrectlySetInKernelGC[*it]=0;
            }
            deviceMemoryAccess(memBlockName,memorySequence,sLength,4);
           
        }
        
    }
    
    
    return deviceMemoryAccess(memBlockName,memorySequence,sLength);
}
template <typename int_doub> mysint GraphicCard::readDeviceMemory(std::string memBlockName,int_doub *memorySequence, myint sLength){
    return deviceMemoryAccess(memBlockName,memorySequence,sLength,1);
}

/*
mysint GraphicCard::writeDeviceMemory(std::string memBlockName,cl_double *memorySequence, myint sLength){
    return deviceMemoryAccess(memBlockName,memorySequence,sLength);
}
mysint GraphicCard::readDeviceMemory(std::string memBlockName,cl_double *memorySequence, myint sLength){
    return deviceMemoryAccess(memBlockName,memorySequence,sLength,1);
}
*/
mysint GraphicCard::freeDeviceMemory(std::string memBlockName){
    myint*memorySequence, sLength;
    return deviceMemoryAccess(memBlockName,memorySequence,sLength,4);
}



GraphicCard::~GraphicCard(){
    Cleanup();
    
}













class TripleGC{
public:
    myint a;
    myint b;
    myint c ;
    myint sorting;
    TripleGC();
    TripleGC(myint, myint, myint);
    TripleGC(myint, myint, myint,myint);
    mysint operator<(const TripleGC&);
    mysint operator>(const TripleGC&);
    mysint operator=(const TripleGC&);
};
TripleGC::TripleGC(){
    a=0;
    b=0;
    c =0;
    sorting=-1;
}
TripleGC::TripleGC(myint p, myint q, myint r){
    a=p;
    b=q;
    c=r;
}
TripleGC::TripleGC(myint p, myint q, myint r,myint s ){
    a=p;
    b=q;
    c=r;
    sorting=s;
}
mysint TripleGC::operator=(const TripleGC& t){
    a=t.a;
    b=t.b;
    c=t.c;
    sorting=t.sorting;
    return 1;
}

mysint TripleGC::operator<(const TripleGC& t){
    if(a<t.a){
        return 1;
    }
    if(a>t.a){
        return 0;
    }
    if(b<t.b){
        return 1;
    }
    if(b>t.b){
        return 0;
    }
    if(c<t.c){
        return 1;
    }
    if(c>t.c){
        return 0;
    }
    return 0;
}

mysint TripleGC::operator>(const TripleGC& t){
    if(a>t.a){
        return 1;
    }
    if(a<t.a){
        return 0;
    }
    if(b>t.b){
        return 1;
    }
    if(b<t.b){
        return 0;
    }
    if(c>t.c){
        return 1;
    }
    if(c<t.c){
        return 0;
    }
    return 0;
}
myint mergeSortTGC(TripleGC *seq, myint lS,  TripleGC *sSeq=nullptr, myint lSS=0){
    //First job is to sort the first sequence seq
    // and the second job is to add the sorted sequence sSeq into seq
    // If the length of the first sequence is 1 or 0, there is no need for sorting.
    if(lS>1){
        //The sorting will be done by dividing the sequence in two parts.
        myint middle=(lS+1)/2;
        mergeSortTGC(seq+middle,lS-middle);
        mergeSortTGC(seq,middle,seq+middle,lS-middle);
    }
    if(lSS>0){
        myint length=lS+lSS;
        TripleGC* help;
        help=new TripleGC[length];
        myint r=0, rL=0,rR=0;
        while(r<length){
            if(rL>=lS){
                help[r]=sSeq[rR];
                rR++;
            }
            else{
                if(rR>=lSS){
                    help[r]=seq[rL];
                    rL++;
                }
                else{
                    if(seq[rL]<sSeq[rR]){
                        help[r]=seq[rL];
                        rL++;
                    }
                    else{
                        help[r]=sSeq[rR];
                        rR++;
                    }
                }
            }
            r++;
        }
        for(r=0;r<length;r++){
            seq[r]=help[r];
        }
        delete[] help;
    }
    return 0;
}



myint powerGC(myint base, myint exponent, myint modulo=0){
    if(exponent==1){
        if(modulo!=0){
            return base%modulo;
        }
        else{
            return base;
        }
    }
    if(exponent==0){
        return 1;
    }
    myint forReturn=1;
    myint forReturnH,exp1;
    if(exponent%2==1){
        forReturn=base;
    }
    exp1=exponent/2;
    forReturnH= powerGC(base, exp1, modulo);
    forReturn*=forReturnH;
    forReturn*=forReturnH;
    if(modulo!=0){
        forReturn%=modulo;
    }
    return forReturn;
    
}





mysint GraphicCard::generatePermutationsGC(myint *seqToFill, myint N, myint r){
    if(uniformLimitGC<r){
        uniformLimitGC=r;
    }
    std::uniform_int_distribution<myint> uInt(0,uniformLimitGC);
    myint i,j;
    TripleGC *permuts;
    permuts=new TripleGC[r];
    
    myint seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    seed+=currentSeedGC;
    std::mt19937 mt_randIF(seed);
    
    for(i=0;i<r;i++){
        permuts[i].a=0;
        permuts[i].b=0;
        permuts[i].c=0;
    }
    
    for(i=0;i<N;i++){
        for(j=0;j<r;j++){
            permuts[j].a=uInt(mt_randIF);
            permuts[j].sorting=j;
        }
        mergeSortTGC(permuts,r);
        for(j=0;j<r;j++){
            seqToFill[i*r+j]=permuts[j].sorting;
        }
    }
    
    currentSeedGC+=uInt(mt_randIF);
    delete[] permuts;
    return 1;
    
    
}


myint GraphicCard::createPascalTriangleGC(myint sizeM ){
    if(pascalTriangleGC!=nullptr){
        delete[] pascalTriangleGC;
    }
    
    myint size=sizeM;
    if(size<5){size=5;}
    *sizePascalTrGC=( (size +1)*(size +2) )/2;
    
    pascalTriangleGC=new myint[*sizePascalTrGC];
    myint index1,index2,index3;
    pascalTriangleGC[*sizePascalTrGC-1]=1;
    pascalTriangleGC[0]=1;
    pascalTriangleGC[1]=1;
    
    for(myint n=2;n <=size;n++){
        index1=( n*(n+1) )/2;
        pascalTriangleGC[index1]=1;
        pascalTriangleGC[index1-1]=1;
        index2=index1-n;
        for(myint k =1;k <n; k++){
            //calculating  n choose k
            //Its index in the sequence is k+ n(n+1)/2
            // We will use that \binom nk=\binom{n-1}{k}+\binom{n-1}{k-1}
            
            pascalTriangleGC[index1+k]=pascalTriangleGC[index2+k]+pascalTriangleGC[index2+k-1];
        }
    }
    
    return 1;
    
}

myint GraphicCard::createBalancedNumbersGC(myint *seqBalNum,  myint *cNumber, myint* remZero, myint *remOne, myint *cCounter){
    if((*remZero==0) && (*remOne==0)){
        seqBalNum[*cCounter]=*cNumber;
        *cCounter+=1;
        
        return 1;
    }
    
    if(*remZero>0){
        *cNumber= (*cNumber) * 2;
        *remZero-=1;
        createBalancedNumbersGC(seqBalNum,cNumber,remZero,remOne,cCounter);
        *remZero+=1;
        *cNumber=(*cNumber)/2;
        
    }
    
    if(*remOne>0){
        *cNumber= (*cNumber) * 2 +1 ;
        *remOne-=1;
        createBalancedNumbersGC(seqBalNum,cNumber,remZero,remOne,cCounter);
        *remOne+=1;
        *cNumber=(*cNumber)/2;
        
        
    }
    
    return 1;
}




mysint GraphicCard::generateRandomNumbers(myint N,myint r,mysint inputIndPGen){
    myint indPGen=inputIndPGen;
    
    myint padding;
    //    padding=0;
    padding=*preferred_workgroup_sizeGC;

   
    std::string *lArgRNK;
    myint numArg=12;
    lArgRNK=new std::string[numArg];
    lArgRNK[0]="rNumGCCL";
    lArgRNK[1]="xAxGCCL";
    lArgRNK[2]="yAxGCCL";
    lArgRNK[3]="permutationsXGCCL";
    lArgRNK[4]="permutationsYGCCL";
    lArgRNK[5]="powersOfTwoGCCL";
    lArgRNK[6]="balancedNumbersGCCL";
    lArgRNK[7]="pascalTGCCL";
    lArgRNK[8]="axisSizeGCCL";
    lArgRNK[9]="binaryLengthGCCL";
    lArgRNK[10]="numBalancedNumbersGCCL";
    lArgRNK[11]="shufflingPrimeGCCL";



    if((indPGen==1)||(N!=*axisSizeGC)||(r!=*lengthInBinaryGC)){
        indPGen=1;
        indFirstRNInitGC=0;
        if(r!=*lengthInBinaryGC){
            delete[]balancedNumbersGC;
            balancedNumbersGC=nullptr;
            *numBalancedNumbersGC=0;
            delete[]pascalTriangleGC;
            pascalTriangleGC=nullptr;
            *sizePascalTrGC=0;
        }
        
        for(myint i=0;i<numArg;i++){
            freeDeviceMemory(lArgRNK[i]);
        }
        
        
        *lengthInBinaryGC=r;
        *axisSizeGC=N;
        if(xAxisGC!=nullptr){
            delete[] xAxisGC;
        }
        xAxisGC=new myint[*axisSizeGC+padding];
        if(yAxisGC!=nullptr){
            delete[] yAxisGC;
        }
        yAxisGC=new myint[*axisSizeGC+padding];
        
        uniformLimitGC=*numBalancedNumbersGC-1;
        if(randNumbersGC==nullptr){
            delete[] randNumbersGC;
        }
        randNumbersGC=new myint[ (*axisSizeGC+padding)*(*axisSizeGC+padding)];
        for(myint i=0;i<(*axisSizeGC+padding)*(*axisSizeGC+padding);i++){
            randNumbersGC[i]=17;
        }
        if(xPermutationsGC!=nullptr){
            delete[] xPermutationsGC;
        }
        xPermutationsGC=new myint[(*axisSizeGC+padding)* (*lengthInBinaryGC)];
        if(yPermutationsGC!=nullptr){
            delete[] yPermutationsGC;
        }
        yPermutationsGC=new myint[(*axisSizeGC+padding)* (*lengthInBinaryGC)];
        
        if(powersOfTwoGC!=nullptr){
            delete[] powersOfTwoGC;
        }
        powersOfTwoGC=new myint[*lengthInBinaryGC];
        powersOfTwoGC[0]=1;
        for(myint k=1;k<*lengthInBinaryGC;k++){
            powersOfTwoGC[k]=2 * powersOfTwoGC[k-1];
        }
        
        
        generatePermutationsGC(xPermutationsGC,*axisSizeGC+padding,*lengthInBinaryGC);
        generatePermutationsGC(yPermutationsGC,*axisSizeGC+padding,*lengthInBinaryGC);
        
    }
    if(balancedNumbersGC==nullptr){
        myint **binMatrix;
        myint bigInt=*lengthInBinaryGC+2;
        binMatrix=new myint*[bigInt];
        for(myint i=0;i<bigInt;i++){
            binMatrix[i]=new myint[bigInt];
            for(myint j=0;j<bigInt;j++){
                binMatrix[i][j]=-1;
            }
        }
        
        createPascalTriangleGC(*lengthInBinaryGC+3);
        
        *numBalancedNumbersGC=pascalTriangleGC[ ((*lengthInBinaryGC)*(*lengthInBinaryGC+1))/2+ (*lengthInBinaryGC)/2];
        
        uniformLimitGC=*numBalancedNumbersGC-1;
        for(myint i=0;i<bigInt;i++){
            delete[] binMatrix[i];
        }
        delete[] binMatrix;
        balancedNumbersGC=new myint[*numBalancedNumbersGC];
        myint cNumber=0;
        myint remZero=(*lengthInBinaryGC)/2;
        myint remOne=remZero;
        myint cCounter=0;
        createBalancedNumbersGC(balancedNumbersGC,  &cNumber, &remZero, &remOne, &cCounter);
        
    }
    
    
    
    
    
    std::uniform_int_distribution<myint> uInt(0,uniformLimitGC);
    myint seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    seed+=currentSeedGC;
    std::mt19937 mt_randIF(seed);
    for(myint i=0;i<*axisSizeGC+padding;i++){
        
        xAxisGC[i]= uInt(mt_randIF);
        yAxisGC[i]= uInt(mt_randIF);
        
    }

    currentSeedGC+=uInt(mt_randIF);
    
    writeDeviceMemory("xAxGCCL",xAxisGC,*axisSizeGC+padding);
    writeDeviceMemory("yAxGCCL",yAxisGC,*axisSizeGC+padding);
    
    if(indPGen==1){
        writeDeviceMemory("permutationsXGCCL",xPermutationsGC,(*axisSizeGC+padding)*(*lengthInBinaryGC));
        writeDeviceMemory("permutationsYGCCL",yPermutationsGC,(*axisSizeGC+padding)*(*lengthInBinaryGC));
        writeDeviceMemory("powersOfTwoGCCL",powersOfTwoGC,*lengthInBinaryGC);
        writeDeviceMemory("balancedNumbersGCCL",balancedNumbersGC,*numBalancedNumbersGC);
        
        writeDeviceMemory("binaryLengthGCCL",lengthInBinaryGC,1);
        writeDeviceMemory("shufflingPrimeGCCL",shufflingPrimeGC,1);
        writeDeviceMemory("numBalancedNumbersGCCL",numBalancedNumbersGC,1);
        writeDeviceMemory("pascalTGCCL",pascalTriangleGC,*sizePascalTrGC);
        
        
        writeDeviceMemory("axisSizeGCCL",axisSizeGC,1);
        
    }
    
    
    
    
    
    if(indFirstRNInitGC==0){
        indPGen=1;
        indFirstRNInitGC=1;
        writeDeviceMemory("rNumGCCL",randNumbersGC,(*axisSizeGC+padding)*(*axisSizeGC+padding));
        
        findAddKernel("genMainRandomMatrixGC");
        

        setKernelArguments("genMainRandomMatrixGC",lArgRNK,numArg);
        
        
    }
    delete[] lArgRNK;
    
    
    executeKernel("genMainRandomMatrixGC", (*axisSizeGC)*(*axisSizeGC));
    
    return 1;
}

mysint GraphicCard::generateNormalBeasleySpringerMoro(myint N,myint r,  myint prec1ReqInput, myint *sampleLength,
                                                      myint *overLoadSeq, myint overLoadLen,//these two are for debugging
                                     cl_double par1, cl_double par2, mysint inputIndPGen){
    generateRandomNumbers(N,r,inputIndPGen);
    mysint inspResFun=0;
    myint precReq=  prec1ReqInput;
    
    myint type=7;
    
    
    *highestPrecisionGC=powerGC(*numBalancedNumbersGC, prec1ReqInput)-1;
 
    *precisionRequestGC=precReq;
    
    *sizeForRejectionSamplingGC=precReq;
 
    *sampleLength= (N * N) / (*sizeForRejectionSamplingGC);
    if((*sampleLength!=*sampleLengthGC)||(lastRandAlgUsedGC!=type)){
        lastRandAlgUsedGC=type;
        
        delete[] randomSampleGC;
        inspectorExecutedGC=0;
        sampleGenKernExecutedGC=0;
        *sampleLengthGC=*sampleLength;
        randomSampleGC=new cl_double[*sampleLengthGC];
        
        writeDeviceMemory("randSampGCCL",randomSampleGC,*sampleLengthGC);
        writeDeviceMemory("sampleLengthGCCL",sampleLengthGC,1);
        writeDeviceMemory("sizeRejectionSamplingGCCL",sizeForRejectionSamplingGC,1);
        writeDeviceMemory("precReqGCCL",precisionRequestGC,1);
        writeDeviceMemory("par1GCCL",&par1,1);
        writeDeviceMemory("par2GCCL",&par2,1);
        writeDeviceMemory("numBalancedNumbersGCCL",numBalancedNumbersGC,1);
        writeDeviceMemory("biggestNumGCCL",highestPrecisionGC,1);
        
        cl_double *a,*b,*c;
        myint numAB=4, numC=9;
        a=new cl_double[numAB];
        b=new cl_double[numAB];
        c=new cl_double[numC];
        a[0]=2.50662823884;
        a[1]=-18.61500062529;
        a[2]=41.39119773534;
        a[3]=-25.44106049637;
        
        b[0]=-8.47351093090;
        b[1]=23.08336743743;
        b[2]=-21.06224101826;
        b[3]=3.13082909833;
        
        c[0]=0.3374754822726147;
        c[1]=0.9761690190917186;
        c[2]=0.1607979714918209;
        c[3]=0.0276438810333863;
        c[4]=0.0038405729373609;
        c[5]=0.0003951896511919;
        c[6]=0.0000321767881768;
        c[7]=0.0000002888167364;
        c[8]=0.0000003960315187;
        writeDeviceMemory("seqAGCCL",a,numAB);
        writeDeviceMemory("seqBGCCL",b,numAB);
        writeDeviceMemory("seqCGCCL",c,numC);
        writeDeviceMemory("numABGCCL",&numAB,1);
        writeDeviceMemory("numCGCCL",&numC,1);
        delete[]a;
        delete[]b;
        delete[]c;
    }
    
   
    std::string *lArgRNK;
    myint numArg=14;
    lArgRNK=new std::string[numArg];
    lArgRNK[0]="randSampGCCL";
    lArgRNK[1]="rNumGCCL";
    lArgRNK[2]="sampleLengthGCCL";
    lArgRNK[3]="sizeRejectionSamplingGCCL";
    lArgRNK[4]="precReqGCCL";
    lArgRNK[5]="par1GCCL";
    lArgRNK[6]="par2GCCL";
    lArgRNK[7]="numBalancedNumbersGCCL";
    lArgRNK[8]="biggestNumGCCL";
    lArgRNK[9]="seqAGCCL";
    lArgRNK[10]="seqBGCCL";
    lArgRNK[11]="seqCGCCL";
    lArgRNK[12]="numABGCCL";
    lArgRNK[13]="numCGCCL";
    
    writeDeviceMemory("rNumGCCL",overLoadSeq,overLoadLen);
    
    if(sampleGenKernExecutedGC==0){
        sampleGenKernExecutedGC=1;
        findAddKernel(normalBSMNameGC);
        setKernelArguments(normalBSMNameGC,lArgRNK,numArg);
        
    }
    
    executeKernel(normalBSMNameGC,*sampleLengthGC);
 
    
    delete[] lArgRNK;
 
    
    return 1;
    
}





mysint GraphicCard::generateExponential(myint N,myint r,  myint prec1ReqInput, myint *sampleLength,
                                                    myint *overLoadSeq, myint overLoadLen,//these two are for debugging
                                                      cl_double lambda, mysint inputIndPGen){
    generateRandomNumbers(N,r,inputIndPGen);
    mysint inspResFun=0;
    myint precReq=  prec1ReqInput;
    
    myint type=8;
    
    
    *highestPrecisionGC=powerGC(*numBalancedNumbersGC, prec1ReqInput)-1;
    
    *precisionRequestGC=precReq;
    
    *sizeForRejectionSamplingGC=precReq;
    
    *sampleLength= (N * N) / (*sizeForRejectionSamplingGC);
    if((*sampleLength!=*sampleLengthGC)||(lastRandAlgUsedGC!=type)){
        lastRandAlgUsedGC=type;
        
        delete[] randomSampleGC;
        sampleGenKernExecutedGC=0;
        *sampleLengthGC=*sampleLength;
        randomSampleGC=new cl_double[*sampleLengthGC];
        
        writeDeviceMemory("randSampGCCL",randomSampleGC,*sampleLengthGC);
        writeDeviceMemory("sampleLengthGCCL",sampleLengthGC,1);
        writeDeviceMemory("sizeRejectionSamplingGCCL",sizeForRejectionSamplingGC,1);
        writeDeviceMemory("precReqGCCL",precisionRequestGC,1);
        writeDeviceMemory("par2GCCL",&lambda,1);
 
        writeDeviceMemory("numBalancedNumbersGCCL",numBalancedNumbersGC,1);
        writeDeviceMemory("biggestNumGCCL",highestPrecisionGC,1);
        
    }
    
    
    std::string *lArgRNK;
    myint numArg=8;
    lArgRNK=new std::string[numArg];
    lArgRNK[0]="randSampGCCL";
    lArgRNK[1]="rNumGCCL";
    lArgRNK[2]="sampleLengthGCCL";
    lArgRNK[3]="sizeRejectionSamplingGCCL";
    lArgRNK[4]="precReqGCCL";
    lArgRNK[5]="par2GCCL";
    lArgRNK[6]="numBalancedNumbersGCCL";
    lArgRNK[7]="biggestNumGCCL";
    
    
    writeDeviceMemory("rNumGCCL",overLoadSeq,overLoadLen);
    
    if(sampleGenKernExecutedGC==0){
        sampleGenKernExecutedGC=1;
        findAddKernel(exponentialDistNameGC);
        setKernelArguments(exponentialDistNameGC,lArgRNK,numArg);
        
    }
    
    executeKernel(exponentialDistNameGC,*sampleLengthGC);
    
    
    delete[] lArgRNK;
    
    
    return 1;
    
}







