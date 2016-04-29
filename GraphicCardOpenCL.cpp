
#include "gcocl_osx_beginning.cpp"


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

   




    if((indPGen==1)||(N!=*axisSizeGC)||(r!=*lengthInBinaryGC)){
        indPGen=1;
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
        
        //std::cout<<"From the graphic card: "<<*numBalancedNumbersGC<< std::endl;
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
        setKernelArguments("genMainRandomMatrixGC",lArgRNK,numArg);
        delete[] lArgRNK;
        
    }
    executeKernel("genMainRandomMatrixGC", (*axisSizeGC)*(*axisSizeGC));
    return 1;
}

