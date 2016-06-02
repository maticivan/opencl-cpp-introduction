
typedef long int mysint;
typedef long int myint;


myint invBalNum(__global myint *seqPascTr, myint needle){
    
    myint numZeroes=0;
    myint numOnes=0;
    myint needleHelp=needle;
    myint currentDigit;
    myint fR=0; myint helpVar;
    myint m=0;

    while(needleHelp>0){
        currentDigit=needleHelp % 2;
        needleHelp = needleHelp/2;
        
        numZeroes+=(1-currentDigit);
        numOnes+=currentDigit;

        // m is a 0-1 indicator that is 0 until the first digit 0 appears. Once
        // the first 0 appears, m becomes 1 and stays 1.
        m+=(1-currentDigit)*(1-m);
        

        // If currentDigit==1 then we have to add \varpsi(numOnes, numZeroes)
        // \varpsi(numOnes, numZeroes)  = \binom{numOnes+numZeroes-1}{numZeroes-1}
        // \binom(n,k)=seqPascTr( n*(n+1)/2+k)
        helpVar=numOnes+numZeroes-1;
        
        fR+= currentDigit* m * seqPascTr[ numZeroes-1 +( helpVar*(helpVar+1)/2)];
            
       
    }
    
    return fR;

    
}


myint properlyShuffleInsideGC(__global myint *permutationX,
                              __global myint *permutationY,
                              __global myint *powersOfTwo,
                              __global myint *seqPascTr,
                              __global myint *seqBalNum,
                              myint numSeqBalNum,
                              myint r, myint p,
                              myint i, myint j, myint xA, myint yA){
    myint fR;
    
    myint xShuffled=0;
    myint yShuffled=0;
    myint nextDigitX,nextDigitY;
    myint k=0;myint positionForNextDigitX, positionForNextDigitY;
    myint biggestPowerOfTwo=powersOfTwo[r-1];
    biggestPowerOfTwo*=2;
    myint x=seqBalNum[xA];myint y=seqBalNum[yA];
    myint jr=j*r; myint ir=i*r;
    for(k=0;k<r;k++){

        nextDigitX=x % 2;
        nextDigitY=y % 2;
        positionForNextDigitX=permutationY[jr+k];
        positionForNextDigitY=permutationX[ir+k];


        xShuffled+=nextDigitX * powersOfTwo[positionForNextDigitX];
        yShuffled+=nextDigitY * powersOfTwo[positionForNextDigitY];
         x=x/2;
        y=y/2;
    }

    xShuffled=invBalNum(seqPascTr, xShuffled);
    yShuffled=invBalNum(seqPascTr,yShuffled);
    
    fR= (xShuffled+yShuffled)% numSeqBalNum;
    
    return fR;
    
}

__kernel void genMainRandomMatrixGC(__global myint *rNumGCCL,
                __global  myint *xAxGCCL,
                __global  myint *yAxGCCL,
                __global  myint *permutationsXGCCL,
                __global  myint *permutationsYGCCL,
                __global myint *powersOfTwoGCCL,
                __global myint *balancedNumbersGCCL,
                __global myint *pascalTGCCL,
                __global myint *axisSizeGCCL,
                __global myint *binaryLengthGCCL,
                __global myint *numBalancedNumbersGCCL,
                __global myint *shufflingPrimeGCCL)
{   myint gid = get_global_id(0);

        myint xCoord=gid % (*axisSizeGCCL);
        myint yCoord=gid / (*axisSizeGCCL);
        rNumGCCL[gid]=*binaryLengthGCCL;
        

        
        myint xH, yH;
        xH=yAxGCCL[0];
        yH=xAxGCCL[0];

        
        rNumGCCL[gid]=properlyShuffleInsideGC(permutationsXGCCL, permutationsYGCCL,powersOfTwoGCCL,
                                              pascalTGCCL, balancedNumbersGCCL,
                                              *numBalancedNumbersGCCL,
                                              *binaryLengthGCCL, *shufflingPrimeGCCL,
                                              xCoord, yCoord, xAxGCCL[xCoord],yAxGCCL[yCoord]);

}



__kernel void inspectorKernelGC(__global myint *inspectorGCCL,
                                    __global  myint *sampleLengthGCCL,
                                    __global  myint *oldLayerGCCL)
{
    myint gid=get_global_id(0);
    myint writingPosition= gid * (*oldLayerGCCL)*2;
    if(writingPosition+*oldLayerGCCL<*sampleLengthGCCL){
        myint theFirst= inspectorGCCL[writingPosition];
        myint theSecond=inspectorGCCL[writingPosition+*oldLayerGCCL];
        theFirst%=2;
        theSecond%=2;
        inspectorGCCL[writingPosition]=1;
        if((theFirst==0)||(theSecond==0)){
            inspectorGCCL[writingPosition]=0;
        }
    }
}


__kernel void rejectionSamplerExponentialGC(__global double *randSampGCCL,
                                            __global myint *rNumGCCL,
                                            __global myint *inspectorGCCL,
                                            __global myint *sampleLengthGCCL,
                                            __global myint *sizeRejectionSamplingGCCL,
                                            __global myint *precReqGCCL,
                                            __global myint *par1GCCL,
                                            __global myint *par2GCCL,
                                            __global myint *numBalancedNumbersGCCL,
                                            __global myint *biggestNumGCCL)
{
    myint gid=get_global_id(0);
    if(gid<*sampleLengthGCCL){
        myint sL=*sizeRejectionSamplingGCCL;
        myint base=*numBalancedNumbersGCCL;
        myint readPos=gid * sL;
        myint lastReadPos=(gid+1)*sL;
        myint oneNumSize= (*precReqGCCL)/2;
        myint success=0;
        double x,y;
        myint j;
        while((success==0)&&(readPos<lastReadPos)){
            x=0.0;y=0.0;
            for(j=0;j<oneNumSize;j++){
                x=x * ((double) base)+((double) rNumGCCL[readPos+1+j]);
                y=y * ((double) base)+((double) rNumGCCL[readPos+1+j+oneNumSize]);
            }
            x= x/ ((double)*biggestNumGCCL);
            y= y/ ((double)*biggestNumGCCL);
            readPos+=sL;
            success=1;
            j=2;
            x=(double)rNumGCCL[readPos+1+j] +0.1;
        }
        randSampGCCL[gid]=x;
        inspectorGCCL[gid]=success;
        
    }
}




__kernel void rejectionSamplerNormalGC(__global double *randSampGCCL,
                                       __global myint *rNumGCCL,
                                       __global myint *inspectorGCCL,
                                       __global myint *sampleLengthGCCL,
                                       __global myint *sizeRejectionSamplingGCCL,
                                       __global myint *precReqGCCL,
                                       __global myint *par1GCCL,
                                       __global myint *par2GCCL,
                                       __global myint *numBalancedNumbersGCCL,
                                       __global myint *biggestNumGCCL)
{
    myint gid=get_global_id(0);
    if(gid<*sampleLengthGCCL){
        myint sL=*sizeRejectionSamplingGCCL;
        myint base=*numBalancedNumbersGCCL;
        myint readPos=gid * sL;
        myint lastReadPos=(gid+1)*sL;
        myint oneNumSize= (*precReqGCCL)/2;
        myint success=0;
        double x,y;
        myint j;
        while((success==0)&&(readPos<lastReadPos)){
            x=0.0;y=0.0;
            for(j=0;j<oneNumSize;j++){
                x=x * ((double) base)+((double) rNumGCCL[readPos+1+j]);
                y=y * ((double) base)+((double) rNumGCCL[readPos+1+j+oneNumSize]);
            }
            x= x/ ((double)*biggestNumGCCL);
            y= y/ ((double)*biggestNumGCCL);
            readPos+=sL;
            success=1;
            
            
        }
        randSampGCCL[gid]=x;
        inspectorGCCL[gid]=success;
        
    }
}













__kernel void normalBSMGC(__global double *randSampGCCL,
                          __global myint *rNumGCCL,
                          __global myint *sampleLengthGCCL,
                          __global myint *sizeRejectionSamplingGCCL,
                          __global myint *precReqGCCL,
                          __global double *par1GCCL,
                          __global double *par2GCCL,
                          __global myint *numBalancedNumbersGCCL,
                          __global myint *biggestNumGCCL,
                          __global double *seqAGCCL,
                          __global double *seqBGCCL,
                          __global double *seqCGCCL,
                          __global myint* numABGCCL,
                          __global myint* numCGCCL)
{
    myint gid=get_global_id(0);
    if(gid<*sampleLengthGCCL){
        myint sL=*sizeRejectionSamplingGCCL;
        myint base=*numBalancedNumbersGCCL;
        myint readPos=gid * sL;
        
        double x;
        myint i ;
        
        x=0.0;
        for(i=0;i<sL;i++){
            x=x * ((double) base)+((double) rNumGCCL[readPos+i]);
            
        }
        x= x/ ((double)*biggestNumGCCL);
        
        //Z = mu+ sigma F^{-1}(X)
        
        double fInvX;
        
        myint numAB=*numABGCCL, numC=*numCGCCL;
        __global double *a,*b,*c;
        a=seqAGCCL;
        b=seqBGCCL;
        c=seqCGCCL;
 
        double y=x - 0.5;
 
        double absy= y,signy=1.0;
        if(y<0){
            absy=-y;
            signy=-1.0;
        }
        
        double r;
        double numerator=0, denominator=0;
        
        if(absy<0.42){
            
            r= absy * absy;
            for(i=0;i<numAB;i++){
                numerator*= r;
                numerator+=a[numAB-i-1];
            }
            numerator*=y;
            for( i =0;i<numAB;i++){
                denominator*=r;
                denominator+=b[numAB-i-1];
            }
            denominator*=r;
            denominator+=1;
            fInvX=0;
            if(denominator!=0){
                fInvX= numerator/denominator;
            }
            
        }
        else{
            r=x;
            if(x>0.5){
                r=1-x;
            }
            fInvX=0;
            if(r!=0){
                r=log(-(log(r)));
                fInvX=0;
                for(i=0;i<numC;i++){
                    fInvX*=r;
                    fInvX+=c[numC-i-1];
                }
                fInvX*=signy;
            }
        }
        

        fInvX*= (*par2GCCL);
        fInvX+=*par1GCCL;
 
        randSampGCCL[gid]=fInvX;
        

        
    }
}




__kernel void exponentialDistGC(__global double *randSampGCCL,
                          __global myint *rNumGCCL,
                          __global myint *sampleLengthGCCL,
                          __global myint *sizeRejectionSamplingGCCL,
                          __global myint *precReqGCCL,
                          __global double *par2GCCL,
                          __global myint *numBalancedNumbersGCCL,
                          __global myint *biggestNumGCCL)
{
    myint gid=get_global_id(0);
    if(gid<*sampleLengthGCCL){
        myint sL=*sizeRejectionSamplingGCCL;
        myint base=*numBalancedNumbersGCCL;
        myint readPos=gid * sL;
        
        double x;
        myint i ;
        
        x=0.0;
        for(i=0;i<sL;i++){
            x=x * ((double) base)+((double) rNumGCCL[readPos+i]);
            
        }
        x= x/ ((double)*biggestNumGCCL);
        
        //Z = -\frac{1}{par2GCCL} log(1-x)
        
        
        double fInvX=0;
        
        if(x!=1){
            fInvX=-log(1-x);
            fInvX/= *par2GCCL;
        }
  
        
        randSampGCCL[gid]=fInvX;
        
        
        
    }
}


