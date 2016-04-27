
typedef int mysint;
typedef int myint;


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

        //positionForNextDigitX=permutationY[ir+k];
        //positionForNextDigitY=permutationX[jr+k];

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
