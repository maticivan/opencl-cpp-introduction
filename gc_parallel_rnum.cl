
typedef int mysint;
typedef int myint;


myint properlyShuffleInsideGC(__global myint *permutationX,
                              __global myint *permutationY,
                              __global myint *powersOfTwo,
                              myint r, myint p,
                              myint i, myint j, myint xA, myint yA, myint xHelp, myint yHelp){
//The function returns \mu_i(x|_p y)+\pi_j(y|_p x)
// x|_p y is the number obtained by taking the lowest p bits of x and the highest r-p bits of y.
    
    myint fR;
    
    myint xShuffled=0;
    myint yShuffled=0;
    myint nextDigitX,nextDigitY;
    myint k=0;myint positionForNextDigitX, positionForNextDigitY;
    myint biggestPowerOfTwo=powersOfTwo[r-1];
    biggestPowerOfTwo*=2;
    myint x=xA;myint y=yA; myint x1=xHelp; myint y1=yHelp;
    //calculation of xShuffled=\mu_i(x|_p y) and yShuffled=p_j(y|_p x)
    
    
    myint jr=j*r; myint ir=i*r;
   // p=r;
    for(k=0;k<r;k++){
        if(k<p){
            nextDigitX=x % 2;
            nextDigitY=y % 2;
        }
        else{
            nextDigitX=x1 % 2;
            nextDigitY=y1 % 2;
        }
        positionForNextDigitX=permutationY[jr+k];
        
        positionForNextDigitY=permutationX[ir+k];
        
        xShuffled+=nextDigitX * powersOfTwo[positionForNextDigitX];
        yShuffled+=nextDigitY * powersOfTwo[positionForNextDigitY];
        //xShuffled+=nextDigitX* (powersOfTwo[k]);
        //yShuffled+=nextDigitY* (powersOfTwo[k]);
        //k++;
        x=x/2;
        y=y/2;
        x1=x1/2;
        y1=y1/2;
    }
    
    //calculation of yShuffled=\mu_i(y)
    

    
    
    fR= (xShuffled+yShuffled)%biggestPowerOfTwo;
    
    return fR;
    
}

__kernel void genMainRandomMatrixGC(__global myint *rNumGCCL,
                __global  myint *xAxGCCL,
                __global  myint *yAxGCCL,
                __global  myint *permutationsXGCCL,
                __global  myint *permutationsYGCCL,
                __global myint *powersOfTwoGCCL,
                __global myint *axisSizeGCCL,
                __global myint *binaryLengthGCCL,
                __global myint *shufflingPrimeGCCL)
{   myint gid = get_global_id(0);
    if(gid<(*axisSizeGCCL)*(*axisSizeGCCL)){
        myint xCoord=gid % (*axisSizeGCCL);
        myint yCoord=gid / (*axisSizeGCCL);
        //rNumGCCL[gid]=*binaryLengthGCCL;
        
        //rNumGCCL[gid]=xAxGCCL[xCoord]+yAxGCCL[yCoord];
        
        myint xH, yH;
        xH=yAxGCCL[0];
        yH=xAxGCCL[0];
        if(yCoord<*axisSizeGCCL-1){
            xH=yAxGCCL[yCoord+1];
        }
        if(xCoord<*axisSizeGCCL-1){
            yH=xAxGCCL[xCoord+1];
        }
        
        rNumGCCL[gid]=properlyShuffleInsideGC(permutationsXGCCL, permutationsYGCCL,powersOfTwoGCCL,
                                              *binaryLengthGCCL, *shufflingPrimeGCCL,
                                              xCoord, yCoord, xAxGCCL[xCoord],yAxGCCL[yCoord],xH,yH);
    }
}
