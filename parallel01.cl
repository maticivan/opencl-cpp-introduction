
__kernel void eachTermGetsDoubled(__global int *a, __global int *n)
{
    int lengthOfTheSequence=*n;
    int index = get_global_id(0);
    if(index<lengthOfTheSequence){
        a[index]=a[index]*2;
    }
}

