
__kernel void X(__global  int *a,
                __global  int *b,
                __global int *n)
{   int index = get_global_id(0);
    if(index<*n){
        if (index < (*n)-1) {
            b[index]=a[index]+2*a[index+1];
        }
        else {
            b[index]=a[index];
        }
    }
}

__kernel void Y(__global  int *a,
                __global  int *b,
                __global int *n)
{   int index = get_global_id(0);
    if(index<*n){
        if (index > 0) {
            b[index] = a[index] + 3*a[index-1];
        }
        else {
            b[index] = a[index];
        }
    }
}
