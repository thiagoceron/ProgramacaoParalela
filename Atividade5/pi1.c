#include <stdio.h>
#include <omp.h>

int main(){
    int n = 1000000;
    double width = 1.0/n;
    double sum = 0.0;
    double x;

    #pragma omp taskloop simd private(x) reduction(+:sum) simdlen(8) grainsize(2048)
    for (int i=0; i<n; i++){
        x = (i + 0.5) * width;
        sum += 4.0/(1.0 + x*x);
    }

    double pi = sum * width;
    printf("Estimated Pi = %.15f\n", pi);

    return 0;
}
