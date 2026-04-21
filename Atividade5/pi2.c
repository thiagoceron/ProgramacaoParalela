#include <stdio.h>
#include <omp.h>

#pragma omp declare simd uniform(width) linear(i)
double integracao(int i, double width) {
    double x = (i + 0.5) * width;
    return 4.0 / (1.0 + x * x);
}

int main(){
    int n = 1000000;
    double width = 1.0/n;
    double sum = 0.0;

    #pragma omp parallel
    {
        #pragma omp single
        {
            #pragma omp taskloop simd reduction(+:sum) simdlen(8) grainsize(2048)
            for (int i=0; i<n; i++){
                sum += integracao(i, width);
            }
        }
    }

    double pi = sum * width;
    printf("Estimated Pi = %.15f\n", pi);

    return 0;
}