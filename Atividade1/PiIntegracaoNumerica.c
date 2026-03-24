#include <stdio.h>
#include <omp.h>

double f(double x){
    return 4.0 / (1.0 + x*x);
}

int main(){

    int n = 100000;
    double a = 0.0, b = 1.0;
    double h = (b-a)/n;
    double pi = 0.0;

    int num_threads = omp_get_max_threads();
    double soma_parcial[num_threads];

    omp_set_num_threads(num_threads);

    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        int nthreads = omp_get_num_threads();

        double soma_local = 0.0;

        for(int i = id; i < n; i += nthreads){
            double x = a + i*h;
            soma_local += f(x);
        }

        soma_parcial[id] = soma_local;
    }

    double soma = 0.0;

    for(int i = 0; i < num_threads; i++){
        soma += soma_parcial[i];
    }

    pi = h * soma;

    printf("Valor aproximado de pi: %.10f\n", pi);

    return 0;
}