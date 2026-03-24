#include <stdio.h>
#include <omp.h>

double calcular_funcao(double x) {
    return 4.0 / (1.0 + x * x);
}

int main() {
    long long num_passos = 1000000000; 
    double passo = 1.0 / num_passos;
    double pi;
    double soma = 0.0;
    
    #pragma omp parallel for reduction(+:soma)
    for (long long i = 0; i < num_passos; i++) {
        double x = (i + 0.5) * passo;
        soma = soma + calcular_funcao(x);
    }

    pi = passo * soma;

    printf("Valor de Pi calculado: %.15f\n", pi);
    printf("Threads usadas: %d\n", omp_get_max_threads());

    return 0;
}