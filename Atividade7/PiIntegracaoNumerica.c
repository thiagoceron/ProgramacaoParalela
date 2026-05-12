#include <stdio.h>
#include <omp.h>

int main() {
    long numero_passos = 1000000000;
    double passo = 1.0 / (double)numero_passos;
    double pi = 0.0;
    double tempo_inicio, tempo_fim;

    tempo_inicio = omp_get_wtime();

    #pragma omp target map(tofrom: pi) map(to: passo)
    #pragma omp teams distribute parallel for reduction(+:pi)
    for (long i = 0; i < numero_passos; i++) {
        double x = (i + 0.5) * passo;
        pi += 4.0 / (1.0 + x * x);
    }
    
    pi *= passo;
    
    tempo_fim = omp_get_wtime();

    printf("Valor de Pi: %.15f\n", pi);
    printf("Tempo de execução: %f segundos\n\n", tempo_fim - tempo_inicio);

    return 0;
}