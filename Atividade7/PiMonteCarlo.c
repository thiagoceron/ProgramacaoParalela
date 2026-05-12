#include <stdio.h>
#include <omp.h>

int main() {
    long numero_tentativas = 1000000000;
    long contagem = 0;
    double tempo_inicio, tempo_fim;

    tempo_inicio = omp_get_wtime();

    #pragma omp target map(tofrom: contagem)
    #pragma omp teams distribute parallel for reduction(+:contagem)
    for (long i = 0; i < numero_tentativas; i++) {
        unsigned int semente = (unsigned int)(i + 1) ^ 0x55555555;
        
        semente ^= semente << 13; semente ^= semente >> 17; semente ^= semente << 5;
        double x = (double)semente / 4294967295.0;
        
        semente ^= semente << 13; semente ^= semente >> 17; semente ^= semente << 5;
        double y = (double)semente / 4294967295.0;

        if (x * x + y * y <= 1.0) {
            contagem++;
        }
    }

    double pi = 4.0 * ((double)contagem / (double)numero_tentativas);
    
    tempo_fim = omp_get_wtime();

    printf("Valor aproximado de Pi: %.15f\n", pi);
    printf("Tempo de execução: %f segundos\n", tempo_fim - tempo_inicio);

    return 0;
}