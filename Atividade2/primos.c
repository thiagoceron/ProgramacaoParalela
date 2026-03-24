#include <stdio.h>
#include <stdbool.h>
#include <omp.h>

#define TOTAL 5000

int main() {
    bool primos[TOTAL + 1];

//divide a inicializacao do vetor em paralelo
    for (int i = 0; i <= TOTAL; i++) {
        primos[i] = true;
    }

    primos[0] = false;
    primos[1] = false;

//cria a regiao paralela
    #pragma omp parallel
    {
        for (int j = 2; j * j <= TOTAL; j++){
             if (primos[j]) {
                //cada thread fica responsavel por marcar uma parte dos multiplos como nao primos
                #pragma omp for
                for (int i = j * j; i <= TOTAL; i += j) {
                    primos[i] = false;
                }
            }
        }
        
           
    }

    printf("Numeros primos de 1 ate %d:\n", TOTAL);
    for (int i = 2; i <= TOTAL; i++) {
        if (primos[i]) {
            printf("%d ", i);
        }
    }
    printf("Threads usadas: %d\n", omp_get_max_threads());

    return 0;
}