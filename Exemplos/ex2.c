#include "omp.h"
#include "stdio.h"

#define tamanho 80
#define maxThreads 32

int main()
{
    int a[tamanho], resultado = 0, resultado_local[maxThreads] = {0}, numero_threads;

    for (int i = 0; i < tamanho; i++)
    {
        a[i] = i;
    }
#pragma omp parallel

    {
        int ID = omp_get_thread_num();
        numero_threads = omp_get_num_threads();

        int resultado_temp = 0;

        int carga_trabalho = (tamanho + numero_threads - 1) / numero_threads;

        int inicio = ID * carga_trabalho;
        int fim = inicio + carga_trabalho;

        for (int i = inicio; (i < fim) && (i < tamanho); i++)
        {
            resultado_temp += a[i];
        }

        resultado_local[ID] = resultado_temp;

    }
    for (int i = 0; i < numero_threads; i++)
    {
        resultado += resultado_local[i];
    }
    printf("Resultado: %d\n", resultado);
}