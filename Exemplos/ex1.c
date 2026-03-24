#include "omp.h"
#include <stdio.h>
#define tamanho 80

int main()
{
    int a[tamanho], b[tamanho], c[tamanho];

    for (int i = 0; i < tamanho; i++)
    {
        a[i] = i;
        b[i] = 2 * i;
    }
#pragma omp parallel

    {
        int ID = omp_get_thread_num();
        int numero_threads = omp_get_num_threads();
        
        int carga_trabalho = (tamanho + numero_threads - 1) / numero_threads;

        int inicio = ID * carga_trabalho;
        int fim = inicio + carga_trabalho;

        for (int i = inicio; (i < fim) && (i < tamanho); i++)
        {
            c[i] = a[i] + b[i];
            printf("Thread %d calculou c[%d] = %d\n", ID, i, c[i]);
        }
    }  
}