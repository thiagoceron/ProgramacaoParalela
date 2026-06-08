#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

int main(int argc, char** argv) {
    int id_processo, total_processos;
    long long pontos_totais = 10000000;
    long long meus_pontos, meus_acertos = 0, acertos_totais = 0;
    double x, y, pi_estimado;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id_processo);
    MPI_Comm_size(MPI_COMM_WORLD, &total_processos);

    meus_pontos = pontos_totais / total_processos;
    srand(time(NULL) + id_processo * 1000);

    for (long long i = 0; i < meus_pontos; i++) {
        x = (double)rand() / RAND_MAX;
        y = (double)rand() / RAND_MAX;
        if (x * x + y * y <= 1.0) {
            meus_acertos++;
        }
    }

    // Comunicacao Coletiva
    MPI_Reduce(&meus_acertos, &acertos_totais, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (id_processo == 0) {
        pi_estimado = 4.0 * ((double)acertos_totais / (double)pontos_totais);
        printf("Pi por Monte Carlo (Coletiva) = %.15f\n", pi_estimado);
    }

    MPI_Finalize();
    return 0;
}