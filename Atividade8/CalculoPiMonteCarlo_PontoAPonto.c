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

    // Comunicacao Ponto a Ponto
    if (id_processo == 0) {
        acertos_totais = meus_acertos;
        for (int i = 1; i < total_processos; i++) {
            long long acertos_temp;
            MPI_Recv(&acertos_temp, 1, MPI_LONG_LONG, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            acertos_totais += acertos_temp;
        }
        pi_estimado = 4.0 * ((double)acertos_totais / (double)pontos_totais);
        printf("Pi por Monte Carlo (Ponto a Ponto) = %.15f\n", pi_estimado);
    } else {
        MPI_Send(&meus_acertos, 1, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}