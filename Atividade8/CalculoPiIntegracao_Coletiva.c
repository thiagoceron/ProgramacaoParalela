#include <stdio.h>
#include <mpi.h>

int main(int argc, char** argv) {
    int id_processo, total_processos;
    long long passos = 1000000;
    double base = 1.0 / (double)passos;
    double soma = 0.0, x, meu_pi, pi_total;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id_processo);
    MPI_Comm_size(MPI_COMM_WORLD, &total_processos);

    for (long long i = id_processo; i < passos; i += total_processos) {
        x = base * ((double)i + 0.5);
        soma += 4.0 / (1.0 + x * x);
    }
    meu_pi = base * soma;

    // Comunicacao Coletiva
    MPI_Reduce(&meu_pi, &pi_total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (id_processo == 0) {
        printf("Pi por Integracao (Coletiva) = %.15f\n", pi_total);
    }

    MPI_Finalize();
    return 0;
}