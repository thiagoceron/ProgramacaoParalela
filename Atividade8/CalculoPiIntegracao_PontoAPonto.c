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

    // Comunicacao Ponto a Ponto
    if (id_processo == 0) {
        pi_total = meu_pi;
        for (int i = 1; i < total_processos; i++) {
            double pi_temp;
            MPI_Recv(&pi_temp, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            pi_total += pi_temp;
        }
        printf("Pi por Integracao (Ponto a Ponto) = %.15f\n", pi_total);
    } else {
        MPI_Send(&meu_pi, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}