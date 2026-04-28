#include <stdio.h>
#include <omp.h>

#define M 2147483647LL
#define A 48271LL
#define THREADS 8
#define PONTOS 10000000

long long modpow(long long base, long long exp)
{
    long long res = 1;
    base = base % M;
    while (exp > 0)
    {
        if (exp % 2 == 1)
            res = (res * base) % M;
        base = (base * base) % M;
        exp /= 2;
    }
    return res;
}

double pi_leapfrog()
{
    int dentro = 0;
    long long semente_base = 12345;
    long long A_leap = modpow(A, THREADS);

#pragma omp parallel for reduction(+ : dentro) num_threads(THREADS)
    for (int t = 0; t < THREADS; t++)
    {
        long long seed = (semente_base * modpow(A, t)) % M;

        for (int i = 0; i < PONTOS / THREADS; i++)
        {
            double x = (double)seed / M;
            seed = (seed * A_leap) % M;
            double y = (double)seed / M;
            seed = (seed * A_leap) % M;
            if (x * x + y * y <= 1.0)
                dentro++;
        }
    }

    return 4.0 * dentro / PONTOS;
}

double pi_modificado()
{
    int dentro = 0;
    long long semente_base = 12345;
    long long tamanho_bloco = (PONTOS / THREADS) * 2;

#pragma omp parallel for reduction(+ : dentro) num_threads(THREADS)
    for (int t = 0; t < THREADS; t++)
    {
        long long seed = (semente_base * modpow(A, t * tamanho_bloco)) % M;

        for (int i = 0; i < PONTOS / THREADS; i++)
        {
            double x = (double)seed / M;
            seed = (seed * A) % M;

            double y = (double)seed / M;
            seed = (seed * A) % M;

            if (x * x + y * y <= 1.0)
                dentro++;
        }
    }

    return 4.0 * dentro / PONTOS;
}

int main()
{
    double pi1 = pi_leapfrog();
    double pi2 = pi_modificado();

    printf("Pi Leapfrog = %.6f\n", pi1);
    printf("Pi Modificado = %.6f\n", pi2);

    return 0;
}