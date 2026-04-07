#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define N 9

int found = 0;

void print(int arr[N][N])
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
            printf("%d ", arr[i][j]);
        printf("\n");
    }
}

int isSafe(int grid[N][N], int row, int col, int num)
{
    for (int x = 0; x <= 8; x++)
        if (grid[row][x] == num)
            return 0;
    for (int x = 0; x <= 8; x++)
        if (grid[x][col] == num)
            return 0;

    int startRow = row - row % 3, startCol = col - col % 3;

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            if (grid[i + startRow][j + startCol] == num)
                return 0;

    return 1;
}

int solveSudoku(int grid[N][N], int row, int col, int depth)
{
    if (found)
        return 1;
    if (row == N - 1 && col == N)
    {
#pragma omp critical
        {
            if (!found)
            {
                found = 1;
                print(grid);
            }
        }
        return 1;
    }
    if (col == N)
    {
        row++;
        col = 0;
    }
    if (grid[row][col] > 0)
        return solveSudoku(grid, row, col + 1, depth);
    for (int num = 1; num <= N; num++)
    {
        if (isSafe(grid, row, col, num) == 1)
        {
            int newGrid[N][N];
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                    newGrid[i][j] = grid[i][j];

            newGrid[row][col] = num;

            // utilizei 5, porque paralelizar apenas os niveis iniciais da arvore ja gera ramificacoes
            // aprofundar mais do que isso faria o OpenMP perder desempenho gerenciando o excesso de tasks
            if (depth < 5)
            {
#pragma omp task firstprivate(newGrid, row, col, depth) final(depth >= 5)
                {
                    solveSudoku(newGrid, row, col + 1, depth + 1);
                }
            }
            else
            {
                solveSudoku(newGrid, row, col + 1, depth + 1);
            }
        }
    }
#pragma omp taskwait

    return found;
}

int main()
{
    int grid[N][N] = {{5, 3, 0, 0, 7, 0, 0, 0, 0},
                      {6, 0, 0, 1, 9, 5, 0, 0, 0},
                      {0, 9, 8, 0, 0, 0, 0, 6, 0},
                      {8, 0, 0, 0, 6, 0, 0, 0, 3},
                      {4, 0, 0, 8, 0, 3, 0, 0, 1},
                      {7, 0, 0, 0, 2, 0, 0, 0, 6},
                      {0, 6, 0, 0, 0, 0, 2, 8, 0},
                      {0, 0, 0, 4, 1, 9, 0, 0, 5},
                      {0, 0, 0, 0, 8, 0, 0, 7, 9}};

#pragma omp parallel
    {
#pragma omp single
        {
            solveSudoku(grid, 0, 0, 0);
        }
    }

    return 0;
}