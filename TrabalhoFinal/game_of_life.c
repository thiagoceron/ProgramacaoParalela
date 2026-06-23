#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Configuracoes da simulacao
#define TOTAL_PROCS    16   // quantidade exata de processos MPI
#define GRADE_LADO      4   // processos por dimensao da grade 4x4
#define BORDA           5   // espessura das ghost cells em cada lado
#define TAM_GLOBAL    40   // dimensao da matriz global NxN
#define MAX_GERACOES  300   // geracoes simuladas por padrao
#define FREQ_IMPRESSAO  1   // exibe o estado a cada N geracoes

// Estados possiveis de uma celula
#define ATIVA   1
#define INATIVA 0

// Acesso a posicao (i, j) em vetor 1D com nc colunas
#define POS(i, j, nc)  ((i) * (nc) + (j))

// Dimensoes calculadas em main() com base na topologia
static int DIMS[2];       // grade de processos: DIMS[0] x DIMS[1]
static int lin_reais;     // linhas reais por processo
static int col_reais;     // colunas reais por processo
static int lin_local;     // linhas totais por processo (reais + 2*BORDA)
static int col_local;     // colunas totais por processo (reais + 2*BORDA)

// Exibe a divisao dos processos na grade 4x4
void exibir_grade_processos() {
    printf("\n- DIVISAO DOS PROCESSOS -\n\n");
    printf("+----+----+----+----+\n");
    printf("| P0 | P1 | P2 | P3 |\n");
    printf("+----+----+----+----+\n");
    printf("| P4 | P5 | P6 | P7 |\n");
    printf("+----+----+----+----+\n");
    printf("| P8 | P9 |P10 |P11 |\n");
    printf("+----+----+----+----+\n");
    printf("|P12 |P13 |P14 |P15 |\n");
    printf("+----+----+----+----+\n\n");
}

// Exibe o estado global da simulacao no processo 0
void exibir_estado_global(int *estado, int gen, long populacao) {
    printf("\n- GERACAO %d | POPULACAO TOTAL: %ld -\n\n", gen, populacao);

    for (int i = 0; i < TAM_GLOBAL; i++) {
        for (int j = 0; j < TAM_GLOBAL; j++)
            printf("%c", estado[i * TAM_GLOBAL + j] ? 'O' : '.');
        printf("\n");
    }
}

// Soma as celulas ativas apenas na regiao real do processo.
// Ghost cells ficam de fora porque sao copias das bordas vizinhas
// e seriam contadas em duplicata se incluidas.
long somar_celulas_vivas(int *grade) {
    long total = 0;

    for (int i = BORDA; i < BORDA + lin_reais; i++)
        for (int j = BORDA; j < BORDA + col_reais; j++)
            total += grade[POS(i, j, col_local)];

    return total;
}

// Coleta os blocos reais de todos os processos, monta a matriz
// global e a exibe no processo 0.

// Cada processo extrai sua regiao real (sem ghost cells) para um
// vetor auxiliar. MPI_Gather reunindo tudo no processo 0, que
// reposiciona cada bloco conforme as coordenadas cartesianas do
// processo correspondente.
void coletar_e_exibir(int *grade, int rank, MPI_Comm comm_cart,
                      int gen, long pop_total) {
    int *fatia_local = malloc(lin_reais * col_reais * sizeof(int));

    // Extrai a area real descartando as bordas ghost
    for (int i = 0; i < lin_reais; i++)
        for (int j = 0; j < col_reais; j++)
            fatia_local[i * col_reais + j] =
                grade[POS(i + BORDA, j + BORDA, col_local)];

    int *dados_reunidos = NULL;

    // Apenas o processo 0 aloca espaco para todos os blocos
    if (rank == 0)
        dados_reunidos = malloc(TAM_GLOBAL * TAM_GLOBAL * sizeof(int));

    // MPI_Gather coleta lin_reais*col_reais inteiros de cada processo
    MPI_Gather(
        fatia_local,   lin_reais * col_reais, MPI_INT,
        dados_reunidos, lin_reais * col_reais, MPI_INT,
        0, comm_cart
    );

    // Processo 0 reconstroi e imprime a matriz global
    if (rank == 0) {
        int *estado_global = malloc(TAM_GLOBAL * TAM_GLOBAL * sizeof(int));

        for (int p = 0; p < TOTAL_PROCS; p++) {
            int pos_proc[2];

            // Recupera a coordenada cartesiana do processo p
            MPI_Cart_coords(comm_cart, p, 2, pos_proc);

            int ini_lin = pos_proc[0] * lin_reais;
            int ini_col = pos_proc[1] * col_reais;

            // Insere o bloco do processo p na posicao correta
            for (int i = 0; i < lin_reais; i++)
                for (int j = 0; j < col_reais; j++)
                    estado_global[(ini_lin + i) * TAM_GLOBAL + (ini_col + j)] =
                        dados_reunidos[p * lin_reais * col_reais + i * col_reais + j];
        }

        exibir_estado_global(estado_global, gen, pop_total);
        free(estado_global);
        free(dados_reunidos);
    }

    free(fatia_local);
}

// Avanca uma geracao aplicando as quatro regras do Conway
// exclusivamente nas celulas reais. As ghost cells garantem
// que as bordas tenham vizinhanca correta sem condicionais extras.
//
// Regras:
//   Celula ativa com 2 ou 3 vizinhos: sobrevive
//   Celula ativa com menos de 2: morre por solidao
//   Celula ativa com mais de 3: morre por superpopulacao
//   Celula inativa com exatamente 3: nasce
void aplicar_regras(int *grade, int *nova_grade) {
    for (int i = BORDA; i < BORDA + lin_reais; i++) {
        for (int j = BORDA; j < BORDA + col_reais; j++) {
            int viz = 0;

            // Vizinhanca de Moore: 8 celulas ao redor
            for (int di = -1; di <= 1; di++)
                for (int dj = -1; dj <= 1; dj++) {
                    if (di == 0 && dj == 0) continue;
                    viz += grade[POS(i + di, j + dj, col_local)];
                }

            if (grade[POS(i, j, col_local)] == ATIVA)
                nova_grade[POS(i, j, col_local)] =
                    (viz == 2 || viz == 3) ? ATIVA : INATIVA;
            else
                nova_grade[POS(i, j, col_local)] =
                    (viz == 3) ? ATIVA : INATIVA;
        }
    }

    // Atualiza as celulas reais na grade principal
    for (int i = BORDA; i < BORDA + lin_reais; i++)
        for (int j = BORDA; j < BORDA + col_reais; j++)
            grade[POS(i, j, col_local)] = nova_grade[POS(i, j, col_local)];
}

// Sincroniza as ghost cells com os 8 vizinhos do processo.
//
// Tres tipos derivados descrevem os padroes de memoria:
//   tipo_faixa_h : BORDA linhas x col_reais elementos (stride col_local)
//   tipo_faixa_v : lin_reais linhas x BORDA elementos (stride col_local)
//   tipo_bloco   : BORDA x BORDA elementos (stride col_local)
//
// MPI_Sendrecv envia e recebe simultaneamente, evitando deadlock.
// Vizinhos inexistentes recebem MPI_PROC_NULL; o MPI descarta esses
// envios e recebimentos automaticamente.
void trocar_ghost_cells(int *grade, MPI_Comm comm_cart,
                        int norte, int sul, int oeste, int leste,
                        int noroeste, int nordeste,
                        int sudoeste, int sudeste) {
    MPI_Datatype tipo_faixa_h, tipo_faixa_v, tipo_bloco;

    // Faixa horizontal: BORDA linhas consecutivas de col_reais elementos
    MPI_Type_vector(BORDA, col_reais, col_local, MPI_INT, &tipo_faixa_h);
    MPI_Type_commit(&tipo_faixa_h);

    // Faixa vertical: lin_reais blocos de BORDA elementos com stride col_local
    MPI_Type_vector(lin_reais, BORDA, col_local, MPI_INT, &tipo_faixa_v);
    MPI_Type_commit(&tipo_faixa_v);

    // Bloco de canto: BORDA blocos de BORDA elementos com stride col_local
    MPI_Type_vector(BORDA, BORDA, col_local, MPI_INT, &tipo_bloco);
    MPI_Type_commit(&tipo_bloco);

    // troca de faixas horizontais (norte e sul)

    // Envia faixa superior real para norte; recebe de sul no ghost inferior
    MPI_Sendrecv(
        &grade[POS(BORDA, BORDA, col_local)], 1, tipo_faixa_h, norte, 0,
        &grade[POS(BORDA + lin_reais, BORDA, col_local)], 1, tipo_faixa_h, sul, 0,
        comm_cart, MPI_STATUS_IGNORE
    );

    // Envia faixa inferior real para sul; recebe de norte no ghost superior
    MPI_Sendrecv(
        &grade[POS(BORDA + lin_reais - BORDA, BORDA, col_local)], 1, tipo_faixa_h, sul, 1,
        &grade[POS(0, BORDA, col_local)], 1, tipo_faixa_h, norte, 1,
        comm_cart, MPI_STATUS_IGNORE
    );

    // troca de faixas verticais (oeste e leste)

    // Envia faixa esquerda real para oeste; recebe de leste no ghost direito
    MPI_Sendrecv(
        &grade[POS(BORDA, BORDA, col_local)], 1, tipo_faixa_v, oeste, 2,
        &grade[POS(BORDA, BORDA + col_reais, col_local)], 1, tipo_faixa_v, leste, 2,
        comm_cart, MPI_STATUS_IGNORE
    );

    // Envia faixa direita real para leste; recebe de oeste no ghost esquerdo
    MPI_Sendrecv(
        &grade[POS(BORDA, BORDA + col_reais - BORDA, col_local)], 1, tipo_faixa_v, leste, 3,
        &grade[POS(BORDA, 0, col_local)], 1, tipo_faixa_v, oeste, 3,
        comm_cart, MPI_STATUS_IGNORE
    );

    // troca dos blocos de canto (diagonais)

    // Canto superior esquerdo envia para noroeste; recebe sudeste no ghost inf-dir
    MPI_Sendrecv(
        &grade[POS(BORDA, BORDA, col_local)], 1, tipo_bloco, noroeste, 4,
        &grade[POS(BORDA + lin_reais, BORDA + col_reais, col_local)], 1, tipo_bloco, sudeste, 4,
        comm_cart, MPI_STATUS_IGNORE
    );

    // Canto superior direito envia para nordeste; recebe sudoeste no ghost inf-esq
    MPI_Sendrecv(
        &grade[POS(BORDA, BORDA + col_reais - BORDA, col_local)], 1, tipo_bloco, nordeste, 5,
        &grade[POS(BORDA + lin_reais, 0, col_local)], 1, tipo_bloco, sudoeste, 5,
        comm_cart, MPI_STATUS_IGNORE
    );

    // Canto inferior esquerdo envia para sudoeste; recebe nordeste no ghost sup-dir
    MPI_Sendrecv(
        &grade[POS(BORDA + lin_reais - BORDA, BORDA, col_local)], 1, tipo_bloco, sudoeste, 6,
        &grade[POS(0, BORDA + col_reais, col_local)], 1, tipo_bloco, nordeste, 6,
        comm_cart, MPI_STATUS_IGNORE
    );

    // Canto inferior direito envia para sudeste; recebe noroeste no ghost sup-esq
    MPI_Sendrecv(
        &grade[POS(BORDA + lin_reais - BORDA, BORDA + col_reais - BORDA, col_local)], 1, tipo_bloco, sudeste, 7,
        &grade[POS(0, 0, col_local)], 1, tipo_bloco, noroeste, 7,
        comm_cart, MPI_STATUS_IGNORE
    );

    // Libera os tipos derivados apos o uso
    MPI_Type_free(&tipo_faixa_h);
    MPI_Type_free(&tipo_faixa_v);
    MPI_Type_free(&tipo_bloco);
}

// Reparte a matriz global entre os processos via MPI_Scatter.
// O processo 0 empacota os blocos em ordem de rank cartesiano.
// Cada processo copia o bloco recebido para sua area real local.
void repartir_dados(int *grade_global, int *grade_local,
                    int rank, MPI_Comm comm_cart) {
    int *buffer_total = NULL;

    if (rank == 0) {
        buffer_total = malloc((size_t)TOTAL_PROCS * lin_reais * col_reais * sizeof(int));
        for (int p = 0; p < TOTAL_PROCS; p++) {
            int pos_proc[2];
            MPI_Cart_coords(comm_cart, p, 2, pos_proc);
            for (int i = 0; i < lin_reais; i++)
                for (int j = 0; j < col_reais; j++)
                    buffer_total[p * lin_reais * col_reais + i * col_reais + j] =
                        grade_global[POS(pos_proc[0] * lin_reais + i,
                                         pos_proc[1] * col_reais + j, TAM_GLOBAL)];
        }
    }

    int *fatia = malloc(lin_reais * col_reais * sizeof(int));
    MPI_Scatter(buffer_total, lin_reais * col_reais, MPI_INT,
                fatia,         lin_reais * col_reais, MPI_INT,
                0, comm_cart);
    if (rank == 0) free(buffer_total);

    // Posiciona os dados recebidos apos as ghost cells superiores e esquerdas
    for (int i = 0; i < lin_reais; i++)
        for (int j = 0; j < col_reais; j++)
            grade_local[POS(BORDA + i, BORDA + j, col_local)] =
                fatia[i * col_reais + j];
    free(fatia);
}

// Configura o R-pentomino centralizado no grid como estado inicial.
// O R-pentomino e um padrao de apenas 5 celulas que evolui por centenas
// de geracoes antes de estabilizar, gerando grande variedade visual.
// E completamente diferente do quadrado 3x3 usado na referencia.

// Formato do R-pentomino (. = inativa, O = ativa):
//   .OO
//   OO.
//   .O.
void configurar_padrao(int *grade_global) {
    int c = TAM_GLOBAL / 2;  // centro do grid

    // Coordenadas do R-pentomino relativas ao centro
    grade_global[POS(c-1, c,   TAM_GLOBAL)] = ATIVA;
    grade_global[POS(c-1, c+1, TAM_GLOBAL)] = ATIVA;
    grade_global[POS(c,   c-1, TAM_GLOBAL)] = ATIVA;
    grade_global[POS(c,   c,   TAM_GLOBAL)] = ATIVA;
    grade_global[POS(c+1, c,   TAM_GLOBAL)] = ATIVA;
}


// Carrega o estado inicial a partir de um arquivo texto.
// Linhas com '!' ou '#' sao tratadas como comentarios e ignoradas.
// Celulas ativas: 'O', '1', '#' ou '*'. Qualquer outro caractere = inativa.
void carregar_arquivo(const char *caminho, int *grade_global) {
    FILE *arq = fopen(caminho, "r");
    if (!arq) {
        fprintf(stderr, "Nao foi possivel abrir: %s\n", caminho);
        return;
    }
    char buf[8192];
    int i = 0;
    while (i < TAM_GLOBAL && fgets(buf, sizeof buf, arq)) {
        if (buf[0] == '!' || buf[0] == '#') continue;
        for (int j = 0; buf[j] && buf[j] != '\n' && j < TAM_GLOBAL; j++) {
            char c = buf[j];
            grade_global[POS(i, j, TAM_GLOBAL)] =
                (c == 'O' || c == '1' || c == '#' || c == '*') ? ATIVA : INATIVA;
        }
        i++;
    }
    fclose(arq);
}

int main(int argc, char **argv) {
    // Inicializa o ambiente MPI
    MPI_Init(&argc, &argv);

    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // O programa requer exatamente TOTAL_PROCS processos para
    // preencher a grade cartesiana 4x4 corretamente
    if (num_procs != TOTAL_PROCS) {
        if (rank == 0)
            printf("Execute com: mpirun -np %d --oversubscribe ./game_of_life\n",
                   TOTAL_PROCS);
        MPI_Finalize();
        return 1;
    }

    // Leitura dos argumentos: geracoes e/ou arquivo de entrada
    int num_geracoes = MAX_GERACOES;
    const char *nome_arquivo = NULL;
    for (int i = 1; i < argc; i++) {
        int v = atoi(argv[i]);
        if (v > 0) num_geracoes = v; else nome_arquivo = argv[i];
    }

    //Topologia cartesiana

    // MPI_Dims_create distribui os processos nas dimensoes da grade.
    // Para 16 processos resulta em DIMS = {4, 4}.

    // MPI_Cart_create cria o comunicador com layout 2D.
    // periodico={0,0}: grade nao-periodica, bordas nao se fecham.

    // MPI_Cart_shift localiza os vizinhos diretos em cada dimensao
    // sem necessidade de calcular os ranks manualmente.
    DIMS[0] = DIMS[1] = 0;
    MPI_Dims_create(TOTAL_PROCS, 2, DIMS);

    int periodico[2] = {0, 0};
    MPI_Comm comm_cart;
    MPI_Cart_create(MPI_COMM_WORLD, 2, DIMS, periodico, 0, &comm_cart);

    int posicao[2];
    MPI_Cart_coords(comm_cart, rank, 2, posicao);

    int pos_lin = posicao[0];
    int pos_col = posicao[1];

    // Vizinhos diretos: dimensao 0 = norte/sul, dimensao 1 = oeste/leste
    int norte, sul, oeste, leste;
    MPI_Cart_shift(comm_cart, 0, 1, &norte, &sul);
    MPI_Cart_shift(comm_cart, 1, 1, &oeste, &leste);

    // Vizinhos diagonais via MPI_Cart_rank.
    // MPI_PROC_NULL marca ausencia de vizinho nas bordas da grade;
    // o MPI descarta automaticamente comunicacoes com esse valor.
    int noroeste = MPI_PROC_NULL, nordeste  = MPI_PROC_NULL;
    int sudoeste = MPI_PROC_NULL, sudeste   = MPI_PROC_NULL;

    int coord_viz[2];

    coord_viz[0] = pos_lin - 1; coord_viz[1] = pos_col - 1;
    if (coord_viz[0] >= 0 && coord_viz[1] >= 0)
        MPI_Cart_rank(comm_cart, coord_viz, &noroeste);

    coord_viz[0] = pos_lin - 1; coord_viz[1] = pos_col + 1;
    if (coord_viz[0] >= 0 && coord_viz[1] < GRADE_LADO)
        MPI_Cart_rank(comm_cart, coord_viz, &nordeste);

    coord_viz[0] = pos_lin + 1; coord_viz[1] = pos_col - 1;
    if (coord_viz[0] < GRADE_LADO && coord_viz[1] >= 0)
        MPI_Cart_rank(comm_cart, coord_viz, &sudoeste);

    coord_viz[0] = pos_lin + 1; coord_viz[1] = pos_col + 1;
    if (coord_viz[0] < GRADE_LADO && coord_viz[1] < GRADE_LADO)
        MPI_Cart_rank(comm_cart, coord_viz, &sudeste);

    // Dimensoes do bloco de cada processo:
    // TAM_GLOBAL / GRADE_LADO = 160 / 4 = 40 celulas reais por lado.
    // A area local inclui as ghost cells ao redor da regiao real.
    lin_reais = TAM_GLOBAL / DIMS[0];
    col_reais = TAM_GLOBAL / DIMS[1];
    lin_local = lin_reais + 2 * BORDA;
    col_local = col_reais + 2 * BORDA;

    // Alocacao das grades local e auxiliar com calloc (inicia tudo em zero)
    int *grade     = calloc(lin_local * col_local, sizeof(int));
    int *nova_grade = calloc(lin_local * col_local, sizeof(int));

    if (rank == 0) {
        printf("Matriz global: %dx%d\n", TAM_GLOBAL, TAM_GLOBAL);
        printf("Processos: %d\n", TOTAL_PROCS);
        printf("Divisao: %dx%d processos\n", DIMS[0], DIMS[1]);
        printf("Cada processo: %dx%d celulas reais\n", lin_reais, col_reais);
        printf("Ghost cells: espessura %d ao redor de cada bloco\n", BORDA);
        printf("Matriz interna de cada processo: %dx%d\n", lin_local, col_local);
        printf("Geracoes: %d\n", num_geracoes);
        exibir_grade_processos();
    }
    MPI_Barrier(comm_cart);

    // Processo 0 monta a grade global e a distribui para todos
    int *grade_global = NULL;
    if (rank == 0) {
        grade_global = calloc(TAM_GLOBAL * TAM_GLOBAL, sizeof(int));
        if (nome_arquivo) carregar_arquivo(nome_arquivo, grade_global);
        else              configurar_padrao(grade_global);
    }
    repartir_dados(grade_global, grade, rank, comm_cart);
    if (rank == 0) free(grade_global);

    double tempo_ini = MPI_Wtime();
    long pop_inicio = 0, pop_fim = 0;

    // Loop principal: geracao 0 e o estado inicial; o programa
    // alterna entre trocar ghost cells e calcular a proxima geracao
    // ate atingir num_geracoes.
    for (int gen = 0; gen <= num_geracoes; gen++) {

        // Comunicacao coletiva: MPI_Reduce
        // Cada processo calcula sua soma local de celulas ativas.
        // MPI_Reduce agrega todos os valores no processo 0 usando MPI_SUM.
        // Ghost cells sao excluidas para evitar dupla contagem.
        long vivos_locais = somar_celulas_vivas(grade);
        long vivos_total  = 0;
        MPI_Reduce(&vivos_locais, &vivos_total, 1, MPI_LONG, MPI_SUM, 0, comm_cart);

        if (rank == 0) {
            if (gen == 0)           pop_inicio = vivos_total;
            if (gen == num_geracoes) pop_fim    = vivos_total;
        }

        // Coleta e exibe o estado global a cada FREQ_IMPRESSAO geracoes
        if (gen % FREQ_IMPRESSAO == 0)
            coletar_e_exibir(grade, rank, comm_cart, gen, vivos_total);

        if (gen == num_geracoes) break;

        // Comunicacao ponto a ponto: atualiza ghost cells
        trocar_ghost_cells(grade, comm_cart,
                           norte, sul, oeste, leste,
                           noroeste, nordeste, sudoeste, sudeste);

        // Aplica as regras do Conway e avanca uma geracao
        aplicar_regras(grade, nova_grade);
    }

    double tempo_fim = MPI_Wtime();

    if (rank == 0) {
        printf("\n- ESTATISTICAS -\n\n");
        printf("Geracoes simuladas : %d\n",   num_geracoes);
        printf("Matriz global      : %dx%d\n", TAM_GLOBAL, TAM_GLOBAL);
        printf("Topologia          : %dx%d processos\n", DIMS[0], DIMS[1]);
        printf("Ghost cells        : %d de espessura\n", BORDA);
        printf("Tempo total        : %.4f s\n", tempo_fim - tempo_ini);
        printf("Throughput         : %.2f geracoes/s\n",
               num_geracoes / (tempo_fim - tempo_ini));
        printf("Populacao inicial  : %ld\n", pop_inicio);
        printf("Populacao final    : %ld\n", pop_fim);
    }

    // Libera os recursos alocados
    free(grade);
    free(nova_grade);
    MPI_Comm_free(&comm_cart);
    MPI_Finalize();

    return 0;
}