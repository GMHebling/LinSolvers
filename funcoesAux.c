#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "funcoesAux.h"

//------------------------------------------------------------------------------
//FUNÇÕES AUXILIARES
//
//------------------------------------------------------------------------------
//Aloca matriz double de tamanho m x n
double **alocaMatriz(int m, int n)
{
    int i, j;
    double **A;

    A = (double **)malloc(m * sizeof(double *));
    for (i = 0; i < m; i++)
        A[i] = (double *)malloc(n * sizeof(double));

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            A[i][j] = 0;
        }
    }
    return (A);
}
//Aloca vetor double de tamanho m
double *alocaVetor(int m)
{
    int i;
    double *A;

    A = (double *)malloc(m * sizeof(double));

    for (i = 0; i < m; i++)
    {
        A[i] = 0;
    }
    return (A);
}
//Imprime na tela matriz double double de tamanho m x n
void imprimeMatriz(double **A, int m, int n)
{
    int i, j;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            printf("%.8f\t", A[i][j]);
        }
        printf("\n");
    }
}
//Imprime na tela vetor double double de tamanho m
void imprimeVetor(double *A, int m)
{
    int i;
    //FILE* mat;
    printf("Vetor: \n");
    for (i = 0; i < m; i++)
    {
        printf("%.8f\n", A[i]);
    }
}
//------------------------------------------------------------------------------

//Leitura da matriz A e do vetor b
void leituraSistemaLinear(double ***A, double **b, int *n)
{
    int i, j, nvar;
    double valor;
    FILE *arquivo = NULL;

    arquivo = fopen("sist_linear.txt", "r");
    if (arquivo == NULL)
    {
        printf("Arquivo não encontrado\n");
    }
    else
    {
        fscanf(arquivo, "%d\n", &nvar);
        *n = nvar;
        *b = alocaVetor(nvar);
        *A = alocaMatriz(nvar, nvar);

        for (i = 0; i < nvar; i++)
        {
            fscanf(arquivo, "%lf\n", &valor);
            (*b)[i] = valor;
        }
        for (i = 0; i < nvar; i++)
        {
            for (j = 0; j < nvar; j++)
            {
                fscanf(arquivo, "%lf\t", &valor);
                (*A)[i][j] = valor;
            }
            fscanf(arquivo, "\n");
        }
        printf("Fim da leitura do arquivo...\n");
    }
}

double dotProduct(double **A, double **B, int n)
{
    double soma = 0.0;
    int i;
    for (i = 0; i < n; i++)
    {
        soma += (*A)[i] * (*B)[i];
    }
    return soma;
}

double *MultVetorMatriz(double ***A, double **B, int n)
{
    double *result = NULL;
    result = alocaVetor(n);
    int i, j;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            result[i] += (*A)[i][j] * (*B)[j];
        }
    }
    return result;
}

double **MultMatriz(double ***A, double ***B, int n)
{
    double **result = NULL;
    result = alocaMatriz(n, n);

    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            result[i][j] += (*A)[i][j] * (*B)[j][i];
        }
    }
    return result;
}

double *sumVector(double **A, double **B, int n){
    int i;
    double *result;

    result = 0;
    for (i=0;i<n;i++){
        result[i] += (*A)[i] + (*B)[i];
    }
    return result;
}
