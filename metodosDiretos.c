#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "metodosDiretos.h"
#include "funcoesAux.h"
//-----------------------------------------------------------------------------------------
// METODOS DIRETOS DE SOLUCAO DE SISTEMAS LINEARES
//-----------------------------------------------------------------------------------------

double *LUFactorization(double ***A, double **b, int n)
{
    double *result = NULL;
    double *Y = NULL;
    int i, j, k;
    double **L = NULL;
    double **U = NULL;

    result = alocaVetor(n);
    Y = alocaVetor(n);
    L = alocaMatriz(n, n);
    U = alocaMatriz(n, n);
    //elementos da diagonal de L iguais a 1
    for (i = 0; i < n; i++)
    {
        L[i][i] = 1;
    }
    //elemento A[0,0] diferente de zero
    if ((*A)[0][0] == 0)
    {
        return 0;
    }
    else
    {
        U[0][0] = (*A)[0][0] / L[0][0];
    }
    for (j = 1; j < n; j++)
    {
        U[0][j] = (*A)[0][j] / L[0][0];
        L[j][0] = (*A)[j][0] / U[0][0];
    }
    for (i = 1; i < n - 1; i++)
    {
        double soma_diag;
        soma_diag = 0;
        for (k = 0; k < i; k++)
        {
            soma_diag += (L[i][k] * U[k][i]);
        }
        U[i][i] = ((*A)[i][i] - soma_diag) / L[i][i];
        if (U[i][i] == 0)
        {
            return 0;
        }
        for (j = 0; j < n; j++)
        {
            double somaU;
            double somaL;
            somaU = 0;
            somaL = 0;
            for (k = 0; k < i; k++)
            {
                somaU += (L[i][k] * U[k][j]);
                somaL += (L[j][k] * U[k][i]);
            }
            U[i][j] = (1 / L[i][i]) * ((*A)[i][j] - somaU);
            L[j][i] = (1 / U[i][i]) * ((*A)[j][i] - somaL);
        }
    }
    double somaN;
    somaN = 0;
    for (k = 0; k < n - 1; k++)
    {
        somaN += L[n - 1][k] * U[k][n - 1];
    }
    //ultimo elemento da diagonal
    U[n - 1][n - 1] = ((*A)[n - 1][n - 1] - somaN) / L[n - 1][n - 1];

    //Forward Substitution
    Y[0] = (*b)[0] / L[0][0];
    for (i = 1; i < n; i++)
    {
        double somaFS;
        somaFS = 0;
        for (j = 0; j < i; j++)
        {
            somaFS += L[i][j] * Y[j];
        }
        Y[i] = (1 / L[i][i]) * ((*b)[i] - somaFS);
    }
    //Backwards Substitution
    result[n - 1] = Y[n - 1] / U[n - 1][n - 1];
    for (i = n - 2; i >= 0; i--)
    {
        double somaBS;
        somaBS = 0;
        for (k = i + 1; k < n; k++)
        {
            somaBS += U[i][k] * result[k];
        }
        result[i] = (1 / U[i][i]) * (Y[i] - somaBS);
    }

    return (result);
}
double *Cholesky(double ***A, double **b, int n)
{
    double *result;
    result = alocaVetor(n);
    int i, j, k;
    double **L;
    double **U;
    U = alocaMatriz(n, n);
    L = alocaMatriz(n, n);
    double *Y;
    Y = alocaVetor(n);

    L[0][0] = sqrt((*A)[0][0]);
    for (j = 1; j < n; j++)
    {
        L[j][0] = (*A)[j][0] / L[0][0];
    }
    for (i = 1; i < n - 1; i++)
    {
        double somaD = 0;
        for (k = 0; k < i; k++)
        {
            somaD = somaD + L[i][k] * L[i][k];
        }
        L[i][i] = sqrt((*A)[i][i] - somaD);
        for (j = i + 1; j < n; j++)
        {
            double somaJ = 0;
            for (k = 0; k < i; k++)
            {
                somaJ = somaJ + L[j][k] * L[i][k];
            }
            L[j][i] = ((*A)[j][i] - somaJ) / L[i][i];
        }
    }
    double somaF = 0;
    for (k = 0; k < n - 1; k++)
    {
        somaF = somaF + L[n - 1][k] * L[n-1][k];
    }
    L[n - 1][n - 1] = sqrt((*A)[n - 1][n - 1] - somaF);
    for (i=0;i<n;i++){
        for (j=0;j<n;j++){
            U[j][i] = L[i][j];
        }
    }
    //Forward Substitution
    Y[0] = (*b)[0] / L[0][0];
    for (i = 1; i < n; i++)
    {
        double somaFS;
        somaFS = 0;
        for (j = 0; j < i; j++)
        {
            somaFS += L[i][j] * Y[j];
        }
        Y[i] = (1 / L[i][i]) * ((*b)[i] - somaFS);
    }
    //Backwards Substitution
    result[n - 1] = Y[n - 1] / U[n - 1][n - 1];
    for (i = n - 2; i >= 0; i--)
    {
        double somaBS;
        somaBS = 0;
        for (k = i + 1; k < n; k++)
        {
            somaBS += U[i][k] * result[k];
        }
        result[i] = (1 / U[i][i]) * (Y[i] - somaBS);
    }

    return (result);
}
double *LUPivot(double ***A, double **b, int n)
{
    double *result = NULL;
    result = alocaVetor(n);
    int i, k;
    double *aux_troca = NULL;
    aux_troca = alocaVetor(n);
    double aux_troca_b = 0;
    double max_value = 0;
    int index_linha = 0;

    //Processo de pivoteamento
    for (k = 0; k < n; k++)
    {
        max_value = (*A)[k][k];
        index_linha = k;
        for (i = k + 1; i < n; i++)
        {
            if ((*A)[i][k] > max_value)
            {
                index_linha = i;
                max_value = (*A)[i][k];
            }
        }

        if (index_linha != k)
        {
            //troca a linha de A
            aux_troca = (*A)[k];
            (*A)[k] = (*A)[index_linha];
            (*A)[index_linha] = aux_troca;

            //trocar o vetor b
            aux_troca_b = (*b)[k];
            (*b)[k] = (*b)[index_linha];
            (*b)[k] = aux_troca_b;
        }
    }
    result = LUFactorization(A, b, n);
    return (result);
}
double *GramSchmidt(double ***A, double **b, int n)
{
    double *result;
    result = alocaVetor(n);
    return (result);
}
double *Householder(double ***A, double **b, int n)
{
    double *result;
    result = alocaVetor(n);
    return (result);
}
double *Givens(double ***A, double **b, int n)
{
    double *result;
    result = alocaVetor(n);
    return (result);
}
