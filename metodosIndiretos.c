#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "metodosIndiretos.h"
#include "funcoesAux.h"
//-----------------------------------------------------------------------------------------
// METODOS INDIRETOS DE SOLUCAO DE SISTEMAS LINEARES
//-----------------------------------------------------------------------------------------

double *Jacobi(double ***A, double **b, int n)
{
    int max_iter = 10;
    double *result = NULL;
    result = alocaVetor(n);
    int i, j, k;
    double soma;
    double *x0 = NULL;
    double *x_atual = NULL;
    x0 = alocaVetor(n);
    x_atual = alocaVetor(n);
    double erro;
    //double tol = 0.001;
    double maxX0;
    double max_Atual;
    int iter = 0;

    for (i = 0; i < n; i++)
    {
        x0[i] = 0.1;
    }

    //verificadao da diagonal dominante
    int ctAux = 0;
    double maxN;
    for (i = 0; i < n; i++)
    {
        double maxDiag;
        maxDiag = fabs((*A)[i][i]);
        for (j = 0; j < n; j++)
        {
            if (i != j)
            {
                maxN = (*A)[i][j];
                if (maxN > maxDiag)
                {
                    ctAux += 1;
                }
            }
        }
    }
    if (ctAux == 0)
    {
        printf("diagonal dominante");
    }
    else
    {
        printf("nao diagonal dominante");
    }
    while (iter < max_iter)
    {
        for (i = 0; i < n; i++)
        {
            // calcula a soma da formula de Jacobi
            soma = 0;
            for (j = 0; j < n; j++)
            {
                if (i != j)
                {
                    soma += (*A)[i][j] * x0[j];
                }
            }
            x_atual[i] = (1 / (*A)[i][i]) * ((*b)[i] - soma);
        }
        // verifica criterio de parada
        maxX0 = fabs(x0[0]);
        max_Atual = fabs(x_atual[0]);
        for (k = 0; k < n; k++)
        {
            if (fabs(x_atual[k]) > max_Atual)
            {
                max_Atual = fabs(x_atual[k]);
            }
            if (fabs(x0[k]) > maxX0)
            {
                maxX0 = fabs(x0[k]);
            }
        }
        erro = (max_Atual - maxX0) / max_Atual;

        for (i = 0; i < n; i++)
        {
            x0[i] = x_atual[i];
        }
        iter = iter + 1;
    }

    return x_atual;
}

double *GaussSeidel(double ***A, double **b, int n)
{
    int max_iter = 20;
    double *result = NULL;
    result = alocaVetor(n);
    int i, j, k;
    double soma, soma2;
    double *x0 = NULL;
    double *x_atual = NULL;
    x0 = alocaVetor(n);
    x_atual = alocaVetor(n);
    double erro;
    double tol = 0.000001;
    double maxX0;
    double max_Atual;
    int iter = 0;

    //inicializa vetores
    for (k = 0; k < n; k++)
    {
        x0[k] = 0;
        x_atual[k] = 0;
    }
    while (iter < max_iter)
    {
        for (i = 0; i < n; i++)
        {
            if (i == 0)
            {
                soma = 0;
                for (j = 0; j < n; j++)
                {
                    if (j != i)
                    {
                        soma = soma + (-1 * (*A)[i][j] * x0[j]);
                    }
                }
                x_atual[i] = (1 / (*A)[i][i]) * (soma + (*b)[i]);
            }
            else
            {
                soma = 0;
                for (j = 0; j < n; j++)
                {
                    if (j != i)
                    {
                        soma = soma + (-1 * (*A)[i][j] * x0[j]);
                    }
                }
                soma2 = 0;
                for (j = 0; j < i; j++)
                {
                    soma2 = soma2 + (-1 * (*A)[i][j] * x_atual[j]);
                }
                x_atual[i] = (1 / (*A)[i][i]) * (soma2 + soma + (*b)[i]);
            }
        }
        // verifica criterio de parada
        maxX0 = fabs(x0[0]);
        max_Atual = fabs(x_atual[0]);
        for (k = 0; k < n; k++)
        {
            if (fabs(x_atual[k]) > max_Atual)
            {
                max_Atual = fabs(x_atual[k]);
            }
            if (fabs(x0[k]) > maxX0)
            {
                maxX0 = fabs(x0[k]);
            }
        }
        erro = (max_Atual - maxX0) / max_Atual;
        if (erro < tol)
        {
            result = x_atual;
            return (result);
        }
        else
        {
            x0 = x_atual;
            iter += 1;
        }
    }
    return x_atual;
}
double *GradientesConjugados(double ***A, double **b, int n)
{
    double *result;
    result = alocaVetor(n);
    double *x0;
    x0 = alocaVetor(n);
    double *x_atual;
    x_atual = alocaVetor(n);
    double passo = 0.0;
    return (result);
}
double *GradPreCondicionados(double ***A, double **b, int n)
{
    double *result;
    result = (double *)malloc(n * sizeof(double));
    return (result);
}