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
    double *result = NULL;
    result = alocaVetor(n);
    double *Ax = NULL;
    double *r0 = NULL;
    double *p1 = NULL;
    double *Ar0 = NULL;
    double *r1 = NULL;
    double *p_new = NULL;
    Ax = alocaVetor(n);
    r0 = alocaVetor(n);
    p1 = alocaVetor(n);
    p_new = alocaVetor(n);
    Ar0 = alocaVetor(n);
    r1 = alocaVetor(n);
    double *x0 = NULL;
    x0 = alocaVetor(n);
    double *x_atual;
    x_atual = alocaVetor(n);
    double q1 = 0.0;
    double alfa = 0.0;
    int iter = 0;
    int i_max = 100;
    double maxX0;
    double max_Atual;
    double erro = 1000;
    double soma;

    int i, j, k;
    //seguindo nomenclatura da pagina 178 do NEIDE
    //Ax0
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            Ax[i] += (*A)[i][j] * x0[j];
        }
    }
    //r0 = Ax + b
    for (i = 0; i < n; i++)
    {
        r0[i] += Ax[i] + (*b)[i];
    }
    //A * r0
    for (i = 0; i < n; i++)
    {
        soma = 0;
        for (j = 0; j < n; j++)
        {
            soma += (*A)[i][j] * r0[j];
        }
        Ar0[i] = soma;
    }

    for (i = 0; i < n; i++)
    {
        p1[i] = -r0[i];
    }
    q1 = (dotProduct(&r0, &r0, n) / dotProduct(&Ar0, &r0, n));

    //calculo de x_atual (v k+1)
    for (i = 0; i < n; i++)
    {
        x_atual[i] = x0[i] + q1 * p1[i];
        x0[i] = x_atual[i];
    }
    //calculo de r1
    for (i = 0; i < n; i++)
    {
        soma = 0;
        for (j = 0; j < n; j++)
        {
            soma += (*A)[i][j] * p1[j];
        }
        r1[i] = r0[i] + (q1 * soma);
    }

    while (iter < i_max)
    {
        //for k>2
        alfa = (dotProduct(&r1, &r1, n) / dotProduct(&r0, &r0, n));
        for (i = 0; i < n; i++)
        {
            p_new[i] = -r1[i] + (alfa * p1[i]);
        }
        //atualiza Ar0 como Ap_new

        for (i = 0; i < n; i++)
        {
            soma = 0;
            for (j = 0; j < n; j++)
            {
                soma += (*A)[i][j] * p_new[j];
            }
            Ar0[i] = soma;
        }
        //q1 new
        q1 = (dotProduct(&r1, &r1, n) / dotProduct(&Ar0, &p_new, n));

        for (i = 0; i < n; i++)
        {
            x_atual[i] = x0[i] + q1 * p_new[i];
        }
        
        for (i = 0; i < n; i++)
        {
            soma = 0;
            for (j = 0; j < n; j++)
            {
                soma += (*A)[i][j] * p_new[j];
            }
            r1[i] = r0[i] + q1 * soma;
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
        //atualiza a iteracao
        for (i = 0; i < n; i++)
        {
            x0[i] = x_atual[i];
            r0[i] = r1[i];
            p1[i] = p_new[i];
        }
        iter += 1;
        printf("iter: %d\n", iter);
        printf("erro: %f\n", erro);
    }

    result = x_atual;
    return (result);
}
double *GradPreCondicionados(double ***A, double **b, int n)
{
    double *result;
    result = (double *)malloc(n * sizeof(double));
    return (result);
}