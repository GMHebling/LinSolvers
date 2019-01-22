/* 
 * File:   main.c
 * Author: Julio
 *
 * Created on 30 de Outubro de 2018, 15:07
 */

#include <stdio.h>
#include <stdlib.h>
#include "funcoesAux.h"
#include "metodosDiretos.h"
#include "metodosIndiretos.h"

/*
 * 
 */
int main(int argc, char** argv) {
    double **A = NULL;
    double *b = NULL; 
    double *x = NULL;
    int n;
    
    leituraSistemaLinear(&A, &b, &n);
    x = alocaVetor(n);

    //x = GradientesConjugados(&A, &b, n);
    x = GradPreCondicionados(&A, &b, n);
    //x = LUFactorization(&A, &b, n);
    //x = Cholesky(&A, &b, n);
    //x = Jacobi(&A, &b, n);
    //x = GaussSeidel(&A, &b, n)    
    
    imprimeVetor(x,n);
    getchar();
    
    return (EXIT_SUCCESS);
}

