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
    double *b = NULL, *x = NULL;
    int n;
    
    leituraSistemaLinear(&A, &b, &n);
    x = alocaVetor(n);
    x = GradientesConjugados(&A, &b, n);    
    //imprimeMatriz(A,n,n);
    //imprimeVetor(b,n);
    imprimeVetor(x,n);
    
    return (EXIT_SUCCESS);
}

