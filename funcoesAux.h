/* 
 * File:   funcoesAux.h
 * Author: Julio
 *
 * Created on 30 de Outubro de 2018, 15:08
 */

#ifndef FUNCOESAUX_H
#define	FUNCOESAUX_H

double **alocaMatriz(int m, int n);
double *alocaVetor(int m);
void imprimeMatriz(double **A,int m,int n);
void imprimeVetor(double *A,int m);
double dotProduct(double **A, double **B, int n);
double *MultVetorMatriz(double ***A, double **B, int n);
double **MultMatriz(double ***A, double ***B, int n);


void leituraSistemaLinear(double ***A, double **b, int *n);

#endif	/* FUNCOESGUSTAVO_H */

