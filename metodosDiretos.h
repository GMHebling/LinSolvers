/* 
 * File:   metodosDiretos.h
 * Author: Gustavo
 *
 * Created on 31 de Outubro de 2018, 09:45
 */

#ifndef METODOSDIRETOS_H
#define	METODOSDIRETOS_H

double *LUFactorization(double ***A, double **b, int n);
double *LUPivot(double ***A, double **b, int n);
double *Cholesky(double ***A, double **b, int n);
double *GramSchmidt(double ***A, double **b, int n);
double *Houserholder(double ***A, double **b, int n);
double *Givens(double ***A, double **b, int n);

#endif	/* METODOSDIRETOS_H */

