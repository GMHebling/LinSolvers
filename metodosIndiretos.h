/* 
 * File:   metodosIndiretos.h
 * Author: Gustavo
 *
 * Created on 31 de Outubro de 2018, 10:02
 */

#ifndef METODOSINDIRETOS_H
#define	METODOSINDIRETOS_H

double *Jacobi(double ***A, double **b, int n);
double *GaussSeidel(double ***A, double **b, int n);
double *GradientesConjugados(double ***A, double **b, int n);
double *GradPreCondicionados(double ***A, double **b, int n);

#endif	/* METODOSINDIRETOS_H */

