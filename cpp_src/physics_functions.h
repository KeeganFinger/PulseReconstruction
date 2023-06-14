//
// Created by finge on 6/5/2023.
//

#ifndef PULSERECONSTRUCTION_PHYSICS_FUNCTIONS_H
#define PULSERECONSTRUCTION_PHYSICS_FUNCTIONS_H

#include <cmath>
#include <complex>
#include <vector>

#define I std::complex<double>(0,1)

using namespace std;

template<typename T>
double sgn(T x);

int Factorial(int n);

complex<double> ErrorFunction (complex<double> z, double a);

double Wigner3J(int j1, int j2, int j3, int m1, int m2, int m3);

#endif //PULSERECONSTRUCTION_PHYSICS_FUNCTIONS_H
