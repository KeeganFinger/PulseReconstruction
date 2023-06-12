//
// Created by finge on 6/5/2023.
//

#include "physics_functions.h"

template<typename T>
double sgn(T x) {
    return (T(0) < x) - (x < T(0));
}


int Factorial(int n) {
    int fact = 1;
    for (int i = 1; i <= n; i++) {
        fact *= i;
    }
    return fact;
}

complex<double> ErrorFunction(complex<double> z, double a) {
    if (abs(z + 2i) < 1.1) {
        z += 2i;
        return exp(a) *
               (-18.564802i + 12.326185 * z + 13.9957261i * z * z - 7.0912286 * z * z * z - 2.1539997i * z * z * z * z +
                0.80057144 * z * z * z * z * z) /
               (1.0 - 2.6545518i * z - 2.9260194 * z * z + 1.665267i * z * z * z + 0.48178231 * z * z * z * z -
                0.054052386i * z * z * z * z * z);
    } else if (abs(z - 2i) < 1.1) {
        z -= 2i;
        return exp(a) *
               (18.564802i + 12.326185 * z - 13.9957261i * z * z - 7.0912286 * z * z * z + 2.1539997i * z * z * z * z +
                0.80057144 * z * z * z * z * z) /
               (1.0 + 2.6545518i * z - 2.9260194 * z * z - 1.665267i * z * z * z + 0.48178231 * z * z * z * z +
                0.054052386i * z * z * z * z * z);
    } else if (abs(z) < 1.4) {
        return sgn(z.real()) * sqrt(exp(2.0 * a) - exp(2.0 * a -
                                                       (1.2732395 + 0.14001229 * z * z) / (1.0 + 0.14001229 * z * z) *
                                                       z * z));
    } else {
        return sgn(z.real()) * exp(a) - exp(a - z * z) * z / (sqrt(M_PI) * (z * z + 0.5 / (1.0 + 1.0 / (1.5 + z * z))));
    }
}

double Wigner3J(int j1, int j2, int j3, int m1, int m2, int m3) {
    //Compute the Wigner 3j symbol using the Racah formula
    int t1 = j2 - m1 - j3;
    int t2 = j1 + m2 - j3;
    int t3 = j1 + j2 - j3;
    int t4 = j1 - m1;
    int t5 = j2 + m2;

    int tmin = max(0, max(t1, t2));
    int tmax = min(t3, min(t4, t5));

    double wigner = 0;
    for (int t = tmin; t <= tmax; t++) {
        wigner += pow(-1, t) /
                  (Factorial(t) * Factorial(t - t1) * Factorial(t - t2) * Factorial(t3 - t) * Factorial(t4 - t) *
                   Factorial(t5 - t));
    }
    return wigner * pow(-1, j1 - j2 - m3) *
           sqrt(Factorial(j1 + j2 - j3) * Factorial(j1 - j2 + j3) * Factorial(-j1 + j2 + j3) /
                Factorial(j1 + j2 + j3 + 1) * Factorial(j1 + m1) * Factorial(j1 - m1) * Factorial(j2 + m2) *
                Factorial(j2 - m2) * Factorial(j3 + m3) * Factorial(j3 - m3));
}