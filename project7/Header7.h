#pragma once
#ifndef HEADER7_H_

#define HEADER7_H_

#include<iostream>
#include<vector>
#include<math.h>
#include<time.h>
#include<stdlib.h>
#include<array>
#include<numeric>
#include<cmath>
#include<fstream>
#include <chrono>
#include <random>

using namespace std;
double multiply_vector(double* X, double* Y, int length) {
	double result=0;
	for (int i = 0; i < length; i++) {
		result += X[i] * Y[i];
	}
	return result;
}
double* derive_F(double* F, double** A, double* B, int size, int m) {
	double* f = new double[size]; 
	
		for (int i = 0; i < size; i++) {
			f[i] = multiply_vector(A[i], F, size) + B[i];
		}
		return f;
	
}
#endif