#pragma once
#ifndef MATMULT_H_
#define MATMULT_H_

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
void matrix_multiply(double* A[], int k,double B[], double* C)
{

	int rowA = k;

	int colA = k;

	int rowB = k;
	

	for (int i = 0; i < rowA; ++i)
	{
		for (int k = 0; k < colA; ++k)
		{
			C[i] += A[i][k] * B[k];
		}
	}
}
#endif