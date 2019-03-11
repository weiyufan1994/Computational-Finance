#ifndef INVERSE_H_
#define INVERSE_H_

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
void inverse_matrix(double* A[], int n, double* C[])
{
	int i, j, k, m = 2 * n;
	double mik, temp;
	double **a = new double*[n];
	double **B = new double*[n];

	for (i = 0; i<n; i++)
	{
		a[i] = new double[2 * n];
		B[i] = new double[n];
	}

	for (i = 0; i<n; i++)
	{
		for (j = 0; j<n; j++)
		{
			if (i == j)
				B[i][j] = 1.0;
			else
				B[i][j] = 0.0;
		}
	}        //initialize B=E

	
	for (i = 0; i<n; i++) {
		for (j = 0; j<n; j++) {
			a[i][j] = A[i][j];  //����A��a������ı�A��ֵ
		}
	}

	for (i = 0; i<n; i++)
		for (j = n; j<m; j++)
			a[i][j] = B[i][j - n];  //����B��a���������

	for (k = 1; k <= n - 1; k++)
	{
		for (i = k + 1; i <= n; i++)
		{
			mik = a[i - 1][k - 1] / a[k - 1][k - 1];
			for (j = k + 1; j <= m; j++)
			{
				a[i - 1][j - 1] -= mik*a[k - 1][j - 1];
			}
		}
	}        //˳���˹��ȥ�������½�Ϊ��

	for (i = 1; i <= n; i++)
	{
		temp = a[i - 1][i - 1];
		for (j = 1; j <= m; j++)
		{
			a[i - 1][j - 1] /= temp;
		}
	}        //��һ��

	for (k = n - 1; k >= 1; k--)
	{
		for (i = k; i >= 1; i--)
		{
			mik = a[i - 1][k];
			for (j = k + 1; j <= m; j++)
			{
				a[i - 1][j - 1] -= mik*a[k][j - 1];
			}
		}
	}        //�����˹��ȥ��������������Ϊ��λ����

	for (i = 0; i<n; i++)
		for (j = 0; j<n; j++)
			B[i][j] = a[i][j + n];  //ȡ��������

	for (i = 0; i<n; i++)
		for (j = 0; j<n; j++)
			if (fabs(B[i][j])<0.0001)
				B[i][j] = 0.0;
	
	for (i = 0; i<n; i++) {
		for (j = 0; j<n; j++) {
			C[i][j] = B[i][j];
		}
	}
}

#endif