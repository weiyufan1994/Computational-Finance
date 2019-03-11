#pragma once
#pragma once
#ifndef HEADER5_H_
#define HEADER5_H_

#include"Header.h"
#include"inverse.h"
#include"MatMult.h"
#include"Polynomial.h"
#include<iostream>
#include<vector>
#include<math.h>
#include<time.h>
#include<stdlib.h>
#include<array>
#include<numeric>
#include<cmath>
#include<fstream>
#include<chrono>
#include<random>
using namespace std;
void transpose(double **A, int N) {
	double temp;
	for (int i = 0; i < N; i++) {
		for (int j = i; j < N; j++) {
			temp = A[i][j];
			A[i][j] = A[j][i];
			A[j][i] = temp;
		}
	}
}
//determinant
double determinant(double **A, int n) {
	double det = 0;
	double **array;
	array = new double *[n];
	for (int i = 0; i < n; i++) {
		array[i] = new double[n];
	}
	if (n == 2)
		return ((A[0][0] * A[1][1]) - (A[1][0] * A[0][1]));
	else {
		for (int x = 0; x < n; x++) {
			int subi = 0;
			for (int i = 1; i < n; i++) {
				int subj = 0;
				for (int j = 0; j < n; j++) {
					if (j == x)
						continue;
					array[subi][subj] = A[i][j];
					subj++;
				}
				subi++;
			}
			det = det + (pow(-1, x) * A[0][x] * determinant(array, n - 1));
		}
	}
	return det;
}
//Adjoint of a Matrix
double **adjoint(double **A, int N) {
	double **adjointA;
	adjointA = new double *[N];
	for (int i = 0; i < N; i++) {
		adjointA[i] = new double[N];
	}

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double **array;
			array = new double *[N - 1];
			for (int m = 0; m < N - 1; m++) {
				array[m] = new double[N - 1];
			}
			int ii = 0;
			for (int m = 0; m < N; m++) {
				int jj = 0;
				if (i == m) continue;
				for (int n = 0; n < N; n++) {
					if (j == n) continue;
					array[ii][jj] = A[m][n];
					jj++;
				}
				ii++;
			}
			adjointA[i][j] = determinant(array, N - 1) * pow(-1, i + j);
		}
	}
	transpose(adjointA, N);
	return adjointA;
}
//Calculate coefficients of Linear Regression, coef = inv(A)*b
double sumprod(double *x, double *y, int N) {
	double sum = 0;
	for (int i = 0; i < N; i++) {
		sum += x[i] * y[i];
	}
	return sum;
}
double *coef(double **A, double *b, int N) {
	double *array = new double[N];
	double det = determinant(A, N);
	double **adjA = adjoint(A, N);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			adjA[i][j] /= det; // A inverse
		}
	}

	for (int i = 0; i < N; i++) {
		array[i] = sumprod(adjA[i], b, N);
	}
	return array;
}
double ysL(double a[],int k, double s) {
	double sum = 0;
	Laguerre l1(1), l2(2), l3(3), l4(4);
	switch (k)
	{
	case 2:
		sum += a[0] * l1.poly(s) + a[1] * l2.poly(s); break;
	case 3:
		sum += a[0] * l1.poly(s) + a[1] * l2.poly(s) + a[2] * l3.poly(s); break;
	case 4:
		sum += a[0] * l1.poly(s) + a[1] * l2.poly(s) + a[2] * l3.poly(s) + a[3] * l4.poly(s); break;
	default:
		break;
	}
	return sum;
}
double ysM(double a[], int k, double s) {
	 double sum = 0;
	 Monomials l1(1), l2(2),l3(3),l4(4);
	 switch (k)
	 {
	 case 2:
		 sum += a[0] * l1.poly(s) + a[1] * l2.poly(s); break;
	 case 3:
		 sum += a[0] * l1.poly(s) + a[1] * l2.poly(s)+a[2] * l3.poly(s); break;
	 case 4:
		 sum += a[0] * l1.poly(s) + a[1] * l2.poly(s) + a[2] * l3.poly(s) + a[3] * l4.poly(s); break;
	 default:
		 break;
	 }
	return sum;
}
double ysH(double a[], int k, double s) {
	double sum = 0;
	Hermite l1(1), l2(2), l3(3), l4(4);
	switch (k)
	{
	case 2:
		sum += a[0] * l1.poly(s) + a[1] * l2.poly(s); break;
	case 3:
		sum += a[0] * l1.poly(s) + a[1] * l2.poly(s) + a[2] * l3.poly(s); break;
	case 4:
		sum += a[0] * l1.poly(s) + a[1] * l2.poly(s) + a[2] * l3.poly(s) + a[3] * l4.poly(s); break;
	default:
		break;
	}
	return sum;
}
void regressMonomial(vector<double> x, vector<double> y, double* ECV, int k) {
	//Laguerre polynomials
	Monomials l1(1),l2(2),l3(3),l4(4);//create a vector to store k Laguerre polynomials,
	double* f1 = new double[x.size()]; double* f2 = new double[x.size()]; double* f3 = new double[x.size()];
	double* f4 = new double[x.size()];
	for (int i = 0; i < x.size(); i++) {
		f1[i] = l1.poly(x[i]);
		f2[i] = (l2.poly(x[i]));
		f3[i] = (l3.poly(x[i]));
		f4[i] = (l4.poly(x[i]));
	}
	double *F[4]; F[0] = f1; F[1] = f2; F[2] = f3; F[3] = f4;

	double* b = new double[k];

	for (int i = 0; i < k; i++) {
		switch (k) {
		case 2:
			b[0] = vector_dot_prod(y, f1);
			b[1] = vector_dot_prod(y, f2); break;
		case 3:
			b[0] = vector_dot_prod(y, f1);
			b[1] = vector_dot_prod(y, f2); b[2] = vector_dot_prod(y, f3);
			break;
		case 4:
			b[0] = vector_dot_prod(y, f1);
			b[1] = vector_dot_prod(y, f2); b[2] = vector_dot_prod(y, f3); b[3] = vector_dot_prod(y, f4);
			break;

		}
	}

	double** A; double** A_inv; A = new double*[k]; A_inv = new double*[k];


	for (int i = 0; i < k; i++) {
		A[i] = new double[k];
		A_inv[i] = new double[k];
	}

	for (int i = 0; i < k; i++) {
		for (int j = 0; j < k; j++) {
			A[i][j] = sumprod(F[i], F[j],x.size());
		}
	}

	//	inverse_matrix(A, k, A_inv);

	double* a = new double[k] {0};

	//matrix_multiply(A_inv, k, b, a);//calculate the a;

	a = coef(A, b, k);

	for (int i = 0; i < k; i++) {
		delete[] A[i];
		delete[] A_inv[i];
	}
	delete[] A; delete[] A_inv;  delete[] b;

	//calculate the expected continous value

	for (int i = 0; i < x.size(); i++) {
		ECV[i]=ysM(a, k, x[i]);
	}
	delete[] a; delete[]f1; delete[]f2; delete[]f3; delete[]f4;
}
void regressLaguerre(vector<double> x, vector<double> y, double* ECV, int k) {
	//Laguerre polynomials
	Laguerre l1(1), l2(2), l3(3), l4(4);//create a vector to store k Laguerre polynomials,
	double* f1 = new double[x.size()]; double* f2 = new double[x.size()]; double* f3 = new double[x.size()];
	double* f4 = new double[x.size()];
	for (int i = 0; i < x.size(); i++) {
		f1[i]=l1.poly(x[i]);
		f2[i]=(l2.poly(x[i]));
		f3[i]=(l3.poly(x[i]));
		f4[i]=(l4.poly(x[i]));
	}
	double *F[4]; F[0] = f1; F[1] = f2; F[2] = f3; F[3] = f4;

	double* b = new double[k];

	for (int i = 0; i < k; i++) {
		switch (k) {
		case 2:
			b[0] = vector_dot_prod(y, f1);
			b[1] = vector_dot_prod(y, f2); break;
		case 3:
			b[0] = vector_dot_prod(y, f1);
			b[1] = vector_dot_prod(y, f2); b[2] = vector_dot_prod(y, f3);
			break;
		case 4:
			b[0] = vector_dot_prod(y, f1);
			b[1] = vector_dot_prod(y, f2); b[2] = vector_dot_prod(y, f3); b[3] = vector_dot_prod(y, f4);
			break;
		
		}
	}

	double** A; double** A_inv; A = new double*[k]; A_inv = new double*[k];


	for (int i = 0; i < k; i++) {
		A[i] = new double[k];
		A_inv[i] = new double[k];
	}

	for (int i = 0; i < k; i++) {
		for (int j = 0; j < k; j++) {
			A[i][j] = sumprod(F[i], F[j], x.size());
		}
	}

	//	inverse_matrix(A, k, A_inv);

	double* a = new double[k] {0};

	//matrix_multiply(A_inv, k, b, a);//calculate the a;

	a = coef(A, b, k);

	for (int i = 0; i < k; i++) {
		delete[] A[i];
		delete[] A_inv[i];
	}
	delete[] A; delete[] A_inv;  delete[] b;

	//calculate the expected continous value

	for (int i = 0; i < x.size(); i++) {
		ECV[i]=ysM(a, k, x[i]);
	}
	delete[] a; delete[]f1; delete[]f2; delete[]f3;
}
void regressHermite(vector<double> x, vector<double> y, double* ECV, int k) {
	//Laguerre polynomials
	Hermite l1(1), l2(2), l3(3), l4(4);//create a vector to store k Laguerre polynomials,
	double* f1 = new double[x.size()]; double* f2 = new double[x.size()]; double* f3 = new double[x.size()];
	double* f4 = new double[x.size()];
	for (int i = 0; i < x.size(); i++) {
		f1[i] = l1.poly(x[i]);
		f2[i] = (l2.poly(x[i]));
		f3[i] = (l3.poly(x[i]));
		f4[i] = (l4.poly(x[i]));
	}

	double *F[4]; F[0] = f1; F[1] = f2; F[2] = f3; F[3] = f4;
	double* b = new double[k];

	for (int i = 0; i < k; i++) {
		switch (k) {
		case 2:
			b[0] = vector_dot_prod(y, f1);
			b[1] = vector_dot_prod(y, f2); break;
		case 3:
			b[0] = vector_dot_prod(y, f1);
			b[1] = vector_dot_prod(y, f2); b[2] = vector_dot_prod(y, f3);
			break;
		case 4:
			b[0] = vector_dot_prod(y, f1);
			b[1] = vector_dot_prod(y, f2); b[2] = vector_dot_prod(y, f3); b[3] = vector_dot_prod(y, f4);
			break;

		}
	}

	double** A; double** A_inv; A = new double*[k]; A_inv = new double*[k];


	for (int i = 0; i < k; i++) {
		A[i] = new double[k];
		A_inv[i] = new double[k];
	}

	for (int i = 0; i < k; i++) {
		for (int j = 0; j < k; j++) {
			A[i][j] = sumprod(F[i], F[j], x.size());
		}
	}

	//	inverse_matrix(A, k, A_inv);

	double* a = new double[k] {0};

	//matrix_multiply(A_inv, k, b, a);//calculate the a;

	a = coef(A, b, k);

	for (int i = 0; i < k; i++) {
		delete[] A[i];
		delete[] A_inv[i];
	}
	delete[] A; delete[] A_inv;  delete[] b;

	//calculate the expected continous value

	for (int i = 0; i < x.size(); i++) {
		ECV[i]=ysM(a, k, x[i]);
	}
	delete[] a; delete[]f1; delete[]f2; delete[]f3; delete[]f4;
}
void generate_simulations(vector<vector<int>> &x,double** simulation,double** option_value, double S0 = 36, double sigma = 0.2, double r = 0.06, int n = 100000, double X = 40, int m = 100, double T = 0.5) {
	//first simulate 100000 paths of St
	vector<int> xi(100, 0); x.push_back(xi); xi.clear();
	for (int i = 0; i < n; i++) {
		simulation[i][0] = S0;
	}
	double delta = T / m;
	for (int i = 1; i <= m; i = i + 1) {
		vector<double> uniform, normal; generate_uniform_distributed_sequence(uniform, 2+(n / 2), i);
		generate_standard_normal_distributed_sequence(normal, uniform);
		for (int j = 0; j < n;j=j+2){
			
			simulation[j][i] = simulation[j][i - 1]+simulation[j][i - 1] * r*delta + simulation[j][i - 1] * sigma*sqrt(delta)*normal[j/2];
			option_value[j][i] = max(X-simulation[j][i],0.0);
			simulation[j+1][i] = simulation[j][i - 1]+simulation[j+1][i - 1] * r*delta + simulation[j+1][i - 1] * sigma*sqrt(delta)*(-normal[j/2]);
			option_value[j+1][i] = max(X-simulation[j+1][i], 0.0);
			if (option_value[j][i]>0)
				xi.push_back(j);
			if (option_value[j + 1][i]>0)
				xi.push_back(j + 1);
		}
		x.push_back(xi);//start from i=1, record the coordinates of S>0;
		xi.clear();
	}
}
void question1L(string p, vector<vector<int>> x, double** simulation, double** option_value, int n = 100000, int m = 100, int k = 3, double T = 0.5) {
	double r = 0.06;
	//The initial parameter
	double delta = T / m;
	double discount;
	discount = exp(-r*delta);
	/*The last period*/
	//When t=T, check if the option will be exercised
	vector<double> xT, yT;
	// create the x and y for regression
	for (int i = 0; i < x[m-1].size(); i++) {
		xT.push_back(simulation[x[m-1][i]][m-1]);//store those exercised paths in a vector
		yT.push_back(discount*option_value[x[m - 1][i]][m]);//store the discounted option value to the y
	}

	double* ECV=new double[xT.size()];
	if (p == "Monomial")
		regressMonomial(xT, yT, ECV, k);
	else if (p == "Laguerre")
		regressLaguerre(xT, yT, ECV, k);
	else if (p == "Hermite")
		regressHermite(xT, yT, ECV, k);
	else
		cout << "re-enter polynomials" << endl;
	vector<int> index_horizontal(n, m);
	//replace the EV with ECV, if not exercise
	for (int i = 0; i < xT.size(); i++) {
		if (option_value[x[m - 1][i]][m - 1]>ECV[i]) {
			index_horizontal[x[m - 1][i]] = m - 1;
		}
	}

	for (int i = m-1; i > 1; i--) {
		 xT.clear(), yT.clear(); 
		// create the x and y for regression
		for (int j = 0; j < x[i-1].size(); j++) {
				discount = (-r)*delta*(index_horizontal[x[i - 1][j]] + 1 - i);
				double yi = option_value[x[i - 1][j]][index_horizontal[x[i - 1][j]]] * exp(discount);// Y value
				xT.push_back(simulation[x[i - 1][j]][i-1]);//store those exercised paths in a vector
				yT.push_back(yi);//store the discounted option value to the y
			}
		int s = xT.size();

		double* ECV=new double[s];
		regressMonomial(xT, yT, ECV, k);		//replace the EV with ECV, if not exercise
		for (int j = 0; j < xT.size(); j++) {
			if (option_value[x[i - 1][j]][i - 1]>ECV[j]) {
				index_horizontal[x[i-1][j]] = i - 1;
			}
		}
	}
	int count = 0; double sum = 0;
	for (int i = 0; i < n; i++) {
		option_value[i][0] = option_value[i][index_horizontal[i]] * exp(-r*delta*(index_horizontal[i]));
		sum += option_value[i][0];
		//   cout << option_value[i][index_horizontal[i]] << " " << pow(discount, m + index_horizontal[i]) << endl;
		//cout << option_value[i][m]<<endl;
	}
	//		cout << count<<endl;
	cout << sum / n;
}


void question2a(double S0 = 65, double sigma = 0.2, double r = 0.06, int n = 10000, double X = 60, int m = 100, double T =1) {

	/*Question 1*/

	//The initial parameters
	double** simulation;
	double** option_value;
	double** option_price;

	double delta = T / m;
	double discount;
	discount = exp(-r*delta);
	/*this vector is used to store the paths which have positive ECV*/
	vector<int>index_vertical;
	/*this vector is used to store the time of exercising option,starting from the end*/
	vector<int>index_horizontal(n, m);
	//first simulate 100000 paths of St
	simulation = new double*[n]; option_value = new double*[n]; option_price = new double*[n];
	//create a n*m matrix to store the simulated paths.
	for (int i = 0; i < n; i++) {
		simulation[i] = new double[m + 1];
		option_value[i] = new double[m + 1];// option value is the exercise value
		option_price[i] = new double[m + 1];// option price is the max of exercise value and ECV
	}
	double sum = 0; int time = 0.2 / delta;
	for (int i = 0; i < n; i = i + 2) {
		double Sx = S0; double Sy = S0;
		vector<double> uniform, normal; int seed = i % 123465 + abs(generate_normal_number(1000, 300));
		simulation[i][0] = S0; simulation[i + 1][0] = S0; //initialize the simulation matrix, option value matrix and option price matrix.
		option_value[i][0] = max(X - Sx, 0.0);
		option_value[i + 1][0] = max(X - Sx, 0.0);
		option_price[i][0] = max(X - Sx, 0.0);
		option_price[i + 1][0] = max(X - Sx, 0.0);
		generate_uniform_distributed_sequence(uniform, m, seed);//generate the normal distribution sequence
		generate_standard_normal_distributed_sequence(normal, uniform);
		for (int j = 1; j <= m; j++) {
			Sx = Sx + r*Sx*delta + sigma*Sx*sqrt(delta)*normal[j - 1];// use antithetic method to simulate the paths
			simulation[i][j] = Sx; option_value[i][j] = max(X - Sx, 0.0); option_price[i][j] = max(X - Sx, 0.0);
			Sy = Sy + r*Sy*delta + sigma*Sy*sqrt(delta)*(-normal[j - 1]);
			simulation[i + 1][j] = Sy; option_value[i + 1][j] = max(X - Sy, 0.0); option_price[i + 1][j] = max(X - Sy, 0.0);
		}
		sum=sum+simulation[i][m] - simulation[i][time];
		sum = sum + simulation[i+1][m] - simulation[i+1][time];
	}
	cout << sum*exp(-r*delta*time) / n;
}

void question2b(double S0 = 65, double sigma = 0.2, double r = 0.06, int n = 10000, double X = 60, int m = 100, double T = 1) {

	/*Question 1*/

	//The initial parameters
	double** simulation;
	double** option_value;
	double** option_price;

	double delta = T / m; int k = 4; int time = 0.2 / delta;
	double discount;
	discount = exp(-r*delta);
	/*this vector is used to store the paths which have positive ECV*/
	vector<int>index_vertical;
	/*this vector is used to store the time of exercising option,starting from the end*/
	vector<int>index_horizontal(n, m);
	//first simulate 100000 paths of St
	
	simulation = new double*[n]; option_value = new double*[n]; option_price = new double*[n];
	//create a n*m matrix to store the simulated paths.
	for (int i = 0; i < n; i++) {
		simulation[i] = new double[m + 1];
		option_value[i] = new double[m + 1];// option value is the exercise value
		option_price[i] = new double[m + 1];// option price is the max of exercise value and ECV
	}
	

	for (int i = 0; i < n; i = i + 2) {
		double Sx = S0; double Sy = S0;
		vector<double> uniform, normal; int seed = i % 123465 + abs(generate_normal_number(1000, 300));
		simulation[i][0] = S0; simulation[i + 1][0] = S0; //initialize the simulation matrix, option value matrix and option price matrix.
		option_value[i][0] = max(X - Sx, 0.0);
		option_value[i + 1][0] = max(X - Sx, 0.0);
		option_price[i][0] = max(X - Sx, 0.0);
		option_price[i + 1][0] = max(X - Sx, 0.0);
		generate_uniform_distributed_sequence(uniform, m, seed);//generate the normal distribution sequence
		generate_standard_normal_distributed_sequence(normal, uniform);
		for (int j = 1; j <= m; j++) {
			Sx = Sx + r*Sx*delta + sigma*Sx*sqrt(delta)*normal[j - 1];// use antithetic method to simulate the paths
			simulation[i][j] = Sx; 
			if (Sx < 0)
				cout << Sx << endl;
			Sy = Sy + r*Sy*delta + sigma*Sy*sqrt(delta)*(-normal[j - 1]);
			simulation[i + 1][j] = Sy; 
		}
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j <= m; j++) {
			option_value[i][j] = max((simulation[i][time] - simulation[i][j]),0.0);
		}
	}
	int count_valid = 0; vector<double> x, y;
	// create the x and y for regression
	for (int i = 0; i < n; i++) {
		if (option_value[i][m]>0) {
			index_vertical.push_back(i);//change the index
			count_valid++;
			x.push_back(simulation[i][m]);//store those exercised paths in a vector
			y.push_back(discount*option_value[i][m + 1]);//store the discounted option value to the y
		}
	}
	double* ECV=new double[x.size()];
	regressMonomial(x, y, ECV, k);


	//replace the EV with ECV, if not exercise
	for (int i = 0; i < index_vertical.size(); i++) {
		if (option_value[index_vertical[i]][m]>ECV[i]) {
			option_price[index_vertical[i]][m] = option_value[index_vertical[i]][m];
			index_horizontal[index_vertical[i]] = m - 1;
		}
		else {
			option_price[index_vertical[i]][m] = ECV[i];
		}
	}
	delete[]ECV;
	for (int i = m; i > time+1; i--) {
		count_valid = 0; x.clear(), y.clear(); index_vertical.clear();
		// create the x and y for regression
		for (int j = 0; j < n; j++) {
			discount = (-r)*delta*(index_horizontal[j] + 1 - i);
			double yi = option_value[j][index_horizontal[j]] * exp(discount);// Y value

			if (option_value[j][i - 1]>0) {
				index_vertical.push_back(j);//change the index
				count_valid++;
				x.push_back(simulation[j][i - 1]);//store those exercised paths in a vector
				y.push_back(yi);//store the discounted option value to the y
			}
		}
		double* ECV = new double[x.size()];
		regressMonomial(x, y, ECV, k);
		//replace the EV with ECV, if not exercise
		for (int j = 0; j < index_vertical.size(); j++) {
			if (option_value[index_vertical[j]][i - 1]>ECV[j]) {
				option_price[index_vertical[j]][i - 1] = option_value[index_vertical[j]][i - 1];
				index_horizontal[index_vertical[j]] = i - 1;
			}
			else {
				option_price[index_vertical[j]][i - 1] = ECV[j];
			}
		}
	}
	int count = 0; double sum = 0.0;
	for (int i = 0; i < n; i++) {
		option_price[i][0] = option_value[i][index_horizontal[i]] * exp(-r*delta*(index_horizontal[i]));	
		if (option_price[i][0]>0)
			count++;
		//   cout << option_value[i][index_horizontal[i]] << " " << pow(discount, m + index_horizontal[i]) << endl;
		sum = sum+option_price[i][0];
	}
	//		cout << count<<endl;
	cout << sum / n;
}

#endif