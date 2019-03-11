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
	vector<Laguerre> l;
	for (int i = 1; i <= k; i++) {
		Laguerre f(i);
		l.push_back(f);
	}
	for (int i = 0; i < k; i++) {
		sum += a[i] * l[i].poly(s);
	}
	return sum;
}
double ysM(double a[], int k, double s) {
	 double sum = 0;
	vector<Monomials> l;
	for (int i=1; i <= k; i++) {
		Monomials f(i);
		l.push_back(f);
	}
	for (int i = 0; i < k; i++) {
		sum += a[i] * l[i].poly(s);
	}
	return sum;
}
double ysH(double a[], int k, double s) {
	 double sum = 0;
	vector<Hermite> l;
	for (int i = 1; i <= k; i++) {
		Hermite f(i);
		l.push_back(f);
	}
	for (int i = 0; i < k; i++) {
		sum += a[i] * l[i].poly(s);
	}
	return sum;
}
void question1L(double S0 = 36, double sigma = 0.2, double r = 0.06, int n = 100000, double X = 40, int m = 100, int k = 3, double T = 0.5) {

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

	for (int i = 0; i < n; i = i + 2) {
		double Sx = S0; double Sy = S0;
		vector<double> uniform, normal; int seed = i+1;
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
	}


	/*The last period*/
	//When t=T, check if the option will be exercised
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


	//Laguerre polynomials
	vector<Laguerre> l;//create a vector to store k Laguerre polynomials,
	for (int i = 1; i <= k; i++) {
		Laguerre li(i);// create a object Li, and there will be L1, L2 ..... Lk
		l.push_back(li);
	}


	double* b = new double[k];

	vector<vector<double>> F;//create basis function
	for (int j = 0; j < k; j++) {
		vector<double> Fj;
		for (int i = 0; i < count_valid; i++) {
			Fj.push_back(l[j].poly(x[i]));
		}
		F.push_back(Fj);
	}

	for (int i = 0; i < k; i++) {
		b[i] = vector_dot_prod(y, F[i]);
	}

	double** A; double** A_inv; A = new double*[k]; A_inv = new double*[k];


	for (int i = 0; i < k; i++) {
		A[i] = new double[k];
		A_inv[i] = new double[k];
	}

	for (int i = 0; i < k; i++) {
		for (int j = 0; j < k; j++) {
			A[i][j] = vector_dot_prod(F[i], F[j]);
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
	vector<double> ECV;
	for (int i = 0; i < count_valid; i++) {
		ECV.push_back(ysM(a, k, x[i]));
	}
	delete[] a;
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

	for (int i = m; i > 1; i--) {
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

		l.clear();
		for (int i = 1; i <= k; i++) {
			Laguerre li(i);
			l.push_back(li);
		}

		double* b = new double[k] {0};

		F.clear();//create basis function
		for (int j = 0; j < k; j++) {
			vector<double> Fj;
			for (int i = 0; i < count_valid; i++) {
				Fj.push_back(l[j].poly(x[i]));
			}
			F.push_back(Fj);
		}
		for (int i = 0; i < k; i++) {
			b[i] = vector_dot_prod(y, F[i]);
		}

		double** A; double** A_inv; A = new double*[k]; A_inv = new double*[k];
		for (int i = 0; i < k; i++) {
			A[i] = new double[k];
			A_inv[i] = new double[k];
		}

		for (int i = 0; i < k; i++) {
			for (int j = 0; j < k; j++) {
				A[i][j] = vector_dot_prod(F[j], F[i]);
			}
		}

		//inverse_matrix(A, k, A_inv);
		double* a = new double[k] {0};
		//	matrix_multiply(A_inv, k, b, a);//calculate the a;

		a = coef(A, b, k);
		for (int i = 0; i < k; i++) {
			delete[] A[i];
			delete[] A_inv[i];
		}
		delete[] A; delete[] A_inv; delete[] b;
		//calculate the expected continous value
		ECV.clear();
		for (int i = 0; i < count_valid; i++) {
			ECV.push_back(ysM(a, k, x[i]));
		}		delete[]a;
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
	int count = 0; double sum = 0;
	for (int i = 0; i < n; i++) {
		option_price[i][0] = option_value[i][index_horizontal[i]] * exp(-r*delta*(index_horizontal[i]));
		sum += option_price[i][0];
		//   cout << option_value[i][index_horizontal[i]] << " " << pow(discount, m + index_horizontal[i]) << endl;
		//cout << option_value[i][m]<<endl;
	}
	//		cout << count<<endl;
	cout << sum / n;
}

void question1H(double S0 = 36, double sigma = 0.2, double r = 0.06, int n = 100000, double X = 40, int m = 100, int k = 3, double T = 0.5) {

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
	}


	/*The last period*/
	//When t=T, check if the option will be exercised
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


	//Laguerre polynomials
	vector<Hermite> l;//create a vector to store k Laguerre polynomials,
	for (int i = 1; i <= k; i++) {
		Hermite li(i);// create a object Li, and there will be L1, L2 ..... Lk
		l.push_back(li);
	}


	double* b = new double[k];

	vector<vector<double>> F;//create basis function
	for (int j = 0; j < k; j++) {
		vector<double> Fj;
		for (int i = 0; i < count_valid; i++) {
			Fj.push_back(l[j].poly(x[i]));
		}
		F.push_back(Fj);
	}

	for (int i = 0; i < k; i++) {
		b[i] = vector_dot_prod(y, F[i]);
	}

	double** A; double** A_inv; A = new double*[k]; A_inv = new double*[k];


	for (int i = 0; i < k; i++) {
		A[i] = new double[k];
		A_inv[i] = new double[k];
	}

	for (int i = 0; i < k; i++) {
		for (int j = 0; j < k; j++) {
			A[i][j] = vector_dot_prod(F[i], F[j]);
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
	vector<double> ECV;
	for (int i = 0; i < count_valid; i++) {
		ECV.push_back(ysM(a, k, x[i]));
	}
	delete[] a;
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

	for (int i = m; i > 1; i--) {
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

		l.clear();
		for (int i = 1; i <= k; i++) {
			Hermite li(i);
			l.push_back(li);
		}

		double* b = new double[k] {0};

		F.clear();//create basis function
		for (int j = 0; j < k; j++) {
			vector<double> Fj;
			for (int i = 0; i < count_valid; i++) {
				Fj.push_back(l[j].poly(x[i]));
			}
			F.push_back(Fj);
		}
		for (int i = 0; i < k; i++) {
			b[i] = vector_dot_prod(y, F[i]);
		}

		double** A; double** A_inv; A = new double*[k]; A_inv = new double*[k];
		for (int i = 0; i < k; i++) {
			A[i] = new double[k];
			A_inv[i] = new double[k];
		}

		for (int i = 0; i < k; i++) {
			for (int j = 0; j < k; j++) {
				A[i][j] = vector_dot_prod(F[j], F[i]);
			}
		}

		//inverse_matrix(A, k, A_inv);
		double* a = new double[k] {0};
		//	matrix_multiply(A_inv, k, b, a);//calculate the a;

		a = coef(A, b, k);
		for (int i = 0; i < k; i++) {
			delete[] A[i];
			delete[] A_inv[i];
		}
		delete[] A; delete[] A_inv; delete[] b;
		//calculate the expected continous value
		ECV.clear();
		for (int i = 0; i < count_valid; i++) {
			ECV.push_back(ysM(a, k, x[i]));
		}		delete[]a;
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
	int count = 0;
	for (int i = 0; i < n; i++) {
		option_price[i][0] = option_value[i][index_horizontal[i]] * exp(-r*delta*(index_horizontal[i]));

		if (option_price[i][0]>0)
			count++;
		//   cout << option_value[i][index_horizontal[i]] << " " << pow(discount, m + index_horizontal[i]) << endl;
		//cout << option_value[i][m]<<endl;
	}
	//		cout << count<<endl;

	double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += option_price[i][0];
	}
	cout << sum / n;
}

void question1M(double S0 = 36, double sigma = 0.2, double r = 0.06, int n = 100000, double X = 40, int m = 100,int k = 3,double T = 0.5) {

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
	}


	/*The last period*/
	//When t=T, check if the option will be exercised
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


	//Laguerre polynomials
	vector<Monomials> l;//create a vector to store k Laguerre polynomials,
	for (int i = 1; i <= k; i++) {
		Monomials li(i);// create a object Li, and there will be L1, L2 ..... Lk
		l.push_back(li);
	}


	double* b = new double[k];

	vector<vector<double>> F;//create basis function
	for (int j = 0; j < k; j++) {
		vector<double> Fj;
		for (int i = 0; i < count_valid; i++) {
			Fj.push_back(l[j].poly(x[i]));
		}
		F.push_back(Fj);
	}

	for (int i = 0; i < k; i++) {
		b[i] = vector_dot_prod(y, F[i]);
	}

	double** A; double** A_inv; A = new double*[k]; A_inv = new double*[k];


	for (int i = 0; i < k; i++) {
		A[i] = new double[k];
		A_inv[i] = new double[k];
	}

	for (int i = 0; i < k; i++) {
		for (int j = 0; j < k; j++) {
			A[i][j] = vector_dot_prod(F[i], F[j]);
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
	vector<double> ECV;
	for (int i = 0; i < count_valid; i++) {
		ECV.push_back(ysM(a, k, x[i]));
	}
	delete[] a;
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

	for (int i = m; i > 1; i--) {
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

		l.clear();
		for (int i = 1; i <= k; i++) {
			Monomials li(i);
			l.push_back(li);
		}

		double* b = new double[k] {0};

		F.clear();//create basis function
		for (int j = 0; j < k; j++) {
			vector<double> Fj;
			for (int i = 0; i < count_valid; i++) {
				Fj.push_back(l[j].poly(x[i]));
			}
			F.push_back(Fj);
		}
		for (int i = 0; i < k; i++) {
			b[i] = vector_dot_prod(y, F[i]);
		}

		double** A; double** A_inv; A = new double*[k]; A_inv = new double*[k];
		for (int i = 0; i < k; i++) {
			A[i] = new double[k];
			A_inv[i] = new double[k];
		}

		for (int i = 0; i < k; i++) {
			for (int j = 0; j < k; j++) {
				A[i][j] = vector_dot_prod(F[j], F[i]);
			}
		}

		//inverse_matrix(A, k, A_inv);
		double* a = new double[k] {0};
		//	matrix_multiply(A_inv, k, b, a);//calculate the a;

		a = coef(A, b, k);
		for (int i = 0; i < k; i++) {
			delete[] A[i];
			delete[] A_inv[i];
		}
		delete[] A; delete[] A_inv; delete[] b;
		//calculate the expected continous value
		ECV.clear();
		for (int i = 0; i < count_valid; i++) {
			ECV.push_back(ysM(a, k, x[i]));
		}		delete[]a;
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
	int count = 0;
	for (int i = 0; i < n; i++) {
		option_price[i][0] = option_value[i][index_horizontal[i]] * exp(-r*delta*(index_horizontal[i]));

		if (option_price[i][0]>0)
			count++;
		//   cout << option_value[i][index_horizontal[i]] << " " << pow(discount, m + index_horizontal[i]) << endl;
		//cout << option_value[i][m]<<endl;
	}
	//		cout << count<<endl;

	double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += option_price[i][0];
	}
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
void regressMonomial(vector<double> x,vector<double> y, vector<double> &ECV,int k) {
	//Laguerre polynomials
	vector<Monomials> l;//create a vector to store k Laguerre polynomials,
	for (int i = 1; i <= k; i++) {
		Monomials li(i);// create a object Li, and there will be L1, L2 ..... Lk
		l.push_back(li);
	}


	double* b = new double[k];

	vector<vector<double>> F;//create basis function
	for (int j = 0; j < k; j++) {
		vector<double> Fj;
		for (int i = 0; i < x.size(); i++) {
			Fj.push_back(l[j].poly(x[i]));
		}
		F.push_back(Fj);
	}

	for (int i = 0; i < k; i++) {
		b[i] = vector_dot_prod(y, F[i]);
	}

	double** A; double** A_inv; A = new double*[k]; A_inv = new double*[k];


	for (int i = 0; i < k; i++) {
		A[i] = new double[k];
		A_inv[i] = new double[k];
	}

	for (int i = 0; i < k; i++) {
		for (int j = 0; j < k; j++) {
			A[i][j] = vector_dot_prod(F[i], F[j]);
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
		ECV.push_back(ysM(a, k, x[i]));
	}
	delete[] a;
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
	vector<double> ECV;
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
		ECV.clear();
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