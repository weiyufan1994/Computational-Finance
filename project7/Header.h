#ifndef HEADER_H_
#define HEADER_H_

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

 void generate_uniform_distributed_sequence(vector<double> &RV, int N, int seed) {

		int m = pow(2, 31) - 1;
		int a = pow(7, 5);
		int b = 0;

		int x_0 = seed; // initialize a seed for the random numbers
		unsigned long int x = (a*x_0 + b) % m;//random number generated
		RV.push_back(x);//Add generated number to the end of the vector
						// generate more random numbers
		for (int i = 1; i < N; i++) {
			x = (a*x + b) % m;
			RV.push_back(x);
		}
		// scale the generated numbers to the uniform distribution variables
		for (int i = 0; i < N; i++) {
			RV[i] = RV[i] / m;
		}
	}
	// scale the generated numbers to the uniform distribution variables


	int generate_Bernoulli_number(double x, double p) {
		if (x <= p)
			return 1;
		else
			return 0;
	}

	//calculate the mean and standard deviation
	void calculate_sum_mean_stdev(vector<double> &RV1,double &sum,double &mean,double &stdev) {
		double mean1 = accumulate(begin(RV1), end(RV1), 0.0) / RV1.size();
		double sum1 = 0;
		for (int i = 0; i < RV1.size(); i++) {
			sum1 += (RV1[i] - mean1)*(RV1[i] - mean1);
		}
		double stdev1 = sqrt(sum1 / RV1.size());
		sum = sum1; mean = mean1; stdev = stdev1;
	}

	void generate_standard_normal_distributed_sequence(vector<double> &normal_distributed_sequence, vector<double> &uniform_sequence) {	
		for (int i = 0; i < uniform_sequence.size()/2; i++) {
			double Z_1 = sqrt(-2 * log(uniform_sequence[i]))*cos(2 * (atan(1) * 4)*uniform_sequence[i + uniform_sequence.size() / 2]);
			double Z_2 = sqrt(-2 * log(uniform_sequence[i]))*sin(2 * (atan(1) * 4)*uniform_sequence[i + uniform_sequence.size() / 2]);
			normal_distributed_sequence.push_back(Z_1);
			normal_distributed_sequence.push_back(Z_2);
		}
	}
	double generate_normal_number(double mean=0, double standardDeviation=1) {
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine generator(seed);
		std::normal_distribution<double> distribution(mean, standardDeviation);
		double x = distribution(generator);
		return x;
	}

	double european_call_option_price_mts(double r, double sigma, double S0, double T, double X, int n = 10000000, int seed = 86534) {
		double total = 0;
		vector<double> UR;//uniform
		vector<double> NR;//normal
		generate_uniform_distributed_sequence(UR, n, seed);
		generate_standard_normal_distributed_sequence(NR, UR);//generate normal sequence
		for (int i = 0; i < n; i++) {
			double e = sigma*NR[i] * sqrt(T) + (r - 0.5*sigma*sigma)*T;
			double St = S0*exp(e);
			double cp = ((St - X) > 0) ? (St - X) : 0;
			double discount_factor = -r*T;
			cp = cp*exp(discount_factor);
			total += cp;
		}
		return total / n;
	}
	//derive call option price through antithetic method
	double european_call_option_price_vr(double r, double sigma, double S0, double T, double X, int n = 1000000, int seed = 618754) {
		double total = 0;
		vector<double> UR;//uniform
		vector<double> NR;//normal

		generate_uniform_distributed_sequence(UR, n, seed);
		generate_standard_normal_distributed_sequence(NR, UR);

		for (int i = 0; i < n; i++) {
			double e1 = sigma*NR[i] * sqrt(T) + (r - 0.5*sigma*sigma)*T;
			double e = sigma*(-NR[i])*sqrt(T) + (r - 0.5*sigma*sigma)*T;
			double St = S0*exp(e);
			double St1 = S0*exp(e1);
			double cp = ((St - X) > 0) ? (St - X) : 0;
			double cm = ((St1 - X) > 0) ? (St1 - X) : 0;
			double discount_factor = -r*T;
			cp = cp*exp(discount_factor); cm = cm*exp(discount_factor);
			total = total + (cp + cm) / 2.0;
		}
		return total / n;
	}
	double european_call_option_price_BSM(double r, double sigma, double S0, double T, double X) {
		double d1 = 1.0 / (sigma*sqrt(T))*(log(S0 / X) + (r + sigma*sigma / 2.0)*T);
		double d2 = d1 - sigma*sqrt(T);
		double Nd1, Nd2;
		if (d1 >= 0) {
			Nd1 = 1.0 - 0.5*pow((1.0 + 0.0498673480*d1 + 0.0211410061*d1*d1 + 0.0032776263*d1*d1*d1 + 0.0000380036*d1*d1*d1*d1 + 0.0000488906*d1*d1*d1*d1*d1 + 0.0000053830*d1*d1*d1*d1*d1*d1), -16);
		}
		else {
			d1 = -d1;
			Nd1 = 1.0 - 1 + 0.5*pow((1.0 + 0.0498673480*d1 + 0.0211410061*d1*d1 + 0.0032776263*d1*d1*d1 + 0.0000380036*d1*d1*d1*d1 + 0.0000488906*d1*d1*d1*d1*d1 + 0.0000053830*d1*d1*d1*d1*d1*d1), -16);
		}
		if (d2 >= 0) {
			Nd2 = 1.0 - 0.5*pow((1.0 + 0.0498673480*d2 + 0.0211410061*d2*d2 + 0.0032776263*d2*d2*d2 + 0.0000380036*d2*d2*d2*d2 + 0.0000488906*d2*d2*d2*d2*d2 + 0.0000053830*d2*d2*d2*d2*d2*d2), -16);
		}
		else {
			d2 = -d2;
			Nd2 = 0.5*pow((1.0 + 0.0498673480*d2 + 0.0211410061*d2*d2 + 0.0032776263*d2*d2*d2 + 0.0000380036*d2*d2*d2*d2 + 0.0000488906*d2*d2*d2*d2*d2 + 0.0000053830*d2*d2*d2*d2*d2*d2), -16);
		}
		double discount_factor = -r*T;
		return Nd1*S0 - Nd2*X*exp(discount_factor);
	}

	double halton_sequence(int n, int base) {
		double f = 1, r = 0;
		while (n > 0) {
			f = f / base;
			r = r + f*(n%base);
			n = n / base;
		}
		return r;
	}

	void set_seed(int &seed) {
		time_t nowtime;
		nowtime = time(NULL); //take time as seed;
		seed = nowtime;
	}

	double call_price_Heston(string method, double T = 0.5, double X = 50.0, double tho = -0.6, double r = 0.03,
		double S0 = 48.0, double V0 = 0.05, double sigma = 0.42, double a = 5.8, double b = 0.0625, double n = 4000) {
		double delta = T / n; double sum = 0; double discount_factor = -r*T;
		// take the method as an input, decide which method should be used
		if (method == "Full Truncation") {
			for (int i = 0; i < n; i++) {
				double St = S0; double Vt = V0;
				vector<double> UR;//uniform
				vector<double> NR;//normal
				int seed = abs(generate_normal_number(1000, 300)) + i; // scramble seed
				generate_uniform_distributed_sequence(UR, n, seed); //generate normal sequence
				generate_standard_normal_distributed_sequence(NR, UR);
				vector<double> UR1;//uniform
				vector<double> NR1;//normal
				seed = abs(generate_normal_number(1000, 300)) + 11 + i;
				generate_uniform_distributed_sequence(UR1, n, seed);
				generate_standard_normal_distributed_sequence(NR1, UR1);
				for (int k = 0; k < n; k++) {
					St = St + r*St*delta + sqrt((Vt>0 ? Vt : 0))*St*sqrt(delta)*NR[k];
					Vt = Vt + a*(b - (Vt>0 ? Vt : 0))*delta + sigma*sqrt(delta)*sqrt((Vt>0 ? Vt : 0))*(tho*NR[k] + sqrt(1.0 - tho*tho)*NR1[k]);
				}
				double c = ((St - X) > 0) ? (St - X) : 0;
				sum += c*exp(discount_factor);
			}
		}
		else if (method == "Partial Truncation") {
			for (int i = 0; i < n; i++) {
				double St = S0; double Vt = V0;
				vector<double> UR;//uniform
				vector<double> NR;//normal
				int seed = abs(generate_normal_number(1000, 300)) + (2 * i) % 17;
				generate_uniform_distributed_sequence(UR, n, seed);
				generate_standard_normal_distributed_sequence(NR, UR);
				vector<double> UR1;//uniform
				vector<double> NR1;//normal
				seed = abs(generate_normal_number(1000, 300)) + i % 2;
				generate_uniform_distributed_sequence(UR1, n, seed);
				generate_standard_normal_distributed_sequence(NR1, UR1);
				for (int k = 0; k < n; k++) {
					St = St + r*St*delta + sqrt((Vt>0 ? Vt : 0))*St*sqrt(delta)*NR[k];
					Vt = Vt + a*(b - Vt)*delta + sigma*sqrt((Vt>0 ? Vt : 0))*sqrt(delta)*(tho*NR[k] + sqrt(1.0 - tho*tho)*NR1[k]);
				}
				double c = ((St - X) > 0) ? (St - X) : 0;

				sum += c*exp(discount_factor);
			}
		}
		else if (method == "Reflection") {
			for (int i = 0; i < n; i++) {
				double St = S0; double Vt = V0;
				vector<double> UR;//uniform
				vector<double> NR;//normal
				int seed = abs(generate_normal_number(1000, 300)) + (i * 5) % 17;
				generate_uniform_distributed_sequence(UR, n, seed);
				generate_standard_normal_distributed_sequence(NR, UR);
				vector<double> UR1;//uniform
				vector<double> NR1;//normal
				seed = abs(generate_normal_number(1000, 300)) + (i * 6) % 13;
				generate_uniform_distributed_sequence(UR1, n, seed);
				generate_standard_normal_distributed_sequence(NR1, UR1);
				for (int k = 0; k < n; k++) {
					St = St + r*St*delta + sqrt(abs(Vt))*St*sqrt(delta)*NR[k];
					Vt = abs(Vt) + a*(b - abs(Vt))*delta + sigma*sqrt(delta)*sqrt(abs(Vt))*(tho*NR[k] + sqrt(1.0 - tho*tho)*NR1[k]);
				}
				double c = ((St - X) > 0) ? (St - X) : 0;

				sum += c*exp(discount_factor);
			}
		}
		else
			cout << "re-enter the name." << endl;
		return sum / n;
	}
	double Halton_integral(int base1, int base2, int n = 10000) {
		double delta = 1 / n;
		double x = 0.0; double y = 0.0;
		double* QMCx = new double[n];
		double* QMCy = new double[n];
		for (int i = 0; i < n; i++) {
			QMCx[i] = halton_sequence(i + 1, base1);
			QMCy[i] = halton_sequence(i + 1, base2);
		}
		double sum = 0;
		for (int i = 0; i < n; i++) {
			sum += exp((-QMCx[i])*QMCy[i])*(sin(6.0 * (atan(1) * 4)*QMCx[i]) + cbrt(cos(2.0 * (atan(1) * 4)*QMCy[i])));
		}

		return sum / n;
	}

	double vector_dot_prod(vector<double> &A, double* B) {
		double sum = 0.0;
		for (int i = 0; i < A.size(); i++) {
			sum += A[i] * B[i];
		}
		return sum;
	}


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
		array = new double* [n];
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
#endif