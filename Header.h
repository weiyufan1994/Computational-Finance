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
#endif