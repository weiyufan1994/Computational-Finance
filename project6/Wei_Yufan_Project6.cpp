#include"Header.h"
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

double call_option_price(double sigma, double r = 0.03, double S0 = 98, double X = 100) {
	double Sx = S0; double delta = 1 / 1000; double Smax=S0;
	vector<double> uniform; vector<double> normal;
	generate_uniform_distributed_sequence(uniform, 1000, 1 + abs(generate_normal_number(100, 30)));
	generate_standard_normal_distributed_sequence(normal, uniform);
	for (int i = 0; i < 1000; i++) {
		Sx = Sx + Sx*r*delta + Sx*sqrt(delta)*sigma*normal[i];
		Smax = (Sx > Smax) ? Sx : Smax;
	}
	double price;
	price = Smax - X;
	return price*exp(-r);
}
double put_option_price(double sigma, double r = 0.03, double S0 = 98, double X = 100) {
	double Sx = S0; double delta = 1 / 1000; double Smin = S0;
	vector<double> uniform; vector<double> normal;
	generate_uniform_distributed_sequence(uniform, 1000, 1 + abs(generate_normal_number(100, 30)));
	generate_standard_normal_distributed_sequence(normal, uniform);
	for (int i = 0; i < 1000; i++) {
		Sx = Sx + Sx*r*delta + Sx*sqrt(delta)*sigma*normal[i];
		Smin = (Sx < Smin) ? Sx : Smin;
	}
	double price;
	price = X-Smin;
	return price*exp(-r);
}

int main() {
	/*Question 1*/
	for (int j = 12; j <= 48; j = j + 4) {
		double sum = 0;
		for (int i = 0; i < 10000; i++) {
			sum += call_option_price(j);
		}
		cout << j << endl;
		cout << sum / 10000 << endl;
	}

	for (int j = 12; j <= 48; j = j + 4) {
		double sum = 0;
		for (int i = 0; i < 10000; i++) {
			sum += put_option_price(j);
		}
		cout << j << endl;
		cout << sum / 10000 << endl;
	}

	/* Question 2*/


}