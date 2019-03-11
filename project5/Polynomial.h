#pragma once
#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

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
int fac(int x){
		int f;
		if (x == 0 || x == 1)
			f = 1;
		else
			f = fac(x - 1)*x;
		return f;
	}


class Laguerre {
private:
	int k;
public:
	Laguerre(int n) {
		k = n;
	};
	//polynomial function
	double poly(double x) {
		if (k == 1) {
			return exp(-x / 2.0);
		}
		else if (k == 2) {
			return exp(-x / 2.0)*(1 - x);
		}
		else if (k == 3) {
			return exp(-x / 2)*(1 - 2.0*x + x*x / 2.0);
		}
		else if (k == 4) {
			return exp(-x / 2.0)*(1 - 3 * x + 3 * x*x / 2.0 - x*x*x / 6.0);
		}
	}
};
class Hermite {
private:
	int k;
public:
	Hermite(int n) {
		k = n;
	};
	//polynomial function
	double poly(double x) {
		if (k == 1) {
			return 1;
		}
		else if (k == 2) {
			return 2.0*x;
		}
		else if (k == 3) {
			return 4*x*x-2.0;
		}
		else if (k == 4) {
			return 8.0*x*x*x-12.0*x;
		}
	}
};

class Monomials {
private:
	int k;
public:
	Monomials(int n) {
		k = n;
	};
	//polynomial function
	double poly(double x) {
		if (k == 1) {
			return 1;
		}
		else if (k == 2) {
			return x;
		}
		else if (k == 3) {
			return x*x;
		}
		else if (k == 4) {
			return x*x*x;
		}
	}
};
#endif