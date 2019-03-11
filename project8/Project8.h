#pragma once
#ifndef PROJECT8_H_
#define PROJECT8_H_
#include<iostream>
#include<ctime>
#include<time.h>
#include<chrono>
#include"Header.h"
#include"Generater.h"
using namespace std;

// pure discount bond
double PDB(double r0, double sigma, double k, double r_, double T, double FV, double N, int step) {
	 double delta; double Ri=0;
	 delta = T / step;
	 Generater* g = new Generater(2);
	for (int j = 0; j < N; j++) {
		double rt = r0;
		double sum_r = rt;
		double* normals = g->generateDefaultNormal(step + 1);
		for (int i = 0; i < step; i++) {
			rt = rt + k*(r_ - rt)*delta + sigma*sqrt(delta)*normals[i];
			sum_r += rt;
		}
		delete[] normals;
		double R;
		R = sum_r*delta;
		Ri += exp(-R);
	}
	double price = Ri / N*FV;
	return price;
}
//pure discount coupon bond
double PDCB(double r0, double sigma, double k, double r_, double T[],int period, double FV,double coupon[], double N, int step) {
	double price = 0;
	for (int i = 0; i < period; i++) {
		double prices;
		prices= PDB(r0, sigma, k, r_, T[i], coupon[i], N, step*(i+1));
		price += prices;
	}
	return price;
}

double PDB_explicit(double r0, double sigma, double k, double r_, double T, double FV) {
	double B;
	B = 1 / k*(1 - exp(-k*T));
	double A;
	A = exp((r_ - sigma*sigma / k / k / 2)*(B - T) - sigma*sigma / 4 / k*B*B);
	double price;
	price = FV*A*exp(-B*r0);
	return price;
}

double European_call_PDB(double maturity, double strike, double r0, double sigma, double k, double r_, double T, double FV) {
	int N = 10000;//simulate 10000 paths
	 int step = maturity * 365; double delta = maturity / step;
	double path[10000]; double path_sum[10000]{ 0 };
	Generater* g = new Generater(2);
	for (int i = 0; i < 10000; i++) {
		double rt=r0; 
		double* normals = g->generateDefaultNormal(step + 1);
		for (int j = 0; j < step; j++) {
			rt = rt + k*(r_ - rt)*delta + sigma*sqrt(delta)*normals[j];
			path_sum[i] += rt;
		}
		path[i] = rt;
		path_sum[i] = path_sum[i] * delta;
	}

	double prices[10000]; double sum = 0;
	for (int i = 0; i < 10000; i++) {
		prices[i]=max(PDB_explicit(path[i], sigma, k, r_, T, FV)-980,0.0);
		prices[i] = exp(-path_sum[i])*prices[i];
		sum += prices[i];
	}
	
	return sum / 10000;
	
}

double European_call_PDCB(double maturity, double strike, double r0, double sigma, double k, double r_, double T[], int period, double FV, double coupon[], double N, int step) {
	 int steps = maturity * 365; double delta = maturity / steps;
	double path[1000]; double path_sum[1000]{ 0 };
	Generater* g = new Generater(2);
	for (int i = 0; i < 1000; i++) {
		double rt=r0; 
		double* normals = g->generateDefaultNormal(steps + 1);
		for (int j = 0; j < steps; j++) {
			rt = rt + k*(r_ - rt)*delta + sigma*sqrt(delta)*normals[j];
			path_sum[i] += rt;
		}
		path[i] = rt;
		path_sum[i] = path_sum[i] * delta;
	}

	double prices[1000]; double sum = 0;
	for (int i = 0; i < 1000; i++) {
		prices[i]=max(PDCB(path[i], sigma, k, r_, T, period,FV,coupon,N,step)-980,0.0);
		prices[i] = exp(-path_sum[i])*prices[i];
		sum += prices[i];
	}
	
	return sum / 1000;
	
}

double PDB_explicit_CIR(double r0, double sigma, double k, double r_, double T, double FV) {
	double h1 = sqrt(k*k + 2 * sigma*sigma);
	double h2 = (k + h1) / 2;
	double h3 = (2 * k*r_) / sigma / sigma;
	double A, B;
	A = h1*exp(h2*T) / (h1 + h2*(exp(h1*T) - 1));
	A = pow(A, h3);
	B = (exp(h1*T) - 1) / (h1 + h2*(exp(h1*T) - 1));
	double price;
	price = A*exp(-B*r0)*FV;
	return price;
}

double European_call_PDB_CIR(double maturity, double strike, double r0, double sigma, double k, double r_, double T, double FV) {
	int N = 10000;//simulate 10000 paths
	int step = maturity * 365; double delta = maturity / step;
	double path[10000]; double path_sum[10000]{ 0 };
	Generater* g = new Generater(2);
	for (int i = 0; i < 10000; i++) {
		double rt = r0;
		double* normals = g->generateDefaultNormal(step + 1);
		for (int j = 0; j < step; j++) {
			rt = rt + k*(r_ - rt)*delta + sqrt(rt)*sigma*sqrt(delta)*normals[j];
			path_sum[i] += rt;
		}
		path[i] = rt;
		path_sum[i] = path_sum[i] * delta;
	}

	double prices[10000]; double sum = 0;
	for (int i = 0; i < 10000; i++) {
		prices[i] = max(PDB_explicit_CIR(path[i], sigma, k, r_, T, FV) - 980, 0.0);
		prices[i] = exp(-path_sum[i])*prices[i];
		sum += prices[i];
	}

	return sum / 10000;

}

void getAB(double r0, double sigma, double k, double r_, double T,double &A,double &B) {
	double h1 = sqrt(k*k + 2 * sigma*sigma);
	double h2 = (k + h1) / 2;
	double h3 = (2 * k*r_) / sigma / sigma;
	A = h1*exp(h2*T) / (h1 + h2*(exp(h1*T) - 1));
	A = pow(A, h3);
	B = (exp(h1*T) - 1) / (h1 + h2*(exp(h1*T) - 1));
}
double PDB_CIR(double r0, double sigma, double k, double r_, double T, double FV, double N, int step) {
	double delta; double Ri = 0;
	delta = T / step;
	Generater* g = new Generater(2);
	for (int j = 0; j < N; j++) {
		double rt = r0;
		double sum_r = rt;
		double* normals = g->generateDefaultNormal(step + 1);
		for (int i = 0; i < step; i++) {
			rt = rt + k*(r_ - rt)*delta + sigma*sqrt(abs(rt))*sqrt(delta)*normals[i];
			sum_r += rt;
		}
		delete[] normals;
		double R;
		R = sum_r*delta;
		Ri += exp(-R);
	}
	double price = Ri / N*FV;
	return price;
}
double European_call_PBD_CIR_MC(double maturity, double strike, double r0, double sigma, double k, double r_, double T, double FV) {
	int steps = maturity * 365; double delta = maturity / steps;
	double path[1000]; double path_sum[1000]{ 0 };
	Generater* g = new Generater(2);
	for (int i = 0; i < 1000; i++) {
		double rt = r0;
		double* normals = g->generateDefaultNormal(steps + 1);
		for (int j = 0; j < steps; j++) {
			rt = rt + k*(r_ - rt)*delta + sqrt(abs(rt))*sigma*sqrt(delta)*normals[j];
			path_sum[i] += rt;
		}
		delete[] normals;
		path[i] = rt;
		path_sum[i] = path_sum[i] * delta;
	}

	double prices[1000]; double sum = 0;
	for (int i = 0; i < 1000; i++) {
		prices[i] = max(PDB_CIR(path[i],sigma,k,r_,T,FV,1000,steps) - 980, 0.0);
		prices[i] = exp(-path_sum[i])*prices[i];
		sum += prices[i];
	}
	return sum / 1000;
}

double G2pp_PBD(double x, double y, double phi, double r0, double tho, double a, double b, double sigma, double n, double T,double FV) {
	double xt, yt; double delta = 1.0 / 365; int step = T / delta;
	  Generater* g = new Generater(2); Generater* g2 = new Generater(5); double Ri = 0;
	for (int j = 0; j < 1000; j++) {
		double sum = 0; double rt = r0; xt = x; yt = y;
		double* normals = g->generateDefaultNormal(step + 1);
		double* normals2 = g2->generateDefaultNormal(step + 1);
		for (int i = 0; i < step; i++) {
			xt = xt - a*xt*delta + sigma*sqrt(delta)*normals[i];
			yt = yt - b*yt*delta + n*sqrt(delta)*(tho*normals[i] + sqrt(1 - tho*tho)*normals2[i]);
			rt = xt + yt + phi;
			sum += rt;
		}
		delete[] normals;
		delete[] normals2;
		double R;
		R = sum*delta;
		Ri += exp(-R);
	}
	double price = Ri /1000 *FV;
	return price;
}

double G2pp_option(double strike,double maturity,double x, double y, double phi, double r0, double tho, double a, double b, double sigma, double n, double T) {
	double path[1000][2]{ 0 }; Generater* g = new Generater(2); Generater* g2 = new Generater(5);
	double path2[1000][2]{ 0 };
	double delta = 1.0 / 365; int step = maturity / delta;
	for (int i = 0; i < 1000; i++) {
		double* normals = g->generateDefaultNormal(step + 1);
		double* normals2 = g2->generateDefaultNormal(step + 1);
		double xt = x; double yt = y; double rt = r0;
		for (int j = 0; j < step; j++) {
			xt = xt - a*xt*delta + sqrt(delta)*sigma*normals[j];
			yt = yt - b*yt*delta + sqrt(delta)*n*(tho*normals[j] + sqrt(1 - tho*tho)*normals2[j]);
			rt = xt + yt + phi;
			path[i][0] = rt;
			path[i][1] += rt;
			path2[i][0] = xt;
			path2[i][1] = yt;
		}
		delete[] normals;
		delete[] normals2;
	}	
	double prices[1000]; double Ri = 0;
	for(int i = 0; i < 1000; i++) {
		prices[i] = max(strike-G2pp_PBD(path2[i][0],path2[i][1], phi, path[i][0], tho, a, b, sigma, n, T-maturity, 1000),0.0);
		prices[i] = prices[i] * exp(-path[i][1] * delta);
		Ri += prices[i];
	}
	return Ri / 1000;
}

#endif