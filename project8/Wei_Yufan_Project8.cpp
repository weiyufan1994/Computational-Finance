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
#include"Header.h"
#include"Project8.h"	
using namespace std;

int main() {
	/*Question 1*/
	double r0 = 0.05; double sigma = 0.18; double k = 0.82; double r_ = 0.05;
	// a)
	double T = 0.5; double FV = 1000;
	double price_a;
	price_a=PDB(r0, sigma, k, r_, T, FV, 10000, 180);
	cout << price_a<<endl;
	
	// b)
	double T_CB[8];
	double Coupon[8];
	for (int i = 0; i < 8; i++) {
		T_CB[i] = (i + 1)*0.5;
		Coupon[i] = 30;
	}
	Coupon[7] = 1030;
	double price_b;
	price_b=PDCB(r0, sigma, k, r_, T_CB, 8, FV, Coupon, 1000, 180);
	cout << price_b<<endl;

	//c)
	
	double option_price;
	option_price=European_call_PDB(0.25, 980, r0, sigma, k, r_, T, FV);
	double option_price2;
	option_price2 = European_call_PDCB(0.25, 980, r0, sigma, k, r_, T_CB, 8, FV, Coupon, 1000, 180);
	cout << option_price << endl;
	cout << option_price2 << endl;

	/*Question 2*/
	// a
	k = 0.92; r_ = 0.055;
	double option_price_2a;
	option_price_2a = European_call_PBD_CIR_MC(0.5, 980, r0, sigma, k, r_, T, FV);
	cout << option_price_2a << endl;
	
	
	// b
	double option_price_2b;
	option_price_2b = European_call_PDB_CIR(0.5, 980, r0, sigma, k, r_, T, FV);
	cout << option_price_2b << endl;
	/*Question 3*/
	double option_price_3;
	option_price_3=G2pp_option(985, 0.5, 0, 0, 0.03, 0.03, 0.7, 0.1, 0.3, 0.03, 0.08, 1);
	cout << option_price_3 << endl;
	system("pause");
	return 0;
}