#include"Header.h"
#include"inverse.h"
#include"MatMult.h"
#include"Header5.h"
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

double** simulation; double** option_value;
int n = 100000; int m = 200;
vector<vector<int>> x;

int main() {
	simulation = new double*[n];
	option_value = new double*[n];
	vector<vector<int>> x;
	//create a n*m matrix to store the simulated paths.
	for (int i = 0; i < n; i++) {
		simulation[i] = new double[m + 1];
		option_value[i] = new double[m + 1];// option value is the exercise value
	}

		/*Question 1*/
	double parameter1a[3][3]={36,40,44,2,3,4,0.5,1,2};
		//The initial parameters
	
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int h = 0; h < 3; h++) {
				cout << "The current price is " << parameter1a[0][i] << endl;
				cout << "The number of polynomials is " << parameter1a[1][j] << endl;
				cout << "The maturity is " << parameter1a[2][h] << endl;
				cout << "The american put option price using Lagerrue is ";
				generate_simulations(x,simulation,option_value,parameter1a[0][i],0.2,0.06,10000,40,100, parameter1a[2][h]);
				question1L("Laguerre",x,simulation, option_value,10000,100, parameter1a[1][j], parameter1a[2][h]);
				cout << endl;
			}
		}
	}
	/*Question 1b*/
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int h = 0; h < 3; h++) {
				cout << "The current price is " << parameter1a[0][i] << endl;
				cout << "The number of polynomials is " << parameter1a[1][j] << endl;
				cout << "The maturity is " << parameter1a[2][h] << endl;
				cout << "The american put option price using Hermite Polynomials is ";
				generate_simulations(x,simulation,option_value,parameter1a[0][i],0.2,0.06,50000,40,100, parameter1a[2][h]);
				question1L("Hermite",x,simulation, option_value,50000,100, parameter1a[1][j], parameter1a[2][h]);
				cout << endl;
			}
		}
	}
	/*Question 1c*/
	
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int h = 0; h < 3; h++) {
				cout << "The current price is " << parameter1a[0][i] << endl;
				cout << "The number of polynomials is " << parameter1a[1][j] << endl;
				cout << "The maturity is " << parameter1a[2][h] << endl;
				cout << "The american put option price using Monominal is ";
				generate_simulations(x,simulation,option_value,parameter1a[0][i],0.2,0.06,50000,40,100, parameter1a[2][h]);
				question1L("Monomial",x,simulation, option_value,50000,100, parameter1a[1][j], parameter1a[2][h]);
				cout << endl;
			}
		}
	}
	
	
	question2a();
	cout << endl;
	question2b();
	system("pause");
	return 0;
}