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
#include <chrono>
#include <random>
using namespace std;
// create a class of node, containing current stock price, the up and down magnitude and the probability of going up.
class node { 
private:
	double StockPrice; double strike_price;
	double up; double down; double Probability; double option_value;
	node *next_up; node *next_down;
public:
	node(double s, double k, double u, double d, double p) {
		StockPrice = s; strike_price = k; up = u; down = d; Probability = p;
		next_up = NULL; next_down = NULL;	
	}
	void call_value() {
		option_value = (StockPrice > strike_price) ? (StockPrice - strike_price) : 0;
	}
	void put_value() {
		option_value = (StockPrice < strike_price) ? (strike_price- StockPrice) : 0;
	}

};

double binomial_option_price(const double &S, const double &K, const double &sigma, const double &r, const double T, const int& n, string question, bool is_call = 1) {
	double delta = T / n; double c, p, d, u;
	double discount_factor = exp(-r*delta);
	if (question == "a") {
		c = 0.5*(exp(-r*delta) + exp((r + sigma*sigma)*delta));
		d = c - sqrt(c*c - 1);
		u = 1 / d;
		p = (exp(r*delta) - d) / (u - d);
	}
	else if (question == "b") {
		d = exp(r*delta)*(1 - sqrt(exp(sigma*sigma*delta) - 1));
		u = exp(r*delta)*(1 + sqrt(exp(sigma*sigma*delta) - 1));
		p = 0.5;
	}
	else if (question == "c") {
		u = exp((r - sigma*sigma*0.5)*delta + sigma*sqrt(delta));
		d = exp((r - sigma*sigma*0.5)*delta - sigma*sqrt(delta));
		p = 0.5;
	}
	else if (question == "d") {
		u = exp(sigma*sqrt(delta));
		d = exp(-sigma*sqrt(delta));
		p = 0.5 + 0.5*((r - sigma*sigma*0.5)*sqrt(delta) / sigma);
	}
	else
		cout << "re-enter the question number" << endl;

	double ud = u / d;
	vector<double> price(n + 1);
	price[0] = S*pow(d, n);
	for (int i = 1; i <= n; ++i) {
		price[i] = ud*price[i - 1];
	}

	vector<double> call_value(n + 1);
	if (is_call == 1){
		for (int i = 0; i <= n; ++i) {
			call_value[i] = max(0.0, price[i] - K);
		}
}
	else {
		for (int i = 0; i <= n; ++i) {
			call_value[i] = max(0.0, K-price[i]);
		}
	}

	for (int step = n - 1; step >= 0; --step) {
		for (int j = 0; j <= step; ++j) {
			call_value[j] = (p*call_value[j + 1] + (1 - p)*call_value[j])*discount_factor;
		}
	}
	return call_value[0];
}

double binomial_american_option_price(const double &S, const double &K, const double &sigma, const double &r, const double T, const int& n, string question, bool is_call = 1) {
	double delta = T / n; double c, p, d, u;
	double discount_factor = exp(-r*delta);
	if (question == "a") {
		c = 0.5*(exp(-r*delta) + exp((r + sigma*sigma)*delta));
		d = c - sqrt(c*c - 1);
		u = 1 / d;
		p = (exp(r*delta) - d) / (u - d);
	}
	else if (question == "b") {
		d = exp(r*delta)*(1 - sqrt(exp(sigma*sigma*delta) - 1));
		u = exp(r*delta)*(1 + sqrt(exp(sigma*sigma*delta) - 1));
		p = 0.5;
	}
	else if (question == "c") {
		u = exp((r - sigma*sigma*0.5)*delta + sigma*sqrt(delta));
		d = exp((r - sigma*sigma*0.5)*delta - sigma*sqrt(delta));
		p = 0.5;
	}
	else if (question == "d") {
		u = exp(sigma*sqrt(delta));
		d = exp(-sigma*sqrt(delta));
		p = 0.5 + 0.5*((r - sigma*sigma*0.5)*sqrt(delta) / sigma);
	}
	else
		cout << "re-enter the question number" << endl;

	double ud = u / d;
	vector<double> price(n + 1);
	price[0] = S*pow(d, n);
	for (int i = 1; i <= n; ++i) {
		price[i] = ud*price[i - 1];
	}

	vector<double> call_value(n + 1);
	if (is_call == 1) {
		for (int i = 0; i <= n; ++i) {
			call_value[i] = max(0.0, price[i] - K);
		}
		for (int step = n - 1; step >= 0; --step) {
			for (int j = 0; j <= step; ++j) {
				call_value[j] = (p*call_value[j + 1] + (1 - p)*call_value[j])*discount_factor;
				price[j] = price[j] / d;
				price[j] = max(call_value[j], price[j]-K);
			}
		}
	}
	else {
		for (int i = 0; i <= n; ++i) {
			call_value[i] = max(0.0, K - price[i]);
		}
		for (int step = n - 1; step >= 0; --step) {
			for (int j = 0; j <= step; ++j) {
				call_value[j] = (p*call_value[j + 1] + (1 - p)*call_value[j])*discount_factor;
				price[j] = price[j] / d;
				call_value[j] = max(call_value[j], K -price[j]);
			}
		}
	}


	return call_value[0];
}

double trinomial_option_price(const double &S, const double &K, const double &sigma, const double &r, const double T, const int& n, string question, bool is_call = 1) {
	double delta = T / n; double c, pu,pm,pd, d, u;
	double discount_factor = exp(-r*delta);
	vector<double> price(2 * n + 1);
	vector<double> call_value(2*n + 1);
	if (question == "a") {
		d = exp(-sigma*sqrt(3 * delta));
		u = 1 / d;
		pd = (r*delta*(1 - u) + (r*delta)*(r*delta) + sigma*sigma*delta) / ((u - d)*(1 - d));
		pu= (r*delta*(1 - d) + (r*delta)*(r*delta) + sigma*sigma*delta) / ((u - d)*(u-1));
		pm = 1 - pu - pd;
		double ud = u / d;
		
		price[0] = S*pow(d, n);
		for (int i = 1; i < 2*n+1; i++){
			price[i] = price[i - 1]/d;
		}

		
		if (is_call == 1) {
			for (int i = 0; i < 2*n+1; i++) {
				call_value[i] = max(0.0, price[i] - K);
			}
		}
		else {
			for (int i = 0; i < 2*n+1; i++) {
				call_value[i] = max(0.0, K - price[i]);
			}
		}

		for (int step = n - 1; step >= 0; --step) {
			for (int j = 0; j < 2*step+1; j++) {
				call_value[j] = (pu*call_value[j + 2] + pd*call_value[j]+pm*call_value[j+1])*discount_factor;
			}
		}
	}
	else if (question == "b") {
		d = -sigma*sqrt(3*delta);
		u = -d; double X = log(S);
		pd=0.5*((sigma*sigma*delta+(r-sigma*sigma*0.5)*(r - sigma*sigma*0.5)*delta*delta)/(u*u)-(r-sigma*sigma*0.5)*delta/u);
		pu= 0.5*((sigma*sigma*delta + (r - sigma*sigma*0.5)*(r - sigma*sigma*0.5)*delta*delta) / (u*u) + (r - sigma*sigma*0.5)*delta / u);
		pm = 1 - pu - pd;
		price[0] = X+n*d;
		for (int i = 1; i < 2*n+1; i++) {
			price[i] = u+price[i - 1];
		}

		if (is_call == 1) {
			for (int i = 0; i < 2*n+1; i++) {
				call_value[i] = max(0.0, exp(price[i]) - K);
			}
		}
		else {
			for (int i = 0; i < 2*n+1; i++) {
				call_value[i] = max(0.0, K - exp(price[i]));
			}
		}

		for (int step = n - 1; step >= 0; --step) {
			for (int j = 0; j <2* step+1; j++) {
				call_value[j] = (pu*call_value[j + 2] + pd*call_value[j] + pm*call_value[j + 1])*discount_factor;
			}
		}
	}

	else
		cout << "re-enter the question number" << endl;


	return call_value[0];
}

double Question6(double S0, double K, double T, double r, double sigma, double N = 10000, double b1 = 1, double b2 = 2) {

	vector<double> H1, H2; int N6 = N;
	for (int i = 0; i < N6; i++) {
		H1.push_back(halton_sequence(i, b1));
		H2.push_back(halton_sequence(i, b2));
	}
	vector<double>Z1, Z2;
	for (int i = 0; i < N6; i++) {
		Z1.push_back(sqrt(-2.0*log(H1[i])*cos(2 * H2[i] * (atan(1) * 4))));
		Z2.push_back(sqrt(-2.0*log(H1[i])*sin(2 * H2[i] * (atan(1) * 4))));
	}
	double sum = 0;
	for (int i = 0; i < N6; i++) {
		sum += max(0.0, S0*exp((r - sigma*sigma*0.5)*T + sigma*sqrt(T)*Z1[i]) - K) + max(0.0, S0*exp((r - sigma*sigma*0.5)*T + sigma*sqrt(T)*Z2[i]) - K);
	}
	double call_price = sum / 2 * N6;
	return call_price;
}

int main() {
	/*Question 1*/
	
	double r = 0.05; double sigma = 0.24; double S0 = 32.0; double K = 30.0; double T = 0.5;
	double n[7] = { 10,20,40,80,100,200,500 };
	for (int i = 0; i < 7; i++) {
		int N = n[i];
		double call_price_a = binomial_option_price(S0, K, sigma, r, T, N,"a");
		cout << "The call option price for question a is " << call_price_a << " when n=" << N << endl;
	}
	for (int i = 0; i < 7; i++) {
		int N = n[i];
		double call_price_b = binomial_option_price(S0, K, sigma, r, T, N, "b");
		cout << "The call option price for question b is " << call_price_b << " when n=" << N << endl;
	}
	for (int i = 0; i < 7; i++) {
		int N = n[i];
		double call_price_c = binomial_option_price(S0, K, sigma, r, T, N, "c");
		cout << "The call option price for question c is " << call_price_c << " when n=" << N << endl;
	}
	for (int i = 0; i < 7; i++) {
		int N = n[i];
		double call_price_d = binomial_option_price(S0, K, sigma, r, T, N, "d");
		cout << "The call option price for question d is " << call_price_d << " when n=" << N << endl;
	}

	/*Question 2*/
	//a)
	
	r = 0.02; K = 1205; S0 = 1145; sigma = 0.23465627965146882;
	T = 1;
	cout<<"The estimated option price is "<<binomial_option_price(S0, K, sigma, r, T, 1000, "a")<<endl;
	//b)
	//The current price of this option is 87.90, so the implied volatility should be smaller than the realized one.
	while (binomial_option_price(S0, K, sigma, r, T, 1000, "a") - 87.90>0.0001) {
		sigma = sigma - 0.0001;
	}
	cout << "The implied sigma is " << sigma << endl;

	/*Question 3*/
	S0 = 49; K = 50; r = 0.03; sigma = 0.2; T = 0.3846;
	double eps = 0.01;
	//i)
	vector<double> Delta; vector<double> Theta; vector<double> Gamma;
	vector<double> Vega; vector<double> Rho;
	for (int i = 20; i <= 80; i=i + 2) {
		double a = (binomial_option_price(i*(1 + eps), K, sigma, r, T, 500, "a") - binomial_option_price(i*(1 - eps), K, sigma, r, T, 500, "a")) / (2 * eps*i);
		Delta.push_back(a);
		double b = (binomial_option_price(i, K, sigma, r, T*(1 + eps), 500, "a") - binomial_option_price(i, K, sigma, r, T*(1 - eps), 500, "a")) / (2 * eps*T);
		Theta.push_back(b);
		double c = (binomial_option_price(i, K, sigma*(1 + eps), r, T, 500, "a") - binomial_option_price(i, K, sigma*(1 - eps), r, T, 500, "a")) / (2 * eps*sigma);
		Vega.push_back(c);
		double d = (binomial_option_price(i, K, sigma, r*(1 + eps), T, 500, "a") - binomial_option_price(i, K, sigma, r*(1 - eps), T, 500, "a")) / (2 * eps*r);
		Rho.push_back(d);
		double e = (binomial_option_price(i*(1 + eps), K, sigma, r, T, 500, "a") - 2*binomial_option_price(i, K, sigma, r, T, 500, "a") + binomial_option_price(i*(1 - eps), K, sigma, r, T, 500, "a")) / (eps*i * eps*i);
		Gamma.push_back(e);
	}
	std::ofstream myfile("question3_delta.txt");
	if (myfile.is_open()) {
		for (int count = 0; count < Delta.size(); count++) {
			myfile << Delta[count];
			if (count != Delta.size() - 1)
				myfile << "" << endl;
		}
		myfile << '\n';
		myfile.close();
	}
	std::ofstream myfile1("question3_theta.txt");
	if (myfile1.is_open()) {
		for (int count = 0; count < Theta.size(); count++) {
			myfile1 << Theta[count];
			if (count != Theta.size() - 1)
				myfile1 << "" << endl;
		}
		myfile1 << '\n';
		myfile1.close();
	}
	std::ofstream myfile2("question3_vega.txt");
	if (myfile2.is_open()) {
		for (int count = 0; count < Vega.size(); count++) {
			myfile2 << Vega[count];
			if (count != Vega.size() - 1)
				myfile2 << "" << endl;
		}
		myfile2 << '\n';
		myfile2.close();
	}
	std::ofstream myfile3("question3_rho.txt");
	if (myfile3.is_open()) {
		for (int count = 0; count < Rho.size(); count++) {
			myfile3 << Rho[count];
			if (count != Rho.size() - 1)
				myfile3 << "" << endl;
		}
		myfile3 << '\n';
		myfile3.close();
	}
	std::ofstream myfile4("question3_gamma.txt");
	if (myfile4.is_open()) {
		for (int count = 0; count < Gamma.size(); count++) {
			myfile4 << Gamma[count];
			if (count != Gamma.size() - 1)
				myfile4 << "" << endl;
		}
		myfile4 << '\n';
		myfile4.close();
	}
	vector<double> Delta_t;
	for (int i = 0; i <= 38; i++) {
		double a = (binomial_option_price(S0*(1 + eps), K, sigma, r, i/100.0, 500, "a") - binomial_option_price(S0*(1 - eps), K, sigma, r, i/100.0, 500, "a")) / (2 * eps*S0);
		Delta_t.push_back(a);
	}
	std::ofstream myfile5("question3_delta_t.txt");
	if (myfile5.is_open()) {
		for (int count = 0; count < Delta_t.size(); count++) {
			myfile5 << Delta_t[count];
			if (count != Delta_t.size() - 1)
				myfile5 << "" << endl;
		}
		myfile5 << '\n';
		myfile5.close();
	}
	/*Question 4*/
	T = 1; r = 0.05; sigma = 0.3; K = 100;
	vector<double> XYZ_E_put, XYZ_A_put;
	for (int S = 80; S <= 120; S=S + 4) {
		XYZ_E_put.push_back(binomial_option_price(S, K, sigma, r, T, 500, "a", false));
		XYZ_A_put.push_back(binomial_american_option_price(S, K, sigma, r, T, 500, "a", false));
	}
	std::ofstream myfile6("question4_1.txt");
	if (myfile6.is_open()) {
		for (int count = 0; count < XYZ_E_put.size(); count++) {
			myfile6 << XYZ_E_put[count];
			if (count != XYZ_E_put.size() - 1)
				myfile6 << "" << endl;
		}
		myfile6 << '\n';
		myfile6.close();
	}
	std::ofstream myfile7("question4_2.txt");
	if (myfile7.is_open()) {
		for (int count = 0; count < XYZ_A_put.size(); count++) {
			myfile7 << XYZ_A_put[count];
			if (count != XYZ_A_put.size() - 1)
				myfile7 << "" << endl;
		}
		myfile7 << '\n';
		myfile7.close();
	}
	/*Question 5*/
	//a)
	
	K = 30; r = 0.05; sigma = 0.24; S0 = 32; T = 0.5;
	int n5[9] = { 10,15,20,40,70,80,100,200,500 };
	for (int i = 0; i < 9; i++) {
		int N = n5[i];
		cout<<trinomial_option_price(S0, K, sigma, r, T, N, "a", true)<<endl;
	}
	//b)
	for (int i = 0; i < 9; i++) {
		int N = n5[i];
		cout << trinomial_option_price(S0, K, sigma, r, T, N, "b", true) << endl;
	}
	
	/*Question 6*/
	//The function is above named Question 6
	/*double Question6(double S0, double K, double T, double r, double sigma, double N = 10000, double b1 = 1, double b2 = 2) {

	vector<double> H1, H2; int N6 = N;
	for (int i = 0; i < N6; i++) {
		H1.push_back(halton_sequence(i, b1));
		H2.push_back(halton_sequence(i, b2));
	}
	vector<double>Z1, Z2;
	for (int i = 0; i < N6; i++) {
		Z1.push_back(sqrt(-2.0*log(H1[i])*cos(2 * H2[i] * (atan(1) * 4))));
		Z2.push_back(sqrt(-2.0*log(H1[i])*sin(2 * H2[i] * (atan(1) * 4))));
	}
	double sum = 0;
	for (int i = 0; i < N6; i++) {
		sum += max(0.0, S0*exp((r - sigma*sigma*0.5)*T + sigma*sqrt(T)*Z1[i]) - K) + max(0.0, S0*exp((r - sigma*sigma*0.5)*T + sigma*sqrt(T)*Z2[i]) - K);
	}
	double call_price = sum / 2 * N6;
	return call_price;
}*/
	system("Pause");
	return 0;
}
