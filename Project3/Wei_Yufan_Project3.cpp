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

double european_call_option_price_mts(double r, double sigma, double S0, double T, double X, int n = 10000000, int seed = 86534) {
	double total = 0;
	vector<double> UR;//uniform
	vector<double> NR;//normal
	generate_uniform_distributed_sequence(UR, n, seed);
	generate_standard_normal_distributed_sequence(NR, UR);//generate normal sequence
	for (int i = 0; i < n; i++) {
		double e = sigma*NR[i]*sqrt(T) + (r - 0.5*sigma*sigma)*T;
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
		double e1 = sigma*NR[i]*sqrt(T) + (r - 0.5*sigma*sigma)*T;
		double e = sigma*(-NR[i])*sqrt(T) + (r - 0.5*sigma*sigma)*T;
		double St = S0*exp(e);
		double St1 = S0*exp(e1);
		double cp = ((St - X) > 0) ? (St - X) : 0;
		double cm = ((St1 - X) > 0) ? (St1 - X) : 0;
		double discount_factor = -r*T;
		cp = cp*exp(discount_factor); cm = cm*exp(discount_factor);
		total = total + (cp+cm)/2.0;
	}
	return total/ n;
}
double european_call_option_price_BSM(double r, double sigma, double S0, double T, double X) {
	double d1 = 1.0 / (sigma*sqrt(T))*(log(S0 / X) + (r + sigma*sigma / 2.0)*T);
	double d2 = d1 - sigma*sqrt(T);
	double Nd1,Nd2;
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

double call_price_Heston(string method,double T = 0.5, double X = 50.0, double tho = -0.6,double r = 0.03,
	double S0 = 48.0, double V0 = 0.05,double sigma = 0.42, double a = 5.8, double b = 0.0625, double n = 4000) {
	double delta = T / n; double sum = 0; double discount_factor = -r*T;
	// take the method as an input, decide which method should be used
	if(method=="Full Truncation"){
		for (int i = 0; i < n; i++) {
			double St = S0; double Vt = V0;
			vector<double> UR;//uniform
			vector<double> NR;//normal
			int seed =abs(generate_normal_number(1000,300))+i; // scramble seed
			generate_uniform_distributed_sequence(UR, n, seed); //generate normal sequence
			generate_standard_normal_distributed_sequence(NR, UR);
			vector<double> UR1;//uniform
			vector<double> NR1;//normal
			seed = abs(generate_normal_number(1000, 300))+11+i;
			generate_uniform_distributed_sequence(UR1, n, seed);
			generate_standard_normal_distributed_sequence(NR1, UR1);
			for (int k = 0; k < n; k++) {
				St = St + r*St*delta + sqrt((Vt>0?Vt:0))*St*sqrt(delta)*NR[k];
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
			int seed = abs(generate_normal_number(1000, 300))+(2*i)%17;
			generate_uniform_distributed_sequence(UR, n, seed);
			generate_standard_normal_distributed_sequence(NR, UR);
			vector<double> UR1;//uniform
			vector<double> NR1;//normal
			seed = abs(generate_normal_number(1000, 300))+i%2;
			generate_uniform_distributed_sequence(UR1, n, seed);
			generate_standard_normal_distributed_sequence(NR1, UR1);
			for (int k = 0; k < n; k++) {
				St = St + r*St*delta + sqrt((Vt>0 ? Vt : 0))*St*sqrt(delta)*NR[k];
				Vt = Vt + a*(b - Vt)*delta + sigma*sqrt((Vt>0?Vt:0))*sqrt(delta)*(tho*NR[k] + sqrt(1.0 - tho*tho)*NR1[k]);
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
			int seed = abs(generate_normal_number(1000, 300))+(i*5)%17;
			generate_uniform_distributed_sequence(UR, n, seed);
			generate_standard_normal_distributed_sequence(NR, UR);
			vector<double> UR1;//uniform
			vector<double> NR1;//normal
			seed = abs(generate_normal_number(1000, 300))+(i*6)%13;
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
double Halton_integral( int base1, int base2, int n = 10000) {
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

	return sum/n;
}

/*Question 1*/
int main() {	
	
	double X0 = 1;
	double Y0 = 0.75;
	double t = 2.0; int n = 2000; double delta = t / n;//number of interval and simulation
	int seed = 154;
	/*dy=2/(1+t)ydt+(1+t^3)/3dt+(1+t^3)/3dZ*/
	
	//simulate 10000 times
	int count = 0;
	for (int k = 0; k < n; k++) {
	double Yk = Y0;
		vector<double> UR;//uniform
		vector<double> NR;//normal
		seed += 11;
		generate_uniform_distributed_sequence(UR, n, seed);
		generate_standard_normal_distributed_sequence(NR, UR);
		for (int i = 0; i < n; i++) {
			Yk += (2.0 / (1 + i*delta)*Yk+(1+ delta*i*delta*i*delta*i)/3.0)*delta + (1 + delta*i*delta*i*delta*i) / 3.0 *sqrt(delta)* NR[i];
		}
		if (Yk > 5)
			count++;
	}
	double p = 100.00*count / n;
	cout << "P(Y2>5) = " << p<< "%"<<endl;

	
	long double sum = 0;
	for (int k = 0; k < n; k++) {
		double Xk = 1.0;
		vector<double> UR;//uniform
		vector<double> NR;//normal
		seed +=1;
		generate_uniform_distributed_sequence(UR, n, seed);
		generate_standard_normal_distributed_sequence(NR, UR);
		for (int i = 0; i < n; i++) {
			Xk += (0.2 - 0.5*Xk)*delta + (2.0 / 3.0) * sqrt(delta)*NR[i];
		}
		sum += cbrt(Xk);
	}
	double Ex = sum / n;
	cout << Ex << endl;
	//EY3
	t = 3;
	delta = t / n;
	sum = 0; seed = 0;
	for (int k = 0; k < n; k++) {
		double Yk = Y0;
		vector<double> UR;//uniform
		vector<double> NR;//normal
		seed += 1;
		generate_uniform_distributed_sequence(UR, n, seed);
		generate_standard_normal_distributed_sequence(NR, UR);
		for (int i = 0; i < n; i++) {
			Yk += (2.0 / (1 + i*delta)*Yk + (1 + delta*i*delta*i*delta*i) / 3.0)*delta + (1 + delta*i*delta*i*delta*i) / 3.0 *sqrt(delta)* NR[i];
		}
		sum +=  Yk;
	}
	double EY3 = sum / n;
	cout << EY3 << endl;
	bool x2;
	sum = 0; t = 2; delta = t / n;
	for (int k = 0; k < n; k++) {
		double Xk = 1.0;
		double Yk = Y0;
		vector<double> UR;//uniform
		vector<double> NR;//normal
		seed += 1;
		generate_uniform_distributed_sequence(UR, n, seed);
		generate_standard_normal_distributed_sequence(NR, UR);
		vector<double> UR1;//uniform
		vector<double> NR1;//normal
		seed += 2;
		generate_uniform_distributed_sequence(UR1, n, seed);
		generate_standard_normal_distributed_sequence(NR1, UR1);
		for (int i = 0; i < n; i++) {
			Xk += (0.2-0.5*Xk)*delta + (2.0 / 3.0) * sqrt(delta)*NR[i];
			Yk += (2.0 / (1 + i*delta)*Yk + (1.0 + delta*i*delta*i*delta*i) / 3.0)*delta + (1.0 + delta*i*delta*i*delta*i) / 3.0 * sqrt(delta)* NR1[i];
		}
		x2 = Xk > 1;
		sum += (Xk)*Yk*x2;
	}	
	double Ex2y2 = sum / n;
	cout << Ex2y2 << endl;
	/*Question 2*/
	t = 3; delta = t / n;
	X0 = 1.0; sum = 0;
	Y0 = 1.0; seed = 123;
	for (int i = 0; i < n; i++) {
		double Xt = X0;
		vector<double> UR;//uniform
		vector<double> NR;//normal
		seed += 1;
		generate_uniform_distributed_sequence(UR, n, seed);
		generate_standard_normal_distributed_sequence(NR, UR);
		for (int k = 0; k < n; k++) {

			Xt += 0.25*Xt*delta + 1.0 / 3.0 * Xt*sqrt(delta)*NR[k] - 0.75*Xt*sqrt(delta)*NR[k];
		}
		sum += cbrt(Xt+1);
	}
	double E3 = sum / n;
	cout << E3 << endl;

	delta = t / n; sum = 0;
	seed = 0;
	for (int i = 0; i < n; i++) {
		vector<double> UR;//uniform
		vector<double> NR;//normal
		seed +=  11;
		generate_uniform_distributed_sequence(UR, n, seed);
		generate_standard_normal_distributed_sequence(NR, UR);
			double e = -0.08*t + 1.0 / 3.0 * NR[i]*sqrt(t) + 3.0 / 4.0 * NR[i]*sqrt(t);
			double Yt = exp(e);
			sum += cbrt(Yt+1);
		}
	
	double Y3 = sum / n;
	cout << Y3 << endl;
	/*Question 3*/
	//a) The function is on the top and also in the header
	//b) The function is on the top and also in the header
	//c)
	
	double X = 20; double sigma = 0.25; double r = 0.04; double T = 0.5;
	double S0 = 15;
	double call_price;
	call_price=european_call_option_price_BSM(r, sigma, S0, T, X);
	cout << call_price<<endl;
	call_price = european_call_option_price_vr(r, sigma, S0, T, X);
	cout << call_price << endl;
	call_price = european_call_option_price_mts(r, sigma, S0, T, X);
	cout << call_price << endl;
		double eps = 0.01;//set the change in each variable as 1%
	double Delta[11]; double Vega[11]; double Theta[11]; double Rho[11]; double Gamma[11];
	for (int i = 15; i <= 25; i++) {
		S0 = i; double S1 = S0 *(1+eps); double S2 = S0 *(1- eps); double r1 = r *(1+eps); double r2 = r *(1- eps);
		double sigma1 = sigma *(1+ eps); double sigma2 = sigma *(1- eps); double T1 = T *(1+ eps); double T2 = T *(1- eps);
		// initialize the change and use variance reduction method to derive call prices after the change
		double call_price1 = european_call_option_price_vr(r, sigma, S1, T, X);
		double call_price2 = european_call_option_price_vr(r, sigma, S2, T, X);
		call_price = european_call_option_price_vr(r, sigma, S0, T, X);

		Gamma[i - 15] = (call_price1 + call_price2 - 2.0 * call_price) / (S0*eps*eps*S0);
		Delta[i - 15] = (call_price1 - call_price2) / (2.0 * eps*S0);
		
		call_price1 = european_call_option_price_vr(r, sigma1, S0, T, X);
		call_price2 = european_call_option_price_vr(r, sigma2, S0, T, X);
		Vega[i-15]= (call_price1 - call_price2) / (2.0 * eps*sigma);
		
		call_price1 = european_call_option_price_vr(r, sigma, S0, T1, X);
		call_price2 = european_call_option_price_vr(r, sigma, S0, T2, X);
		Theta[i - 15] = (call_price1 - call_price2) / (2.0 * eps*T);

		call_price1 = european_call_option_price_vr(r1, sigma, S0, T, X);
		call_price2 = european_call_option_price_vr(r2, sigma, S0, T, X);
		Rho[i - 15] = (call_price1 - call_price2) / (2.0 * eps*r);		
	}
	cout << Delta << endl; cout << Gamma << endl; cout << Vega << endl; cout << Theta << endl; cout << Rho << endl;
	/*Question 4*/
	T = 0.5; X = 50; double tho = -0.6; r = 0.03;
	S0 = 48; double V0 = 0.05; sigma = 0.42; double a = 5.8; double b = 0.0625;
	n = 10000; 
	cout<<call_price_Heston("Full Truncation")<<endl;
	cout<<call_price_Heston("Partial Truncation")<<endl;
	cout << call_price_Heston("Reflection")<<endl;
	
	/*Question 5*/
	//a)
	vector<double> RV1;
	vector<double> RV2;
	generate_uniform_distributed_sequence(RV1, 100, 423);
	generate_uniform_distributed_sequence(RV2, 100, 1234);
	std::ofstream myfile3("question5_uni1.txt");
	if (myfile3.is_open()) {
		for (int count = 0; count < 100; count++) {
			myfile3 << RV1[count];
			if (count != 100 - 1)
				myfile3 << "" << endl;
		}
		myfile3 << '\n';
		myfile3.close();
	}
	std::ofstream myfile4("question5_uni2.txt");
	if (myfile4.is_open()) {
		for (int count = 0; count < 100; count++) {
			myfile4 << RV2[count];
			if (count != 100 - 1)
				myfile4 << "" << endl;
		}
		myfile4 << '\n';
		myfile4.close();
	}
	//b)
	double QMC1[100];
	double QMC2[100];
	for (int i = 0; i < 100; i++) {
		QMC1[i] = halton_sequence(i+1, 2);
	}
	for (int i = 0; i < 100; i++) {
		QMC2[i] = halton_sequence(i + 1, 7);
	}
	std::ofstream myfile("question5_1.txt");
	if (myfile.is_open()) {
		for (int count = 0; count < 100; count++) {
			myfile << QMC1[count];
			if (count != 100 - 1)
				myfile << "" << endl;
		}
		myfile << '\n';
		myfile.close();
	}
	std::ofstream myfile2("question5_2.txt");
	if (myfile2.is_open()) {
		for (int count = 0; count < 100; count++) {
			myfile2 << QMC2[count];
			if (count != 100 - 1)
				myfile2 << "" << endl;
		}
		myfile2 << '\n';
		myfile2.close();
	}
	//c)
	double QMC3[100];
	double QMC4[100];
	for (int i = 0; i < 100; i++) {
		QMC3[i] = halton_sequence(i + 1, 2);
	}
	for (int i = 0; i < 100; i++) {
		QMC4[i] = halton_sequence(i + 1, 4);
	}
	std::ofstream myfile5("question5_3.txt");
	if (myfile5.is_open()) {
		for (int count = 0; count < 100; count++) {
			myfile5 << QMC3[count];
			if (count != 100 - 1)
				myfile5 << "" << endl;
		}
		myfile5 << '\n';
		myfile5.close();
	}
	std::ofstream myfile6("question5_4.txt");
	if (myfile6.is_open()) {
		for (int count = 0; count < 100; count++) {
			myfile6 << QMC4[count];
			if (count != 100 - 1)
				myfile6 << "" << endl;
		}
		myfile6 << '\n';
		myfile6.close();
	}
	//d)
	double I;
	I = Halton_integral(2, 4);
	cout << I << endl;

	I = Halton_integral(2, 7);
	cout << I << endl;

	I = Halton_integral(5, 7);
	cout << I << endl;

	//output the data above to txt files, plot in R
/*	std::ofstream myfile("question3_1.txt");
	if (myfile.is_open()) {
		for (int count = 0; count < 11; count++) {
			myfile << Delta[count];
			if (count != 11 - 1)
				myfile << "" << endl;
		}
		myfile << '\n';
		myfile.close();
	}
	std::ofstream myfile2("question3_2.txt");
	if (myfile2.is_open()) {
		for (int count = 0; count < 11; count++) {
			myfile2 << Gamma[count];
			if (count != 11 - 1)
				myfile2 << "" << endl;
		}
		myfile2 << '\n';
		myfile2.close();
	}
	std::ofstream myfile3("question3_3.txt");
	if (myfile3.is_open()) {
		for (int count = 0; count < 11; count++) {
			myfile3 << Theta[count];
			if (count != 11 - 1)
				myfile3 << "" << endl;
		}
		myfile3 << '\n';
		myfile3.close();
	}
	std::ofstream myfile4("question3_4.txt");
	if (myfile4.is_open()) {
		for (int count = 0; count < 11; count++) {
			myfile4 << Vega[count];
			if (count != 11 - 1)
				myfile4 << "" << endl;
		}
		myfile4 << '\n';
		myfile4.close();
	}
	std::ofstream myfile5("question3_5.txt");
	if (myfile5.is_open()) {
		for (int count = 0; count < 11; count++) {
			myfile5 << Rho[count];
			if (count != 11 - 1)
				myfile5 << "" << endl;
		}
		myfile5 << '\n';
		myfile5.close();
	}
	*/
	system("pause");
	return 0;
}