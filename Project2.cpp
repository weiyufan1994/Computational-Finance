#include<iostream>
#include<vector>
#include<math.h>
#include<time.h>
#include<cstdlib>
#include<stdlib.h>
#include<array>
#include<cmath>
#include<numeric>
#include<cmath>
#include<fstream>
#include"Header.h"
using namespace std;
double Q3_a_simulation(vector<double> &standard_normal_sequence) {
	double Ea = 0;
	for (int i = 0; i < standard_normal_sequence.size(); i++) {
		Ea += standard_normal_sequence[i] * standard_normal_sequence[i] * 5 + sin(standard_normal_sequence[i] * sqrt(5));
	}
	return Ea/standard_normal_sequence.size();
}
double Q3_a_simulation(vector<double> &standard_normal_sequence, double t) {
	double Ea = 0;
	double e = t / 2;
	for (int i = 0; i < standard_normal_sequence.size(); i++) {
		Ea += exp(e)*cos(standard_normal_sequence[i] * sqrt(t));
	}
	return Ea / standard_normal_sequence.size();
}

//antithetical method
double Q3_b_simulation(vector<double> &standard_normal_sequence) {
	double Ea = 0;
	for (int i = 0; i < standard_normal_sequence.size(); i++) {
		Ea += 0.5*(standard_normal_sequence[i] * standard_normal_sequence[i] * 10 + sin(standard_normal_sequence[i] * sqrt(5)) + sin((-standard_normal_sequence[i]) * sqrt(5)));
	}
	return Ea / standard_normal_sequence.size();
}
//antithetical method
double Q3_b_simulation(vector<double> &standard_normal_sequence, double t) {
	double Ea = 0;
	double e = t / 2;
	for (int i = 0; i < standard_normal_sequence.size(); i++) {
		Ea += 0.5*exp(e)*cos(standard_normal_sequence[i] * sqrt(t))+ 0.5*exp(e)*cos((-standard_normal_sequence[i]) * sqrt(t));
	}
	return Ea / standard_normal_sequence.size();
}

double european_call_option_price_mts(double r, double sigma,double S0, double T,double X,int n=10000) {
	double total = 0;
	for (int i = 0; i < n; i++) {
		double e = sigma*generate_normal_number()*sqrt(T) + (r - 0.5*sigma*sigma)*T;
		double St = S0*exp(e);
		if (St > X) {
			total += St-X;
		}
		else
			total = total;
	}
	double discount_factor = -r*T;
	return total*exp(discount_factor) / n;
}
double european_call_option_price_vr(double r, double sigma, double S0, double T, double X, int n = 10000) {
	double total = 0;
	for (int i = 0; i < n; i++) {
		double random=generate_normal_number();
		double e1 = sigma*random*sqrt(T) + (r - 0.5*sigma*sigma)*T;
		double e = sigma*(-random)*sqrt(T) + (r - 0.5*sigma*sigma)*T;
		double St = S0*exp(e);
		double St1 = S0*exp(e1);
		if (St > X) {
			St = St - X;
		}
		else
			St = 0;
		if (St1 > X) {
			St1 = St1 - X;
		}
		else
			St1 = 0;
		
		total = total+St*0.5+St1*0.5;
	}
	double discount_factor = -r*T;
	return total*exp(discount_factor) / n;
}



int main() {
	int seed=16573;
	//Question 1
	// Generate a sequence of uniform distributed sequence
	vector<double> uniform_sequence;
	int n = 1000;
	double a = -0.7;
	//int seed=1;
	//take seed for distribution, n times, and covariance a.
	cout << "Enter an number as the seed for distribution, n and a respectively" << endl;
	//cin >> seed;
	//cin >> n;
	//cin >> a;
	generate_uniform_distributed_sequence(uniform_sequence, 2*n,seed);//Generate 2000 uniform distributed random numbers

	vector<double> standard_normal_sequence;//initiate a sequence used to store normal distributed numbers
	generate_standard_normal_distributed_sequence(standard_normal_sequence, uniform_sequence);// n random numbers will be generated
	//Since the mean vector is (0,0) and diag of covariance matrix is (1,1), x and y are standard normal distributed
	vector<double> X, Y;

	for (int i = 0; i < n; i++) {
		X.push_back(standard_normal_sequence[i]);
		Y.push_back(standard_normal_sequence[i] * (a) + standard_normal_sequence[i + n] * sqrt(1 - a*a));
	}
	double X_sum=0, X_mean = 0, X_std = 0;
	double Y_sum = 0, Y_mean = 0, Y_std = 0;
	calculate_sum_mean_stdev(X, X_sum, X_mean, X_std);
	calculate_sum_mean_stdev(Y, Y_sum, Y_mean, Y_std);
	cout << Y_sum <<" "<< Y_mean <<" "<< Y_std;

	double cov_sum=0;
	double x=0;
	double y = 0;
	
	for (int i = 0; i < n; i++) {
		cov_sum =cov_sum + (X[i] - X_mean)*(Y[i] - Y_mean);
		x += (X[i] - X_mean)*(X[i] - X_mean);
		y += (Y[i] - Y_mean)*(Y[i] - Y_mean);
	}

	double tho_a = (cov_sum/(n-1)) / (sqrt(x/(n-1))*sqrt(y/(n-1)));
	
	cout << x<<endl;
	cout << "The simulation is "<< tho_a<<endl;

	//Question 2
	// construct normal distributed sequence X and Y, based on the standard normal sequence from question 1
	vector<double> Q2_X,Q2_Y;
	double Q2_tho = 0.6;
	cout << "Enter the correlation: ";
	//cin>>Q2_tho;
	for (int i = 0; i < n; i++) {
		Q2_X.push_back(standard_normal_sequence[i]);
		Q2_Y.push_back(standard_normal_sequence[i] * (Q2_tho)+standard_normal_sequence[i + n] * sqrt(1 - Q2_tho*Q2_tho));
	}
	double Q2_sum = 0;
	for (int i = 0; i < n; i++) {
		if (Q2_X[i] * Q2_X[i] * Q2_X[i] + sin(Q2_Y[i]) + Q2_X[i] * Q2_X[i] * Q2_Y[i]>0)
			Q2_sum += Q2_X[i] * Q2_X[i] * Q2_X[i] + sin(Q2_Y[i]) + Q2_X[i] * Q2_X[i] * Q2_Y[i];
	}
	cout << "The value E is " << Q2_sum / n << endl;

	//Question 3
	//use standard normal sequence from question 1 as base
	double t1 = 5; 
	double t2 = 0.5;
	double t3 = 3.2;
	double t4 = 6.5;
	//cin<<t;
	double Ea1;
	double Ea2;
	double Ea3; double Ea4;
	Ea1=Q3_a_simulation(standard_normal_sequence);
	Ea2 = Q3_a_simulation(standard_normal_sequence, t2); Ea3 = Q3_a_simulation(standard_normal_sequence, t3); 
	Ea4 = Q3_a_simulation(standard_normal_sequence, t4);
	cout << "Ea1 = "<<Ea1<<endl; 
	cout << "Ea2 = " << Ea2  << endl; 
	cout << "Ea3 = " << Ea3  << endl; 
	cout << "Ea4 = " << Ea4 << endl;

	//b)
	double Eb1;
	double Eb2;
	long double Eb3;long double Eb4;
	uniform_sequence.erase(uniform_sequence.begin(), uniform_sequence.end());//clean all the elements of uniform sequence from q1
	
	
	generate_uniform_distributed_sequence(uniform_sequence,n,seed);//Generate 20000 uniform distributed random numbers
	standard_normal_sequence.erase(standard_normal_sequence.begin(),standard_normal_sequence.end());//clean elements from q1
	generate_standard_normal_distributed_sequence(standard_normal_sequence, uniform_sequence);
	
	//solve Eb1,2,3,4 through antithetical method, but these functions are not monotone
	//Eb1 = Q3_b_simulation(standard_normal_sequence);
	//Eb2 = Q3_b_simulation(standard_normal_sequence, t2); Eb3 = Q3_b_simulation(standard_normal_sequence, t3);
	//Eb4 = Q3_b_simulation(standard_normal_sequence, t4);

	//control variate method
	// For E(W_5^2+sin(W_5), I choose Z = T*W^2+sin(U),U~[0,pi/2], E[Z]=T*E[w^2]+E[cos(u)]=T-1+pi/2
	//T=W^2*T+sin(W*sqrt(t))-T*W-cos(U)+0.5*pi
	vector<double> U_sequence_pi;
	generate_uniform_distributed_sequence(U_sequence_pi, n,seed);
	for (int i = 0; i < n; i++) {
		U_sequence_pi[i] *= (0.5*(atan(1) * 4));
	}

	Eb1 = 0;
	vector<double> test;//test sequence to evaluate the variance
		for (int i = 0; i < standard_normal_sequence.size(); i++) {
			Eb1 += standard_normal_sequence[i] * standard_normal_sequence[i] * 5 + sin(standard_normal_sequence[i] * sqrt(5))-(-4+5*standard_normal_sequence[i] * standard_normal_sequence[i] +cos(U_sequence_pi[i])-0.5*(atan(1) * 4));
			//store the generated number in the test sequence
			test.push_back(standard_normal_sequence[i] * standard_normal_sequence[i] * 5 + sin(standard_normal_sequence[i] * sqrt(5)) - (-4+5*standard_normal_sequence[i] * standard_normal_sequence[i] + cos(U_sequence_pi[i]) - 0.5*(atan(1) * 4)));
		}
	Eb1=Eb1 / standard_normal_sequence.size();
	double test_mean, test_sum, test_std;
	//calculate the variance of Eb1
	calculate_sum_mean_stdev(test, test_sum, test_mean, test_std);

	Eb2 = Eb3 = Eb4 = 0;//reset Eb234
	//for E[e^t/2*cos(W_t],choose Z=e^t/2*cos(U),U~[0,pi/2],E[Z]=e^t/2*(pi/2-1)
	vector<double> test2; vector<double> test3; vector<double> test4;
	t2 = t2 / 2; t3 = t3 / 2; t4 = t4 / 2;
	for (int i = 0; i < standard_normal_sequence.size(); i++) {
		Eb2 += exp(t2)*cos(standard_normal_sequence[i] * sqrt(t2)) - exp(t2)* (cos(U_sequence_pi[i]) - (0.5*(atan(1) * 4) - 1));
		Eb3 += exp(t3)*cos(standard_normal_sequence[i] * sqrt(t3)) - 3.2*exp(t3)*( cos(U_sequence_pi[i]) - (0.5*(atan(1) * 4) - 1));
		Eb4 += exp(t4)*cos(standard_normal_sequence[i] * sqrt(t4)) - 2*exp(t4)* (cos(U_sequence_pi[i]) - (0.5*(atan(1) * 4) - 1));
		test2.push_back(exp(t2)*cos(standard_normal_sequence[i] * sqrt(t2)) - exp(t2)* ( cos(U_sequence_pi[i]) - (0.5*(atan(1) * 4) - 1)));
		test3.push_back(exp(t3)*cos(standard_normal_sequence[i] * sqrt(t3)) - 3.2*exp(t3)* ( cos(U_sequence_pi[i]) - (0.5*(atan(1) * 4) - 1)));
		test4.push_back(exp(t4)*cos(standard_normal_sequence[i] * sqrt(t4)) - 2*exp(t4)*( cos(U_sequence_pi[i]) - (0.5*(atan(1) * 4) - 1)));
	}

	Eb2 = Eb2 / standard_normal_sequence.size(); Eb3 = Eb3 / standard_normal_sequence.size(); Eb4 = Eb4 / standard_normal_sequence.size();
	double test_mean2, test_sum2, test_std2; double test_mean3, test_sum3, test_std3; double test_mean4, test_sum4, test_std4;
	calculate_sum_mean_stdev(test2, test_sum2, test_mean2, test_std2); calculate_sum_mean_stdev(test3, test_sum3, test_mean3, test_std3); calculate_sum_mean_stdev(test4, test_sum4, test_mean4, test_std4);
	cout << "Eb1 = " << Eb1 << ", the variance = " << test_std*test_std << endl;
	cout << "Eb2 = " << Eb2 << ", the variance = " << test_std2*test_std2 << endl;
	cout << "Eb3 = " << Eb3 << ", the variance = " << test_std3*test_std3 << endl;
	cout << "Eb4 = " << Eb4 << ", the variance = " << test_std4*test_std4 << endl;

//Question 4
	double r = 0.04;
	double sigma = 0.2;
	double S0 = 88;
	//cin >> r >> sigma >> S0;
	double T = 5;
	double strike = 100;
	//cin>>T>>X;
	double c;
	c = european_call_option_price_mts(r,sigma,S0,T,strike);
	cout << "The empirical call option price trhough Monte Carlo simulation is " << c << endl;
	//b)
	c = european_call_option_price_vr(r, sigma, S0, T, strike);
	cout<<"The empirical call option price trhough Variance Reduction simulation is " << c << endl;
	
//Question 5
	//a)
	//vector<int> n_5(10);
	//std::iota(n_5.begin(), n_5.end(), 1);
	vector<double> St;
	St.push_back(88);
	r = 0.04; sigma = 0.18; S0 = 88; n = 1000;
	for (int i = 0; i < 10; i++) {
		double total = 0;
		for (int j = 0; j < n; j++) {
			double e = sigma*generate_normal_number() + (r - 0.5*sigma*sigma)*i;
			double St = S0*exp(e);
			total += St;
		}
		total = total / n;
		St.push_back(total);
	}

	//b)
	S0 = 88;
	double delta = 0.01;
	vector<vector<double>> path;
	for (int i=0; i < 6; i++) {
		vector<double> S_step;
		double temp = 88;
		S_step.push_back(88);
		for (int j = 0; j < n;j++) {
			vector<double> random;
			vector<double> uniform;
			generate_uniform_distributed_sequence(uniform, n,seed);
			generate_standard_normal_distributed_sequence(random, uniform);
				double e = (sigma*random[j]*sqrt(delta) + (r - 0.5*sigma*sigma)*delta);
				double S = temp*exp(e);
				temp = S;
				S_step.push_back(S);
			}
		path.push_back(S_step);
		}

	//c)
	std::ofstream myfile("question5_a.txt");
	if (myfile.is_open()) {
		for (int count = 0; count < St.size(); count++) {
			myfile << St[count];
			if (count != St.size() - 1)
				myfile << "" << endl;
		}
		myfile << '\n';
		myfile.close();
	}
	std::ofstream myfile2("question5_b.txt");
	if (myfile2.is_open()) {
		for (int i = 0; i < 6; i++) {
			std::vector<double> current = path[i];
			for (int count = 0; count < current.size(); count++) {
				myfile2 << current[count];
				if (count != current.size()-1)
					myfile2 << "	" << endl;
			}
			myfile2 << '\n';
		}
		myfile2.close();
	}
	//d)
	sigma = 0.35; S0 = 88;
	vector<double> St_d;
	St_d.push_back(88);
	for (int i = 0; i < 10; i++) {
		double total = 0;
		for (int j = 0; j < n; j++) {
			double e = sigma*generate_normal_number() + (r - 0.5*sigma*sigma)*i;
			double St = S0*exp(e);
			total += St;
		}
		total = total / n;
		St_d.push_back(total);
	}
	S0 = 88;
	vector<vector<double>> path_d;
	for (int i = 0; i < 6; i++) {
		vector<double> S_step_d;
		double temp = 88;
		S_step_d.push_back(88);
		for (int j = 0; j < n; j++) {
			vector<double> random;
			vector<double> uniform;
			generate_uniform_distributed_sequence(uniform, n,seed);
			generate_standard_normal_distributed_sequence(random, uniform);
			double e = (sigma*random[j] * sqrt(delta) + (r - 0.5*sigma*sigma)*delta);
			double S = temp*exp(e);
			temp = S;
			S_step_d.push_back(S);
		}
		path_d.push_back(S_step_d);
	}
	cout << path_d.size();
	std::ofstream myfile3("question5_d.txt");
	if (myfile3.is_open()) {
		for (int count = 0; count < St_d.size(); count++) {
			myfile3 << St_d[count];
			if (count != St_d.size()-1 )
				myfile3 << "" << endl;
		}
		myfile3 << '\n';
		myfile3.close();
	}
	std::ofstream myfile4("question5_d_2.txt");
	if (myfile4.is_open()) {
		for (int i = 0; i < 6; i++) {
			std::vector<double> current = path_d[i];
			for (int count = 0; count < current.size(); count++) {
				myfile4 << current[count];
				if (count != current.size() -1)
					myfile4 << "	" << endl;
			}
			myfile4 << '\n';
		}
		
		myfile4.close();
	}


	//Question 6
	//a)
	double x0 = 0;
	n = 1000;
	delta = 0.001;
	double integral = 0;
	for (int i = 0; i < n; i++){
	integral+=	4 * (sqrt(1 - x0*x0))*delta;
	x0 += delta;
	}
	cout << "pi equals " << integral<<endl;
	//b)
	vector<double> Q6b;
	integral = 0;
	generate_uniform_distributed_sequence(Q6b, n,seed);
	for (int i = 0; i < n; i++) {
		integral += 4 * (sqrt(1 - Q6b[i] * Q6b[i]));
	}
	cout << "Pi equals " << integral/n << " through Monte Carlo Simulation." << endl;
	//c)
	vector<double> Q6c;
	generate_uniform_distributed_sequence(Q6c, n,seed);
	integral = 0;
	for (int i = 0; i < n; i++) {
		integral += 4 * sqrt(1 - Q6c[i] * Q6c[i]) / ((1 - 0.74*Q6c[i] * Q6c[i]) / (1 - 0.74 / 3));
	}
	cout<<"Pi equals "<<integral / n<<" through importance sampling method"<<endl;
	std::system("pause");
	return 0;
}