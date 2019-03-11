#include<iostream>
#include<vector>
#include<math.h>
#include<time.h>
#include<stdlib.h>
#include<array>
#include<numeric>
#include<cmath>
#include<fstream>

using namespace std;

int generate_Bernoulli_number(double x, double p);//generate Bernoulli number for given number and probability
void generate_uniform_distributed_sequence(vector<double> &RV, int N);//generate a sequence of uniform distributed number given the size, 
//the seed is generated by built-in function
void calculate_mean_stdev(vector<double> &RV1);//calculate the mean and standard deviation for given sequence
int main() {

/*Question 1*/
//a)
	int N = 10000; // times of simulation
	// The LGM seed
	int m = pow(2,31)-1;
	int a = pow(7,5);
	int b = 0;
	vector<float> RV;//create an array to store the random numbers
	int x_0 = 1; // initialize a seed for the random numbers
	unsigned long int x=(a*x_0+b)%m;//random number generated
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
	//store the sequence in a txt file
	std::ofstream myfile("question1_a.txt");
	if (myfile.is_open()) {
		for (int count = 0; count < RV.size(); count++) {
			myfile << RV[count];
			if (count != RV.size() - 1)
				myfile << "" << endl;
		}
		myfile << '\n';
		myfile.close();
	}
	//calculate the mean and the standard deviation
	double mean=accumulate(begin(RV),end(RV),0.0)/RV.size();
	double sum = 0;
	for (int i = 0; i < N; i++) {
		sum += (RV[i] - mean)*(RV[i] - mean);
	}
	double stdev = sqrt(sum / RV.size());
	cout << "The empirical mean and standard deviation are "<<mean<<" "<< stdev << endl;

// use built-in function to generate random variables
	vector<float> RV1;//create a vector to store the generated numbers
	// generate random numbers through rand()
	for (int i = 0; i < N; i++) {
		int random_number = 0 + rand() % 1000;
		double rn = random_number*0.001;
		RV1.push_back(rn);
	}
	//calculate the mean and standard deviation
	double mean1 = accumulate(begin(RV1), end(RV1), 0.0) / RV1.size();
	double sum1 = 0;
	for (int i = 0; i < N; i++) {
		sum1 += (RV1[i] - mean1)*(RV1[i] - mean1);
	}
	double stdev1 = sqrt(sum1 / RV1.size());
	cout << "The built-in function mean and standard deviation are " << mean1 << " " << stdev1 << endl;
	
	//Question 2
	// a)
	vector<double> Bou_dist;
	int X;
	for (int i = 0; i < N; i++) {
		if (0<=RV[i] && RV[i] < 0.3) {
			X = -1;
		}
		if (0.3 <= RV[i] && RV[i] < 0.65) {
			X = 0;
		}
		if (0.65 <= RV[i] && RV[i]< 0.85) {
			 X = 1;
		 }
		if(0.85 <= RV[i] && RV[i] <= 1) {
			X = 2;
		}
		Bou_dist.push_back(X);
	}

// b)
	std::ofstream myfile2("question2_b.txt");
		if (myfile2.is_open()) {
			for (int count = 0; count < Bou_dist.size(); count++) {
				myfile2 << Bou_dist[count];
				if (count != Bou_dist.size() - 1)
					myfile2 << "" << endl;
			}
			myfile2 << '\n';
			myfile2.close();
		}
	
/*Question 3*/
	// a)

	//use function to generate a uniform distribution
	vector<double> Binomial_sequence;
	int j = 0;//create 1000 binomial distributed number
	while(j<1000){
		//generate a sequence of uniform distribution
		vector<double> uniform_sequence;
		int X = 0;
		generate_uniform_distributed_sequence(uniform_sequence, 44);
			for (int i = 0; i < 44; i++) {
				X = X + generate_Bernoulli_number(uniform_sequence[i], 0.64);
			}
			j++;
			Binomial_sequence.push_back(X);
	}

	// b)
	// count how many times the P is greater than 40
	int count_p=0;
	for (int i = 0; i < 1000; i++) {
		if (Binomial_sequence[i] >= 40) {
			count_p++;
		}
	}
	double p = count_p*0.001;
	cout << "The probability that X is at least 40: P(X>=40) is " << p << endl;
	std::ofstream myfile3("question3_b.txt");
	if (myfile3.is_open()) {
		for (int count = 0; count < Binomial_sequence.size(); count++) {
			myfile3 << Binomial_sequence[count];
			if (count != Binomial_sequence.size() - 1)
				myfile3 << "" << endl;
		}
		myfile3 << '\n';
		myfile3.close();
	}
	/*Question 4*/
	//a)
	vector<double> exp_distributed_sequence;
	for (int i = 0; i < RV.size(); i++) {
		double Y = -1/1.5*log(1 - RV[i]);
		exp_distributed_sequence.push_back(Y);
	}

	// b)
	//compute P(X>=1)
	int p1 = 0; int p2 = 0;//define how many times have the 2 events happened in 10000 tiems.
	for (int i = 0; i < exp_distributed_sequence.size(); i++) {
		if (exp_distributed_sequence[i] >= 1) {
			p1++;
		}
		if (exp_distributed_sequence[i] >= 4) {
			p2++;
		}
	}

	cout << "P(X>=1)= " << p1*0.0001 << endl;
	cout << "P(X>=4)= " << p2*0.0001  << endl;
	// c)
	//calculate the mean and the standard deviation
	double exp_mean = accumulate(begin(exp_distributed_sequence), end(exp_distributed_sequence), 0.0) / exp_distributed_sequence.size();
	double exp_sum = 0;
	for (int i = 0; i < N; i++) {
		exp_sum += (exp_distributed_sequence[i] - exp_mean)*(exp_distributed_sequence[i] - exp_mean);
	}
	double exp_stdev = sqrt(exp_sum / exp_distributed_sequence.size());
	cout << "The empirical mean and standard deviation of Exponential distribution are " << exp_mean << " " << exp_stdev << endl;
	
	// cout the numbers to a txt file and draw the histogram in R
	std::ofstream myfile4("question4_c.txt");
	if (myfile4.is_open()) {
		for (int count = 0; count < exp_distributed_sequence.size(); count++) {
			myfile4 << exp_distributed_sequence[count];
			if (count != exp_distributed_sequence.size() - 1)
				myfile4 << "" << endl;
		}
		myfile4 << '\n';
		myfile4.close();
	}
/*Question 5*/
	// a)
	vector<double> uniform_sequnce_5;
	generate_uniform_distributed_sequence(uniform_sequnce_5, 5000);
	//b)
	vector<double> normal_distributed_sequence;
		for (int i = 0; i < 2500; i++) {
				double Z_1 = sqrt(-2 * log(uniform_sequnce_5[i]))*cos(2 * (atan(1) * 4)*uniform_sequnce_5[i + 2500]);
				double Z_2 = sqrt(-2 * log(uniform_sequnce_5[i]))*sin(2 * (atan(1) * 4)*uniform_sequnce_5[i + 2500]);
				normal_distributed_sequence.push_back(Z_1);
				normal_distributed_sequence.push_back(Z_2);
		}	
	cout << normal_distributed_sequence.size()<<" normal distributed numbers have been generated by Box-Muller Method"<<endl;

	//c) calculate the mean and standard deviation
	cout << "For Box-Muller Method, ";
	calculate_mean_stdev(normal_distributed_sequence);

	//d)
	
	vector<double> PMM_normal_distributed_sequence;
		for (int i = 0; i < 2500; i++) {
				double v1 = 2 * uniform_sequnce_5[i] - 1, v2 = 2 * uniform_sequnce_5[i + 2500] - 1;
				double w = pow(v1, 2) + pow(v2, 2);
				if (w <= 1) {
					double z1 = v1*sqrt((-2 * log(w) / w));
					double z2 = v2*sqrt((-2 * log(w) / w));
					PMM_normal_distributed_sequence.push_back(z1);
					PMM_normal_distributed_sequence.push_back(z2);
				}
		}
	
	cout << PMM_normal_distributed_sequence.size()<<" normal distributed numbers have been generated by Polar-Marsaglia Method."<<endl;
	
	//e) compute the empirical mean and standard deviation
	cout << "For Polar-Marsaglia Method, ";
	calculate_mean_stdev(PMM_normal_distributed_sequence);

	//f)
	clock_t start_time, end_time;
	start_time = clock();
	vector<double> normal_distributed_sequence_new;
	while (normal_distributed_sequence_new.size() < 20000) {
		int k = 1;
		for (int i = 0; i < uniform_sequnce_5.size() - k; i++) {
			if (normal_distributed_sequence_new.size() < 20000) {
				double Z_1 = sqrt(-2 * log(uniform_sequnce_5[i]))*cos(2 * (atan(1) * 4)*uniform_sequnce_5[i + k]);
				double Z_2 = sqrt(-2 * log(uniform_sequnce_5[i]))*sin(2 * (atan(1) * 4)*uniform_sequnce_5[i + k]);
				normal_distributed_sequence_new.push_back(Z_1);
				normal_distributed_sequence_new.push_back(Z_2);
			}
			else
				break;
		}
	}

	end_time = clock();
	cout << normal_distributed_sequence_new.size() << " normal distributed numbers have been generated by Box-Muller Method" << endl;
	double BM_run_time = -start_time + end_time;
	cout << "The run time of BM method is " << BM_run_time << endl;
	start_time = clock();

	vector<double> PMM_normal_distributed_sequence_new;
	while (PMM_normal_distributed_sequence_new.size() < 20000) {     //check if the sequence is smaller than 20000, if not add 1 to k to 
		int k = 1;
		for (int i = 0; i < uniform_sequnce_5.size() - k; i++) {
			if (PMM_normal_distributed_sequence_new.size() < 20000) {
				double v1 = 2 * uniform_sequnce_5[i] - 1, v2 = 2 * uniform_sequnce_5[i + k] - 1;
				double w = pow(v1, 2) + pow(v2, 2);
				if (w <= 1) {
					double z1 = v1*sqrt((-2 * log(w) / w));
					double z2 = v2*sqrt((-2 * log(w) / w));
					PMM_normal_distributed_sequence_new.push_back(z1);
					PMM_normal_distributed_sequence_new.push_back(z2);
				}
			}
			else
				break;
		}
	}

	end_time = clock();
	cout << PMM_normal_distributed_sequence_new.size() << " normal distributed numbers have been generated by Polar-Marsaglia Method." << endl;
	double PMM_run_time = end_time - start_time;//compute the run time
	cout << "The run time of PM method is " << PMM_run_time << endl;
	system("pause");
	return 0;
}

void generate_uniform_distributed_sequence(vector<double> &RV, int N) {
	
	int m = pow(2, 31) - 1;
	int a = pow(7, 5);
	int b = 0;
	
	int x_0 = rand(); // initialize a seed for the random numbers
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

int generate_Bernoulli_number(double x,double p) {
	if (x <= p)
		return 1;
	else
		return 0;
}

//calculate the mean and standard deviation
void calculate_mean_stdev(vector<double> &RV1) {
	double mean1 = accumulate(begin(RV1), end(RV1), 0.0) / RV1.size();
	double sum1 = 0;
	for (int i = 0; i < RV1.size(); i++) {
		sum1 += (RV1[i] - mean1)*(RV1[i] - mean1);
	}
	double stdev1 = sqrt(sum1 / RV1.size());
	cout << "The empirical mean and standard deviation are " << mean1 << " " << stdev1 << endl;
}
