#include "Generater.h"
#define _USE_MATH_DEFINES
#include <math.h> 
#include <ctime>
#include <cstdlib>
#include <random>
#include <iostream>>

Generater::Generater(int x)
{
	seed = x;
}


Generater::~Generater()
{
}

void Generater::setRandomSeed()
{
	srand((unsigned int)time(0));
	seed = (seed + rand() + 1);
}

double* Generater::generateUniforms(int size)
{
	setRandomSeed();
	long long a = (long long)pow(7, 5);
	long long m = (long long)pow(2, 31) - 1;

	long long* nums = new long long[size]; // an array to store all values calculated
	nums[0] = (a * seed) % m;
	for (int i = 1; i < size; i++)
	{
		nums[i] = (a * nums[i - 1]) % m;       
	}

	double* uniformNums = new double[size];
	for (int i = 0; i < size; i++)
	{
		// Divide by m to get the Uniforms btw 0 and 1
		uniformNums[i] = (double)nums[i] / (double)m;
	}
	// Clean the data
	delete[] nums;
	nums = NULL;
	return uniformNums;
}

double* Generater::generateNormals(int size)
{
	double* uniform_1 = generateUniforms(size*2);
	double* uniform_2 = generateUniforms(size*2);
	double* normal_PM = new double[size];
	int counter = 0;
	int i = 0;
	while (counter < size)
	{
		double U1 = uniform_1[i];
		double U2 = uniform_2[i];
		double V1 = 2 * U1 - 1;
		double V2 = 2 * U2 - 1;
		double W = V1 * V1 + V2 * V2;
		if (W < 1)
		{
			double Z1 = V1 * sqrt((-2.0 * log(W)) / W);
			double Z2 = V2 * sqrt((-2.0 * log(W)) / W);
			normal_PM[counter] = Z1;
			normal_PM[counter + 1] = Z2;
			counter = counter + 2;
		}
		i++;
	}

	delete[]uniform_1;
	uniform_1 = NULL;
	delete[]uniform_2;
	uniform_2 = NULL;
	return normal_PM;
}


double Halton_Seq(int index, int base)
{
	double f = 1, r = 0;
	while (index > 0) {
		f = f / base;
		r = r + f * (index% base);
		index = index / base;
	}
	return r;
}

double** generateUniformWithHalton(int size, int base1, int base2)
{
	double** Halton2D = new double*[2];
	Halton2D[0] = new double[size];
	Halton2D[1] = new double[size];
	for (int i = 1; i <= size; i++)
	{
		Halton2D[0][i - 1] = Halton_Seq(i, base1);
		Halton2D[1][i - 1] = Halton_Seq(i, base2);
	}
	return Halton2D;
}


double* Generater::generateBoxMullerNormals(int size, int base1, int base2)
{
	double* normal_BoxMuller = new double[size];
	double** uniform = generateUniformWithHalton(size, base1, base2);
	for (int i = 0; i < size; i++)
	{
		double U1 = uniform[0][i];
		double U2 = uniform[1][i];
		double Z1 = sqrt(-2 * log(U1)) * cos(2 * M_PI * U2);
		double Z2 = sqrt(-2 * log(U1)) * sin(2 * M_PI * U2);
		normal_BoxMuller[i] = Z1;
	}
	for (int i = 0; i < 2; i++) {
		delete[] uniform[i];
	}
	delete[] uniform;
	return normal_BoxMuller;
}


double Generater::calculateMean(double nums[], int total)
{
	// Function to calculate the mean of an array, given the size
	double mean = 0;
	for (int i = 0; i < total; i++)
	{
		mean = mean + nums[i];
	}
	mean = mean / total;
	return mean;
}

double Generater::calculateSD(double nums[], int total, double mean)
{
	// Function to calculate the sd of an array, given the size and mean
	double sd = 0;
	for (int i = 0; i < total; i++)
	{
		sd = sd + (nums[i] - mean) * (nums[i] - mean);
	}
	sd = sqrt(sd / total);
	return sd;
}

double Generater::cov(double x[], double y[], int total)
{
	double c = 0;
	double x_mean = calculateMean(x, total);
	double y_mean = calculateMean(y, total);
	for (int i = 0; i < total; i++)
	{
		c = c + (x[i] - x_mean) * (y[i] - y_mean);
	}
	c = c / double(total);
	return c;
}


double* Generater::generateDefaultNormal(int size)
{
	//setRandomSeed();
	double* uniformNums = new double[size];
	// random device class instance, source of 'true' randomness for initializing random seed
	std::random_device rd;

	// Mersenne twister PRNG, initialized with seed from previous random device instance
	std::mt19937 gen(rd());


	for (int i = 0; i < size; i++)
	{
		std::normal_distribution<double> d(0, 1.0);
		uniformNums[i] = d(gen);
	}
	return uniformNums;
}


double* Generater::generateDefaultUniform(int size)
{
	setRandomSeed();
	double* uniformNums = new double[size];
	// random device class instance, source of 'true' randomness for initializing random seed
	std::random_device rd;

	// Mersenne twister PRNG, initialized with seed from previous random device instance
	std::mt19937 gen(rd());


	for (int i = 0; i < size; i++)
	{
		std::uniform_real_distribution <double> d(0, 1.0);
		uniformNums[i] = d(gen);
	}
	return uniformNums;
}


double* Generater::generateExponential(int size,double lambda)
{
	double* uniforms = generateDefaultUniform(size);
	double* ExpoNums = new double[size];
	for (int i = 0; i < size; i++)
	{
		ExpoNums[i] = -1 / lambda * log(1 - uniforms[i]);
	}
	delete[] uniforms;
	return ExpoNums;
}