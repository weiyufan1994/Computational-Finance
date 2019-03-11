#pragma once

class Generater
{
public:
	Generater(int x);
	~Generater();
	
	double* generateUniforms(int size);
	double* generateExponential(int size, double lambda);
	double* generateNormals(int size);
	double* generateBoxMullerNormals(int size, int base1, int base2);
	double* generateDefaultUniform(int size);

	double calculateMean(double nums[], int total);
	double calculateSD(double nums[], int total, double mean);

	double cov(double x[], double y[], int total);
	int seed = 0;

	double* generateDefaultNormal(int size);

private:
	void setRandomSeed();
};

