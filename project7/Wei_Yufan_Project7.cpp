#include<iostream>
#include"Header.h"
#include"Header7.h"
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
#include<Eigen/Dense>

using namespace Eigen;

using namespace std;



double sigma = 0.2;
double r = 0.04;
double K = 10;
double T = 0.5;
double S0 = 10;
double Smin = 4; double Smax = 16;
double delta_t = 0.002;
double delta_x(int i) {
	switch (i) {
	case 1:
		return sigma*sqrt(delta_t);
		break;
	case 2:
		return sigma*sqrt(3 * delta_t);
		break;
	case 3:
		return sigma*sqrt(4 * delta_t);
		break;
	}
}

double* EFD(int select,int N) {
	double dx = delta_x(select);
	
	//double ds = (Smax - Smin) / N;
	double* Fi = new double[2*N + 1]{ 0 };
	int m = T / delta_t;
	double* FT = new double[2*N + 1]{ 0 }; 
	double* ST = new double[2 * N + 1]{ 0 }; ST[N] = S0; 
	FT[N] = max(K - ST[N], 0.0);
	for (int i = 1; i < N + 1; i++) {
		ST[N+i] = exp(log(ST[N+i-1]) - dx);
		ST[N - i] = exp(log(ST[N-i+1]) + dx);
		FT[N+i] = max(K - ST[N+i], 0.0);
		FT[N - i] = max(K - ST[N - i], 0.0);
	}
	double pu = delta_t*(sigma*sigma / (2 * dx*dx) + (r - sigma*sigma / 2) / (2 * dx));
	double pd = delta_t*(sigma*sigma / (2 * dx*dx) - (r - sigma*sigma / 2) / (2 * dx));
	double pm = 1 - delta_t*sigma*sigma / (dx*dx) - r*delta_t;
	double** A = new double*[2*N + 1];
	for (int i = 0; i < 2*N + 1; i++) {
		A[i] = new double[2*N + 1]{ 0 };
	}
	A[0][0] = pu; A[0][1] = pm; A[0][2] = pd;
	A[2*N][2*N - 2] = pu; A[2*N][2*N-1] = pm; A[2*N][2*N] = pd;
	for (int i = 1; i < 2*N; i++) {
		A[i][i-1] = pu; A[i][i] = pm; A[i][i + 1] = pd;
	}
	double* B = new double[2*N + 1]{ 0 };
	B[2*N] = -(ST[2*N]-ST[2*N-1]);

	for (int i = 0; i < 2*N+1; i++) {
		Fi[i] = FT[i];
	//	cout << B[i]<<" ";
	//	cout << Fi[i] << endl;
	}

	for (int i = m; i > 0; i--) {
		double* Fi_temp = new double[2*N + 1];
		Fi_temp = derive_F(Fi, A, B, 2*N + 1, m);
		for (int i = 0; i < 2*N + 1; i++) {
			Fi[i] = Fi_temp[i];
			//cout << Fi_temp[i]<<" ";
			//cout << Fi[i] << endl;
		}
		delete[]Fi_temp;
	}
	//cout << ST[N];
	for (int i = 0; i < 2 * N + 1; i++) {
		cout << ST[i] << " " << Fi[i] << endl;
	}
	return Fi;
}
void IFD(int select, int N,VectorXd& F) {
	double dx = delta_x(select);

	//double ds = (Smax - Smin) / N;
	int m = T / delta_t;
	
	VectorXd FT(2 * N + 1);
	VectorXd ST(2 * N + 1);
	ST(N) = S0;
	FT(N) = max(K - ST[N], 0.0);
	for (int i = 1; i < N + 1; i++) {
		ST[N + i] = exp(log(ST[N + i - 1]) - dx);
		ST[N - i] = exp(log(ST[N - i + 1]) + dx);
		FT[N + i] = max(K - ST[N + i], 0.0);
		FT[N - i] = max(K - ST[N - i], 0.0);
	}
	//cout << ST;
	//cout<<FT;
	double pu = -0.5*delta_t*(sigma*sigma/(dx*dx)+ (r - sigma*sigma / 2) /dx);
	double pd = -0.5*delta_t*(sigma*sigma / (dx*dx) - (r - sigma*sigma / 2) / dx);
	double pm = 1 + delta_t*sigma*sigma / (dx*dx) + r*delta_t;
	MatrixXd A(2 * N + 1,2*N+1);
	for (int i = 0; i < 2 * N + 1; i++) {
		for (int j = 0; j < 2 * N + 1; j++) {
			A(i, j) = 0;
		}
	}
	A(0,0)= 1; A(0,1) = -1;
	A(2 * N,2 * N - 1) = 1; A(2 * N,2 * N) = -1;
	for (int i = 1; i < 2 * N; i++) {
		A(i,i - 1) = pu; A(i,i) = pm; A(i,i + 1) = pd;
	}
	MatrixXd A_inv(2 * N + 1, 2 * N + 1);
	//cout << A;
	A_inv = A.inverse();
	//cout << A_inv;
	VectorXd B(2 * N + 1); B(0) = 0;
	B[2 * N] = -ST[2 * N] + ST[2 * N - 1];
	//cout << ST[2 * N] << " " << ST[2 * N - 1];
	for (int i = 1; i < 2 * N; i++) {
		B[i] = FT[i];
	}
	cout << B;
	VectorXd Fi(2 * N + 1);
	for (int i = m; i > 0; i--) {
		VectorXd Fi_temp(2 * N + 1);
		Fi_temp = A_inv*B;
		Fi(0) = Fi_temp(0); Fi(2 * N) = Fi_temp(2 * N);
		for (int i = 1; i < 2 * N; i++) {
			B[i] = Fi_temp[i];
			Fi(i) = Fi_temp(i);
			//cout << Fi_temp[i]<<" ";
			//cout << Fi[i] << endl;
		}
//		cout << A_inv*B;
//		cout << Fi_temp;
//		cout << Fi;
	}
	//cout << ST[N];
	for (int i = 0; i < 2 * N + 1; i++) {
		F(i) = Fi(i);
		cout << ST(i) << " " << Fi(i)<<endl;
	}
	
	
}
void CNFD(int select, int N, VectorXd& F) {
	double dx = delta_x(select);

	//double ds = (Smax - Smin) / N;
	int m = T / delta_t;

	VectorXd FT(2 * N + 1);
	VectorXd ST(2 * N + 1);
	ST(N) = S0;
	FT(N) = max(K - ST[N], 0.0);
	for (int i = 1; i < N + 1; i++) {
		ST[N + i] = exp(log(ST[N + i - 1]) - dx);
		ST[N - i] = exp(log(ST[N - i + 1]) + dx);
		FT[N + i] = max(K - ST[N + i], 0.0);
		FT[N - i] = max(K - ST[N - i], 0.0);
	}
	//cout << ST;
	//cout<<FT;
	double pu = -0.25*delta_t*(sigma*sigma / (dx*dx) + (r - sigma*sigma / 2) / dx);
	double pd = -0.25*delta_t*(sigma*sigma / (dx*dx) - (r - sigma*sigma / 2) / dx);
	double pm = 1 + delta_t*sigma*sigma / (2*dx*dx) +0.5* r*delta_t;
	MatrixXd A(2 * N + 1, 2 * N + 1);
	for (int i = 0; i < 2 * N + 1; i++) {
		for (int j = 0; j < 2 * N + 1; j++) {
			A(i, j) = 0;
		}
	}
	A(0, 0) = 1; A(0, 1) = -1;
	A(2 * N, 2 * N - 1) = 1; A(2 * N, 2 * N) = -1;
	for (int i = 1; i < 2 * N; i++) {
		A(i, i - 1) = pu; A(i, i) = pm; A(i, i + 1) = pd;
	}
	MatrixXd A_inv(2 * N + 1, 2 * N + 1);
	//cout << A;
	A_inv = A.inverse();
	//cout << A_inv;
	VectorXd B(2 * N + 1); B(0) = 0;
	B[2 * N] = -ST[2 * N] + ST[2 * N - 1];
	//cout << ST[2 * N] << " " << ST[2 * N - 1];
	for (int i = 1; i < 2 * N; i++) {
		B[i] = -pu*FT(i-1)-(pm-2)*FT(i)-pd*FT(i+1);
	}
	VectorXd Fi(2 * N + 1);
	for (int i = m; i > 0; i--) {
		VectorXd Fi_temp(2 * N + 1);
		Fi_temp = A_inv*B;
		Fi(0) = Fi_temp(0); Fi(2 * N) = Fi_temp(2 * N);
		for (int i = 1; i < 2 * N; i++) {
			B[i] = -pu*Fi_temp(i - 1) - (pm - 2)*Fi_temp(i) - pd*Fi_temp(i + 1);
			Fi(i) = Fi_temp(i);
			//cout << Fi_temp[i]<<" ";
			//cout << Fi[i] << endl;
		}
		//		cout << A_inv*B;
		//		cout << Fi_temp;
		//		cout << Fi;
	}
	//cout << ST[N];
	for (int i = 0; i < 2 * N + 1; i++) {
		F(i) = Fi(i);
		cout << ST(i) << " " << Fi(i) << endl;
	}


}

void EFD_american(string type,VectorXd& F) {
	double ds=0.25;
	int N = 17/ds;
	//double ds = (Smax - Smin) / N;
	VectorXd Fi(N+1);
	for (int i = N; i >= 0; i--) {
		Fi(i) = 0;
	}
	int m = T / delta_t;
	VectorXd FT(N + 1);
	for (int i = N; i >= 0; i--) {
		Fi(i) = 0;
	}
	VectorXd ST(N + 1);
	for (int i = N; i >= 0; i--) {
		ST(i) = 17-i*ds;
		if (type == "call")
			FT[i] = max(ST[i] - K, 0.0);
		else
			FT[i] = max(K - ST[i], 0.0);
	}

	MatrixXd A(N+1, N+1);
	for (int i = 0; i < N + 1; i++) {
		for (int j = 0; j < N + 1; j++) {
			A(i, j) = 0;
		}
	}
	A(0,0) = delta_t*(r*ST(N-1)/(2*ds*ds)+sigma*sigma*ST(N-1)/ds*ST(N-1)/ds*0.5);
	A(0,1) = 1-delta_t*(sigma*sigma*ST(N-1)/ds*ST(N-1)/ds+r);
	A(0,2) = delta_t*(-r*ST(N - 1) / (2*ds) + sigma*sigma*ST(N - 1)/ds*ST(N - 1)/ds*0.5);
	A(N,N - 2) = delta_t*(r / 2 + sigma*sigma*0.5);
	A(N,N - 1) = 1 - delta_t*(sigma*sigma + r);
	A(N,N) = delta_t*(-r / 2 + sigma*sigma*0.5);
	for (int i = 1; i < N; i++) {
		A(i,i - 1) = delta_t*(r*ST(i)/ds / 2 + sigma*sigma*ST(i)*ST(i)/ds/ds*0.5);
		A(i,i) = 1 - delta_t*(sigma*sigma*ST(i)*ST( i)/ds/ds + r);
		A(i,i + 1) = delta_t*(-r*ST(i) /ds/2 + sigma*sigma*ST(i)*ST(i)/ds/ds*0.5);
	}
	VectorXd B(N + 1);
	for (int i = 0; i < N + 1; i++) {
		B(i) = 0;
	}
	if (type == "put")
		B[0] = -(ST[N] - ST[N - 1]);
	else
		B[N] = ST[N] - ST[N - 1];
	for (int i = 0; i < N-1; i++) {
		Fi[i] = FT[i];
		//	cout << B[i]<<" ";
		//	cout << Fi[i] << endl;
	}

	for (int i = m; i > 0; i--) {
		VectorXd Fi_temp(N + 1);
		Fi_temp = A*Fi - B;
		for (int i = 0; i < N-1; i++) {
			if (Fi_temp[i] > FT[i])
				Fi[i] = Fi_temp[i];
			else
				Fi[i] = FT[i];
			//cout << Fi_temp[i]<<" ";
			//cout << Fi[i] << endl;
		}
	}
	//cout << ST[N];
	for (int i = 0; i < N + 1; i++) {
		F[i] = Fi[i];
	}
}
void IFD_american(string type, VectorXd& F) {
	double ds = 0.25;
	int N = 17 / ds;
	//double ds = (Smax - Smin) / N;
	VectorXd Fi(N + 1);
	for (int i = N; i >= 0; i--) {
		Fi(i) = 0;
	}
	int m = T / delta_t;
	VectorXd FT(N + 1);
	for (int i = N; i >= 0; i--) {
		Fi(i) = 0;
	}
	VectorXd ST(N + 1);
	for (int i = N; i >= 0; i--) {
		ST(i) = 17 - i*ds;
		if (type == "call")
			FT[i] = max(ST[i] - K, 0.0);
		else
			FT[i] = max(K - ST[i], 0.0);
	}
	MatrixXd A=MatrixXd::Zero(N+1,N+1);

	A(0, 0) = -1;
	A(0, 1) = 1;
	A(0, 2) = 0;
	A(N, N - 2) =0;
	A(N, N - 1) = 1;
	A(N, N) = -1;
	for (int i = 1; i < N; i++) {
		A(i, i - 1) = -delta_t*(r*ST(i) / ds / 2 + sigma*sigma*ST(i)*ST(i) / ds / ds*0.5);
		A(i, i) = 1 + delta_t*(sigma*sigma*ST(i)*ST(i) / ds / ds + r);
		A(i, i + 1) = -delta_t*(-r*ST(i) / ds / 2 + sigma*sigma*ST(i)*ST(i) / ds / ds*0.5);
	}
	MatrixXd A_inv(N + 1,N + 1);
//	cout << A;
	A_inv = A.inverse();
//	cout << A_inv;
	VectorXd B(N + 1); 
	if (type == "put") {
		B(0) = 0;
		B[N] = -ST[N] + ST[N - 1];
	}
	else {
		B(0) = ST(N) - ST(N - 1);
		B(N) = 0;
	}//cout << ST[2 * N] << " " << ST[2 * N - 1];
	for (int i = 1; i < N+1; i++) {
		B[i] = FT[i];
	}

	for (int i = m; i > 0; i--) {
		VectorXd Fi_temp(N + 1);
		Fi_temp = A_inv*B;
		//cout << Fi_temp;
			Fi(0) = Fi_temp(0)>FT(0)? Fi_temp(0):FT(0);
			Fi(N) = Fi_temp(N)>FT(0) ? Fi_temp(N) : FT(N);
		
		for (int i = 1; i < N+1; i++) {
			B[i] = Fi_temp[i]>FT(i)?Fi_temp(i):FT(i);
			Fi(i) = Fi_temp(i)>FT(i)?Fi_temp(i) : FT(i);
			//cout << Fi_temp[i]<<" ";
			//cout << Fi << endl;
		}
		//		cout << A_inv*B;
		//		cout << Fi_temp;
		//		cout << Fi;
	}
	//cout << ST[N];
	for (int i = 0; i <N + 1; i++) {
		F(i) = Fi(i);
	}
}
void CNFD_american(string type, VectorXd& F) {
	double ds = 0.25;
	int N = 17 / ds;
	//double ds = (Smax - Smin) / N;
	VectorXd Fi(N + 1);
	for (int i = N; i >= 0; i--) {
		Fi(i) = 0;
	}
	int m = T / delta_t;
	VectorXd FT(N + 1);
	for (int i = N; i >= 0; i--) {
		Fi(i) = 0;
	}
	VectorXd ST(N + 1);
	for (int i = N; i >= 0; i--) {
		ST(i) = 17 - i*ds;
		if (type == "call")
			FT[i] = max(ST[i] - K, 0.0);
		else
			FT[i] = max(K - ST[i], 0.0);
	}
	MatrixXd A = MatrixXd::Zero(N + 1, N + 1);

	A(0, 0) = -1;
	A(0, 1) = 1;
	A(0, 2) = 0;
	A(N, N - 2) = 0;
	A(N, N - 1) = 1;
	A(N, N) = -1;
	for (int i = 1; i < N; i++) {
		A(i, i - 1) = -0.5*delta_t*(r*ST(i) / ds / 2 + sigma*sigma*ST(i)*ST(i) / ds / ds*0.5);
		A(i, i) = 1 + 0.5*delta_t*(sigma*sigma*ST(i)*ST(i) / ds / ds + r);
		A(i, i + 1) = -0.5*delta_t*(-r*ST(i) / ds / 2 + sigma*sigma*ST(i)*ST(i) / ds / ds*0.5);
	}
	MatrixXd A_inv(N + 1, N + 1);
	//	cout << A;
	A_inv = A.inverse();
	//	cout << A_inv;
	VectorXd B(N + 1);
	if (type == "put") {
		B(0) = 0;
		B[N] = -ST[N] + ST[N - 1];
	}
	else {
		B(0) = ST(N) - ST(N - 1);
		B(N) = 0;
	}//cout << ST[2 * N] << " " << ST[2 * N - 1];
	for (int i = 1; i < N + 1; i++) {
		B[i] = FT[i];
	}

	for (int i = m; i > 0; i--) {
		VectorXd Fi_temp(N + 1);
		Fi_temp = A_inv*B;
		//cout << Fi_temp;
		Fi(0) = Fi_temp(0)>FT(0) ? Fi_temp(0) : FT(0);
		Fi(N) = Fi_temp(N)>FT(0) ? Fi_temp(N) : FT(N);

		for (int i = 1; i < N; i++) {
			double S = ST[i];
			double Pu = -0.25 * delta_t * ((sigma*sigma * S * S) / (ds * ds) + (r * S) / ds);
			double Pm = 1 + (delta_t * (sigma * sigma * S * S)) / (2 * ds * ds) + r * delta_t / 2.0;
			double Pd = -0.25 * delta_t * ((sigma*sigma * S * S) / (ds * ds) - (r * S) / ds);
			B(i) = -Pu * Fi(i - 1) - (Pm - 2) * Fi(i) - Pd * Fi(i + 1);
			Fi(i) = Fi_temp(i)>FT(i) ? Fi_temp(i) : FT(i);
			//cout << Fi_temp[i]<<" ";
			//cout << Fi << endl;
		}
		//		cout << A_inv*B;
		//		cout << Fi_temp;
		//		cout << Fi;
	}
	//cout << ST[N];
	for (int i = 0; i <N + 1; i++) {
		F(i) = Fi(i);
	}
}

int main() {
	int N = 100;
	double* F0 = new double[2*N + 1];
	//F0=EFD(1,N);
	//cout << F0[0]<<endl;
	VectorXd prices(2*N+1);
	IFD(1, N,prices);
	//cout << "q2"<<prices;
	VectorXd prices_1c(2 * N + 1);
	CNFD(1, N, prices_1c);
	//cout << "q3" << prices_1c;

	/*Question 2*/
	
	
	VectorXd F0_american(69);
	EFD_american("call", F0_american);
	cout << F0_american;
	cout << F0_american[28];

	IFD_american("call", F0_american);
	cout << F0_american;
	cout << F0_american[28];

	
	CNFD_american("call",F0_american);
	cout<<F0_american;
	cout << F0_american[28];

	EFD_american("put", F0_american);
	cout << F0_american;
	cout << F0_american[28];

	IFD_american("put", F0_american);
	cout << F0_american;
	cout << F0_american[28];


	CNFD_american("put", F0_american);
	cout << F0_american;
	cout << F0_american[28];

	system("pause");
	return 0;
	
}