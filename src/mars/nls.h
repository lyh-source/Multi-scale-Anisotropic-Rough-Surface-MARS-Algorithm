#pragma once
#include<Eigen/core>
class NonLinearSystem {
public:
	NonLinearSystem() { ff = NULL; xx = NULL; };
	~NonLinearSystem() {};
	void setValue(const double& epsilon_, const int& max_iter_);
	Eigen::MatrixXd nls_init(const Eigen::MatrixXd& Rpq);
	Eigen::MatrixXd nls_rhs(const Eigen::MatrixXd& alpha, const Eigen::MatrixXd& rhs);
	Eigen::MatrixXd nls_solve(const Eigen::MatrixXd& guess, int max_iter, double tol, const Eigen::MatrixXd& rhs);
	//void nls_rhs2(int n, double x[], double fx[]);

	double epsilon;
	int max_iter;
	int success;
	Eigen::MatrixXd init;
	Eigen::MatrixXd coeff;
	Eigen::MatrixXd solve;


	double* ff;
	double* xx;
	//int rhs_n, rhs_m;
};