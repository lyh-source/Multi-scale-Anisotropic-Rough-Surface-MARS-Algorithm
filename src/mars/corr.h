#pragma once
#include<Eigen/core>
class Correlation {
public:
	Correlation() {};
	~Correlation() {};
	void setValue(const int& n_, const int& m_,const double& gamma_);

	int n, m;
	double gamma;
	Eigen::MatrixXd rhs;
};