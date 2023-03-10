#pragma once
#include<Eigen/core>
#include<vector>
class HeightMap {
public:
	HeightMap() {};
	~HeightMap() {};
	void setValue(const double& Lx_, const double& Ly_, const int& Nx_, const int& Ny_, const int& N_);


	double Lx,Ly;
	int Nx, Ny,N;
	Eigen::MatrixXd X, Y;
	std::vector<Eigen::MatrixXd > raw, filtered;
};
