#include"mars.h"
#include"fsolve.h"
#include<iostream>
#include<Eigen/core>
void mars_algorithm::MARS::mars_input(const double& Lx_, const double& Ly_, const int& Nx_, const int& Ny_, const int& N_,
	const int& n_, const int& m_, const double& gamma_,
	const double& epsilon_, const int& max_iter_,
	const int& option_, const int& cutoff_) {


	std::cout << "input\n";

	hmap.setValue(Lx_, Ly_, Nx_, Ny_, N_);
	corr.setValue(n_, m_,gamma_);
	filt.setValue(option_, cutoff_);

	nls.setValue(epsilon_, max_iter_);
	nls.success = 0;

}


void meshgrid(const Eigen::VectorXd& x, const Eigen::VectorXd& y, Eigen::MatrixXd& xx, Eigen::MatrixXd& yy)
{
	Eigen::MatrixXd m_x = x.matrix();
	Eigen::MatrixXd m_y = y.matrix();
	int Nx = m_x.rows();
	int Ny = m_y.rows();
	m_x.transposeInPlace();
	xx = m_x.replicate(Ny, 1);
	yy = m_y.replicate(1, Nx);
}

void mars_algorithm::MARS::mars_step_one() {
	std::cout << "Step 1\n";
	Eigen::VectorXd xpg = Eigen::VectorXd::LinSpaced(hmap.Nx, 0, hmap.Lx);
	Eigen::VectorXd ypg = Eigen::VectorXd::LinSpaced(hmap.Ny, 0, hmap.Ly);
	//std::cout << xpg.size()<<" "<<xpg.rows()<<" "<<xpg.matrix() << "\n";

	meshgrid(ypg, xpg, hmap.Y, hmap.X);
	hmap.raw.clear();
	hmap.filtered.clear();

	for (int i = 0; i < hmap.N; i++) {
		hmap.raw.push_back(Eigen::MatrixXd::Zero(hmap.Nx, hmap.Ny));
		hmap.filtered.push_back(Eigen::MatrixXd::Zero(hmap.Nx, hmap.Ny));
	}
	
	//std::cout << hmap.raw.cols() << " \n";

}


void mars_algorithm::MARS::mars_step_two() {
	std::cout << "Step 2\n";
	int n0p1 = corr.n, m0p1 = corr.m;



} 

void mars_algorithm::MARS::mars_step_three() {
	std::cout << "Step 3\n";
}

void mars_algorithm::MARS::mars_step_four() {
	std::cout << "Step 4\n";
}
void mars_algorithm::MARS::mars_step_five(int i) {
	std::cout << "Step 5\n";
}

void mars_algorithm::MARS::solve() {
	mars_step_one();

	nls.success = 0;
	while (nls.success == 0) {
		mars_step_two();
		mars_step_three();

		if (nls.success == 0) {
			corr.gamma = corr.gamma * 0.5;
		}
		else {
			std::cout << "nls success.\n";
		}

		nls.success = 1;////////////////////////////////////////////////////////////////////////////////
	}
	

	for (int i = 0; i < hmap.N; i++) {
		mars_step_four();

		if (filt.option == 1) {
			mars_step_five(i);
		}
		else {
			hmap.filtered[i] = hmap.raw[i];
		}
	}
	
	
}