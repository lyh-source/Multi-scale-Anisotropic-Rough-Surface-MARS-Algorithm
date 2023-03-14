#include"mars.h"
#include"fsolve.h"
#include<iostream>
#include<Eigen/core>
#include <random>
#include <unsupported/Eigen/FFT>
#include<opencv2/core/eigen.hpp>

#include<fstream>
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
	//corr.gamma = 0.0064;///////////////////////////////////////////////////////////////////////////////////////////
	int n = floor((double)corr.n * (log(corr.gamma) / -2.3));
	int m = floor((double)corr.m * (log(corr.gamma) / -2.3));
	std::cout << n<<" " << m << "\n";
	Eigen::MatrixXd tmp = Eigen::MatrixXd::Zero(n, m);

	for (int p = 0; p < n; p++) {
		for (int q = 0; q < m; q++) {
			double arg = -2.3 * sqrt(pow(((double)p / (double)n0p1), 2.0) + pow((double)q / (double)m0p1, 2.0));
			//std::cout << arg << "\n";
			tmp(p, q) = exp(arg);
		}
	}


	//std::cout << tmp << "\n";
	corr.rhs = tmp;


} 

void mars_algorithm::MARS::mars_step_three() {
	std::cout << "Step 3\n";
	
	nls.init = nls.nls_init(corr.rhs);
	//std::cout << nls.init;
	nls.nls_solve(nls.init, nls.max_iter, nls.epsilon, corr.rhs);



	//Evaluate RHS of the non - linear system using coeffcients
	nls.solve = nls.nls_rhs(nls.coeff, corr.rhs);
	//nls.nls_rhs(nls)


}




double getStd(const Eigen::MatrixXd& a) {

	int n = a.rows();
	int m = a.cols();
	Eigen::VectorXd etaR = Eigen::VectorXd::Zero(n * m);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			etaR(i * m + j) = a(i, j);
		}
	}


	double mean = etaR.mean();
	for (int i = 0; i < etaR.rows(); i++)
	{
		etaR(i) = etaR(i) - mean;
	}

	double std = sqrt((etaR.transpose() * etaR)[0] / (double)(etaR.size()));
	return std;
	//std::cout << "std£º " << std << " ]]]]]]]]]]]]]]]]]]]\n";
}

void mars_algorithm::MARS::mars_step_four(int idx) {
	std::cout << "Step 4\n";

	int N = hmap.Nx;
	int M = hmap.Ny;

	int n = nls.coeff.rows();
	int m = nls.coeff.cols();

	Eigen::MatrixXd eta = Eigen::MatrixXd::Zero(N, M);


	//Gaussian random number matrix
	std::random_device rd{};
	std::mt19937 gen{ rd() };
	std::normal_distribution<double> distribution(0.0, 1);
	/*for (int n = 0; n != 10000; ++n)
		std::cout << distribution(gen) << std::endl;*/

	/*eta.setRandom(distribution(gen));*/

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			eta(i, j) = distribution(gen);
		}
	}
	//eta = eta - eta.mean();


	double mean = eta.mean();
	Eigen::MatrixXd mm = Eigen::MatrixXd::Zero(N, M);
	mm.setConstant(mean);
	
	eta = eta - mm;
	
	double std = getStd(eta);
	double fac = 1.0 / std;
	eta = fac * eta;


	//std::cout << "std£º " << eta.row(1) << " ]]]]]]]]]]]]]]]]]]]\n";
	Eigen::MatrixXd tmp = Eigen::MatrixXd::Zero(N, M);

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			for (int k = 0; k < n; k++) {
				for (int l = 0; l < m; l++) {
					tmp(i, j) = tmp(i, j) + nls.coeff(k, l) * eta((i + k+1) % N, (j + l+1) % M);
				}
			}
		}
	}


	mean = tmp.mean();
	mm.setConstant(mean);

	tmp = tmp - mm;
	std = getStd(tmp);
	fac = 1.0 / std;
	tmp = fac * tmp;

	std::cout << "std£º " << getStd(tmp) <<" "<<tmp.rows()<<" " << tmp.cols() << " ]]]]]]]]]]]]]]]]]]]\n";
	hmap.raw[idx] = tmp;




}
void mars_algorithm::MARS::mars_step_five(int idx) {
	std::cout << "Step 5\n";
	int mmax = filt.cutoff;
	double nmax = ((double)hmap.Ly) / ((double)hmap.Lx) * (double)filt.cutoff;

	int set_ellipse = 1;
	cv::Mat fIn,fOut,fTmp;

	cv::eigen2cv(hmap.raw[idx], fIn);
	//cv::eigen2cv(nls.coeff, fIn);
	//std::cout << fIn.row(0) << " kkkkkkkkkkkk\n";
	//std::cout << hmap.raw[idx].rows()<<" "<< hmap.raw[idx].cols() << "  size\n";
	//std::cout << fIn.rows<< "  size\n";



	cv::transpose(fIn, fTmp);
	fIn = fTmp;
	dft(fIn, fOut,  cv::DFT_COMPLEX_OUTPUT);
	cv::transpose(fOut, fTmp);
	fOut = fTmp;

	std::cout << fOut.size << "  sizeo\n";

	int M = fOut.rows;
	int N = fOut.cols;

	Eigen::MatrixXd mask = Eigen::MatrixXd::Zero(M, N);

	for (int i = 0; i < mmax + 1; i++) {
		for (int j = 0; j < nmax + 1; j++) {
			if (i == 0 && j == 0)
				continue;

			if (set_ellipse == 1) {
				double kel = pow(((double)i - 1.0 + 1.0) / (double)mmax, 2) + pow(((double)j - 1.0 + 1.0) / (double)nmax, 2);

				if (kel > 1)
					continue;
			}
			/*cout << "*********************************************\n";

			std::cout << M - i << " " << N - j << "\n";
			std::cout << i << " " <<  j << "\n";
			std::cout << M << " " <<  N << "\n";

			cout << "*********************************************\n";*/
			mask(i, j) = 1;
			if (i > 0 && j > 0)
				mask(M - i, N - j) = 1;
			else if (i > 0)
				mask(M - i, j) = 1;
			else if (j > 0)
				mask(i, N - j) = 1;

		}
	}

	for (int i = 0; i > -mmax; i--) {
		//cout << i << " " << -mmax << "\n";

		for (int j = 1; j < nmax + 1; j++) {
			if (set_ellipse == 1) {
				double kel = pow(((double)i - 1.0) / (double)mmax, 2) + pow(((double)j - 1.0 + 1.0) / (double)nmax, 2);
				if (kel > 1)
					continue;
			}



			mask(M + i - 1, j) = 1;
			mask(-i + 1, N - j) = 1;


		}
	}

	//cout << fOut << "\n";

	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			fOut.at<std::complex<double>>(i, j) *= (double)mask(i, j);
		}
	}
	
	std::cout << mask.rows()<<" "<< mask.cols() << " \nllllllllllllll";

	cv::dft(fOut, fIn, cv::DFT_SCALE|cv::DFT_INVERSE | cv::DFT_COMPLEX_INPUT | cv::DFT_REAL_OUTPUT);
	
	
	Eigen::MatrixXd tmpE;
	cv::cv2eigen(fIn, tmpE);
	double std = getStd(tmpE);
	double fac = 1.0 / std;
	hmap.filtered[idx] = fac * tmpE;
	std::cout << hmap.filtered[idx].row(0) << " \nllllllllllllll";

	std::ofstream file;
	file.open("c.txt");
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++) {
			//std::cout << a1[i * rhs_m + j] << "\b    ";
			file << hmap.filtered[idx](i,j) << " ";

		}
		std::cout << "\n";
		file << "\n";




	}
	file.close();


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
	
	mars_step_four(0);
	mars_step_five(0);
	/*for (int i = 0; i < hmap.N; i++) {
		mars_step_four();

		if (filt.option == 1) {
			mars_step_five(i);
		}
		else {
			hmap.filtered[i] = hmap.raw[i];
		}
	}*/
	
	
}