#include"nls.h"
#include<iostream>
#include"fsolve.h"
//#include<Eigen/NonLinearOptimization>
void NonLinearSystem::setValue(const double& epsilon_, const int& max_iter_) {
	epsilon = epsilon_;
	max_iter = max_iter_;
	ff = new double[3 * 3];
	xx = new double[3 * 3];
}
Eigen::MatrixXd NonLinearSystem::nls_init(const Eigen::MatrixXd& rpq) {
	int n = rpq.rows();
	int m = rpq.cols();

	Eigen::MatrixXd c = Eigen::MatrixXd::Zero(n, m);

	//initial guess
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			c(i, j) = rpq(i, j) / ((double)(n - i) * (double)(m - j));
		}
	}

	double c2 = 0;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			c2 = c2 + pow(c(i, j), 2);
		}
	}

	double s2 = rpq(0, 0) / c2;
	Eigen::MatrixXd a0 = sqrt(s2) * c;


	xx = new double[n * m];
	for (int p = 0; p < n; p++) {
		for (int q = 0; q < m; q++) {
			xx[p * m + q] = a0(p, q);
		}
	}

	std::cout << a0<<"\n";
	return a0;


}

int rhs_n, rhs_m;
Eigen::MatrixXd corr_rhs;
Eigen::MatrixXd NonLinearSystem::nls_rhs(const Eigen::MatrixXd& alpha, const Eigen::MatrixXd& rpq) {

	int n = rpq.rows();
	int m = rpq.cols();
	rhs_n = n;
	rhs_m = m;
	corr_rhs = rpq;

	//std::cout << alpha << "\n";
	Eigen::MatrixXd f = Eigen::MatrixXd::Zero(n, m);

	//Assemble non-linear
	for (int p = 0; p < n; p++) {
		for (int q = 0; q < m; q++) {
			for (int k = 1; k < n - p + 1;k++) {
				for (int l = 1; l < m - q + 1; l++) {
					f(p, q) = f(p, q) + (alpha(k - 1, l - 1) * alpha(k - 1 + p, l - 1 + q));
				}
			}
		}
	}


	
	f = f - rpq;

	//delete[] ff;
	ff = new double[n * m];
	for (int p = 0; p < n; p++) {
		for (int q = 0; q < m; q++) {
			ff[p * m + q] = f(p, q);
			//std::cout <<p<<" "<<q << " oooooooooooooooooooooooooooooooooooooooooooooooo\n";
		}
	}
	

	
	std::cout << f << "\n";
	return f;

}






void nls_rhs2(int n_, double x[], double fx[]) {

	int n = rhs_n;
	int m = rhs_m;
	//std::cout << alpha << "\n";
	Eigen::MatrixXd f = Eigen::MatrixXd::Zero(n, m);
	for (int i = 0; i < n * m; i++)
		fx[i] = 0;
	//Assemble non-linear
	for (int p = 0; p < n; p++) {
		for (int q = 0; q < m; q++) {
			

			for (int k = 1; k < n - p + 1; k++) {
				for (int l = 1; l < m - q + 1; l++) {
					fx[p * m + q] = fx[p * m + q] + (x[(k - 1) * m + l - 1] * x[(k - 1 + p) * m + l - 1 + q]);
				}
			}
		}
	}


	
	for (int p = 0; p < n; p++) {
		for (int q = 0; q < m; q++) {
			fx[p * m + q] = fx[p * m + q]-corr_rhs(p,q);
			//std::cout << corr_rhs.size()<<"     ppppppppppppppppppppppppppppppp\n";
		}
	}


	//std::cout << f << "\n";

}
#include<string>
#include<fstream>
void r8vec2_print(int n, double a1[], double a2[], std::string title)
{
  

    std::cout << "\n";
	std::cout << title << "\n";
	std::cout << "\n";
  /*  for (int i = 0; i <= n - 1; i++)
    {


		std::cout  << i
            << ": "  << a1[i]
            << "  "  << a2[i] << "\n";
    }*/
	std::ofstream file;
	file.open("a.txt");
	for (int i = 0; i < rhs_n; i++)
	{
		for (int j = 0; j < rhs_m; j++) {
			std::cout << a1[i*rhs_m+j] << "\b    ";
			file << a1[i*rhs_m+j] << " ";
				
		}
		std::cout << "\n";
		file << "\n";



		
	}
	file.close();

    return;
}


Eigen::MatrixXd NonLinearSystem::nls_solve(const Eigen::MatrixXd& guess, int max_iter, double tol, const Eigen::MatrixXd& rhs) {
	Eigen::MatrixXd f = nls_rhs(guess, rhs);
	//std::cout << f << "\n";
	Eigen::MatrixXd sol;

	double* fx;
	int i;
	int info;
	int lwa;
	int n = guess.rows() * guess.cols();

	double* wa;
	double* x;

	lwa = (n * (3 * n + 13)) / 2;
	fx = new double[n];
	wa = new double[lwa];
	x = xx;

	std::cout << "\n";
	std::cout << "fsolve_test4():\n";
	std::cout << "  fsolve() solves a nonlinear system of "<<n<<" equations.\n";

	fx = ff;
	nls_rhs2(n, x, fx);
	//r8vec2_print(n, x, fx, "  inil X, FX");
	info = fsolve(nls_rhs2, n, x, fx, tol, wa, lwa);

	std::cout << "\n";
	std::cout << "  Returned value of INFO = " << info << "\n";
	r8vec2_print(n, x, fx, "  Final X, FX");
	//
	//  Free memory
	//
	sol = f;
	for (int i = 0; i < rhs_n; i++)
	{
		for (int j = 0; j < rhs_m; j++) {
			sol(i, j) = x[i * rhs_m + j];
		}
	}
	
	
	delete[] fx;
	delete[] wa;
	//delete[] x;
	coeff = sol;
	std::cout << " aaaaaaaaaaaaaaaaa\n";
	return sol;
}