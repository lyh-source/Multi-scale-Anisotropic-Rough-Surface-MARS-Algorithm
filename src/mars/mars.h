#pragma once
#include"hmap.h"
#include"nls.h"
#include"corr.h"
#include"fourierfilter.h"
namespace mars_algorithm {
	class MARS
	{
	public:
		MARS() {};
		~MARS() {};
		void mars_input(const double& Lx_, const double& Ly_, const int& Nx_, const int& Ny_, const int& N_, 
			const int& n_, const int& m_, const double& gamma_,
			const double& epsilon_, const int& max_iter_,
			const int& option_, const int& cutoff_);

		void mars_step_one();
		void mars_step_two();
		void mars_step_three();
		void mars_step_four(int idx);
		void mars_step_five(int idx);
		void solve();

		//var
		HeightMap hmap;
		NonLinearSystem nls;
		Correlation corr;
		FourierFilter filt;

	private:

	};

}