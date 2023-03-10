#pragma once
class NonLinearSystem {
public:
	NonLinearSystem() {};
	~NonLinearSystem() {};
	void setValue(const double& epsilon_, const int& max_iter_);
	
	
	double epsilon;
	int max_iter;
	int success;
};