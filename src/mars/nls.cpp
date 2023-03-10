#include"nls.h"
void NonLinearSystem::setValue(const double& epsilon_, const int& max_iter_) {
	epsilon = epsilon_;
	max_iter = max_iter_;
}