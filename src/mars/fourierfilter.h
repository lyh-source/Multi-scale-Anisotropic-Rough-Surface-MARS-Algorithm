#pragma once

class FourierFilter {
public:
	FourierFilter() {};
	~FourierFilter() {};
	void setValue(const int& option_, const int& cutoff_);
	int option, cutoff;
};