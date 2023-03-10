#include<iostream>
#include<Eigen/core>
#include"mars/mars.h"
int main(){
	mars_algorithm::MARS mars;
	mars.mars_input(6.0, 3.0, 512, 256, 2,
		44, 3, 0.1,
		0.00001, 8,
		1,40);
	mars.solve();
	
} 
