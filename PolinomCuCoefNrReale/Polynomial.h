#pragma once
#include <vector>
#include <iostream>

class Polynomial
{
public:
	Polynomial(std::vector<double> = std::vector<double>(0));
	~Polynomial();

	const int get_degree();
	void set_degree(int);
protected:
private:
	std::vector<double> m_coeff;
	int m_degree;
};

