#pragma once
#include <vector>
#include <iostream>

class Polynomial
{
public:
	Polynomial(std::vector<double> = std::vector<double>(1, 0));
	~Polynomial();

	//Getters + Setters
	const int get_degree();
	void set_degree(const int&);

	//Unary operators
	Polynomial operator+();
	Polynomial operator-();

	//This will be removed after implementing the << operator
	void Print()
	{
		for (int i = 0; i < m_coeff.size(); i++) {
			std::cout << m_coeff[i] << " ";
		}
		std::cout << "\n";
	}
protected:
private:
	std::vector<double> m_coeff;
	int m_degree;
};

