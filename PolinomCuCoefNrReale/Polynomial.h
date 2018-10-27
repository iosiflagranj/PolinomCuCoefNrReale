#pragma once
#include <vector>
#include <iostream>

class Polynomial
{
public:
	Polynomial(std::vector<double> = std::vector<double>(1, 0));
	Polynomial(const Polynomial &);
	~Polynomial();

	//Getters + Setters
	const int get_degree();

	//Unary operators
	Polynomial operator+();
	Polynomial operator-();

	//Binary Operators
	friend Polynomial operator+(const Polynomial &, const Polynomial &);
	friend Polynomial operator+(double, const Polynomial &);
	friend Polynomial operator+(const Polynomial &, double);

	friend Polynomial operator-(const Polynomial &, const Polynomial &);
	friend Polynomial operator-(double, const Polynomial &);
	friend Polynomial operator-(const Polynomial &, double);

	friend Polynomial operator*(const Polynomial &, const Polynomial &);
	friend Polynomial operator*(double, const Polynomial &);
	friend Polynomial operator*(const Polynomial &, double);

	friend Polynomial operator/(const Polynomial &, const Polynomial &);
	friend Polynomial operator/(double, const Polynomial &);
	friend Polynomial operator/(const Polynomial &, double);

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

