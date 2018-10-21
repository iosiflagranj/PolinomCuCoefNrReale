#include "pch.h"
#include "Polynomial.h"

Polynomial::Polynomial(std::vector<double> coeff) :
	m_coeff(coeff),
	m_degree(coeff.size()-1)
{
}

Polynomial::~Polynomial()
{
}

const int Polynomial::get_degree()
{
	return m_degree;
}

void Polynomial::set_degree(const int& degree)
{
	m_degree = degree;
}

Polynomial Polynomial::operator+()
{
	return *this;
}

Polynomial Polynomial::operator-()
{
	for (int i = 0; i < this->m_coeff.size(); i++) {
		this->m_coeff[i] *= -1;
	}
	return *this;
}

Polynomial operator+(const Polynomial & a, const Polynomial & b)
{
	Polynomial c;
	int maxDeg = a.m_degree > b.m_degree ? a.m_degree : b.m_degree;
	c.m_degree = maxDeg;
	c.m_coeff.resize(maxDeg + 1);
	
	int j = a.m_coeff.size() - 1; int k = b.m_coeff.size() - 1;
	for (int i = c.m_coeff.size() - 1; i >= 0; i--) {
		c.m_coeff[i] = 0;
		if (j >= 0) {
			c.m_coeff[i] += a.m_coeff[j--];
		}
		if (k >= 0) {
			c.m_coeff[i] += b.m_coeff[k--];
		}
	}

	return c;
}

Polynomial operator+(double x, const Polynomial & a)
{
	Polynomial b = a;
	for (int i = 0; i < b.m_coeff.size(); i++) {
		b.m_coeff[i] += x;
	}

	return b;
}

Polynomial operator+(const Polynomial & a, double x)
{
	Polynomial b = a;
	for (int i = 0; i < b.m_coeff.size(); i++) {
		b.m_coeff[i] += x;
	}

	return b;
}
