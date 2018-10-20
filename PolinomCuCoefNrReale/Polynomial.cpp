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

