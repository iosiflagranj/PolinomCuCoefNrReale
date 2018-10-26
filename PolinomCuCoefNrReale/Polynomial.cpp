#include "pch.h"
#include "Polynomial.h"

Polynomial::Polynomial(std::vector<double> coeff) :
	m_coeff(coeff),
	m_degree(coeff.size() - 1)
{
}

Polynomial::Polynomial(const Polynomial & p)
{
	this->m_coeff = p.m_coeff;
	this->m_degree = p.m_degree;
}

Polynomial::~Polynomial()
{
}

const int Polynomial::get_degree()
{
	return m_degree;
}

Polynomial Polynomial::operator+()
{
	Polynomial p = *this;
	return p;
}

Polynomial Polynomial::operator-()
{
	/*for (int i = 0; i < this->m_coeff.size(); i++) {
		this->m_coeff[i] *= -1;
	}
	return *this;*/

	Polynomial p = *this;
	for (int i = 0; i < p.m_coeff.size(); i++) {
		p.m_coeff[i] *= -1;
	}
	return p;
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

Polynomial operator-(const Polynomial & a, const Polynomial & b)
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
			c.m_coeff[i] -= b.m_coeff[k--];
		}
	}

	return c;
}

Polynomial operator-(double x, const Polynomial & a)
{
	Polynomial b = a;
	for (int i = 0; i < b.m_coeff.size(); i++) {
		b.m_coeff[i] = -b.m_coeff[i];
		if (i == b.m_coeff.size() - 1) {
			b.m_coeff[i] += x;
		}
	}

	return b;
}

Polynomial operator-(const Polynomial & a, double x)
{
	Polynomial b = a;
	b.m_coeff[b.m_coeff.size() - 1] -= x;

	return b;
}

Polynomial operator*(const Polynomial & a, const Polynomial & b)
{
	Polynomial c;
	int deg = a.m_degree + b.m_degree;
	c.m_coeff.resize(deg + 1);
	c.m_degree = deg;

	std::cout << "deg: " << c.m_degree << "\n\n";

	for (int i = 0; i < c.m_coeff.size(); i++) {
		// degree of x: m_degree - i
		// we need: a.m_degree - j + b.m_degree - k = m_deg - i

		c.m_coeff[i] = 0;
		for (int j = 0; j < a.m_coeff.size(); j++) {
			// j - index of element in a
			// k - index of elemnt in b

			int k = i - c.m_degree + a.m_degree - j + b.m_degree;


			if (k >= 0 && k <= b.m_degree) {
				c.m_coeff[i] += a.m_coeff[j] * b.m_coeff[k];
			}
		}
	}

	return c;

}

Polynomial operator*(double x, const Polynomial & p)
{
	Polynomial r = p;
	for (int i = 0; i < r.m_coeff.size(); i++) {
		r.m_coeff[i] *= x;
	}

	return r;
}

Polynomial operator*(const Polynomial & p, double x)
{
	Polynomial r = p;
	for (int i = 0; i < r.m_coeff.size(); i++) {
		r.m_coeff[i] *= x;
	}

	return r;
}
