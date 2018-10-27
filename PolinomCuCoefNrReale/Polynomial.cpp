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

Polynomial operator/(const Polynomial & a, const Polynomial & b)
{
	Polynomial q, dividend = a, divisor = b;
	if (divisor.m_degree > dividend.m_degree) {
		return q; //Return 0, the division is impossible
	}

	q.m_degree = dividend.m_degree - divisor.m_degree;
	q.m_coeff.resize(q.m_degree + 1);

	for (int i = 0; i < q.m_coeff.size(); i++) {
		q.m_coeff[i] = dividend.m_coeff[0] / divisor.m_coeff[0];

		int degreeOffset = dividend.m_degree - divisor.m_degree;
		divisor.m_degree += degreeOffset;
		for (int j = 0; j < degreeOffset; j++) {
			divisor.m_coeff.push_back(0);
		}
		
		dividend = dividend - (divisor * q.m_coeff[i]);

		divisor.m_degree -= degreeOffset;
		for (int j = 0; j < degreeOffset; j++) {
			divisor.m_coeff.pop_back();
		}
		dividend.m_degree -= 1;
		dividend.m_coeff.erase(dividend.m_coeff.begin());
	}

	return q;
}

Polynomial operator/(double x, const Polynomial & p)
{
	Polynomial q;
	if (p.m_degree != 0) {
		return q;
	}
	
	q.m_coeff[0] = x / p.m_coeff[0];
	return q;
}

Polynomial operator/(const Polynomial & p, double x)
{
	Polynomial r = p;
	for (int i = 0; i < r.m_coeff.size(); i++) {
		r.m_coeff[i] /= x;
	}

	return r;
}


//TODO: To ask about remainder division
Polynomial operator%(const Polynomial & a, const Polynomial & b)
{
	Polynomial q, dividend = a, divisor = b;
	if (divisor.m_degree > dividend.m_degree) {
		return dividend; //Return dividend, the division is impossible
	}

	q.m_degree = dividend.m_degree - divisor.m_degree;
	q.m_coeff.resize(q.m_degree + 1);

	for (int i = 0; i < q.m_coeff.size(); i++) {
		q.m_coeff[i] = dividend.m_coeff[0] / divisor.m_coeff[0];

		int degreeOffset = dividend.m_degree - divisor.m_degree;
		divisor.m_degree += degreeOffset;
		for (int j = 0; j < degreeOffset; j++) {
			divisor.m_coeff.push_back(0);
		}

		dividend = dividend - (divisor * q.m_coeff[i]);

		divisor.m_degree -= degreeOffset;
		for (int j = 0; j < degreeOffset; j++) {
			divisor.m_coeff.pop_back();
		}
		dividend.m_degree -= 1;
		dividend.m_coeff.erase(dividend.m_coeff.begin());
	}

	return dividend;
}

Polynomial operator%(double x, const Polynomial & p)
{
	Polynomial r;

	if (p.m_degree > 0) {
		r.m_coeff[0] = x;
		return r;
	}

	if ((int)x != x || (int)p.m_coeff[0] != p.m_coeff[0]) {
		return r;
	}

	r.m_coeff[0] = (int)x % (int)p.m_coeff[0];
	return r;
}

Polynomial operator%(const Polynomial & p, double x)
{
	Polynomial r;

	if (p.m_degree > 0) {
		return r;
	}



	if ((int)x != x || (int)p.m_coeff[0] != p.m_coeff[0]) {
		return r;
	}

	r.m_coeff[0] = (int)p.m_coeff[0] % (int)x;
	return r;
}

Polynomial operator^(const Polynomial & p, int pow)
{
	Polynomial r = p;
	for (int i = 1; i < pow; i++) {
		r = r * p;
	}

	return r;
}

