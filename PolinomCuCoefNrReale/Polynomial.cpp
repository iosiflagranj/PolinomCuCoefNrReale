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
	Polynomial p = *this;
	for (int i = 0; i < p.m_coeff.size(); i++) {
		p.m_coeff[i] *= -1;
	}
	return p;
}

Polynomial & Polynomial::operator+=(const Polynomial & p)
{
	*this = *this + p;
	return *this;
}

Polynomial & Polynomial::operator+=(const double & x)
{
	this->m_coeff[this->m_coeff.size() - 1] += x;
	return *this;
}

Polynomial & Polynomial::operator-=(const Polynomial & p)
{
	*this = *this - p;
	return *this;
}

Polynomial & Polynomial::operator-=(const double & x)
{
	this->m_coeff[this->m_coeff.size() - 1] -= x;
	return *this;
}

Polynomial & Polynomial::operator*=(const Polynomial & p)
{
	*this = *this * p;
	return *this;
}

Polynomial & Polynomial::operator*=(const double & x)
{
	*this = *this * x;
	return *this;
}

Polynomial & Polynomial::operator/=(const Polynomial & p)
{
	*this = *this / p;
	return *this;
}

Polynomial & Polynomial::operator/=(const double & x)
{
	*this = *this / x;
	return *this;
}

Polynomial & Polynomial::operator%=(const Polynomial & p)
{
	*this = *this % p;
	return *this;
}

Polynomial & Polynomial::operator%=(const double & x)
{
	*this = *this % x;
	return *this;
}

double Polynomial::operator()(const double & x)
{
	double r = this->m_coeff[this->m_coeff.size() - 1];
	for (int i = 0; i < this->m_coeff.size() - 1; i++) {
		r += this->m_coeff[i] * pow(x, this->m_degree - i);
	}

	return r;
}

double & Polynomial::operator[](const int & deg)
{
	return this->m_coeff.at(this->m_degree - deg);
}

Polynomial::operator std::string() const
{
	Polynomial p = *this;
	return p.ToString();
}

std::string Polynomial::ToString()
{
	std::string string = "";
	if (this->m_degree == 0 && this->m_coeff[0] == 0) {
		string = "0";
		return string;
	}

	// The constant term
	if (this->m_coeff[this->m_coeff.size() - 1] != 0) {
		string += std::to_string(this->m_coeff[this->m_coeff.size() - 1]);
	}

	for (int i = this->m_coeff.size() - 2; i >= 0; i--) {
		if (this->m_coeff[i] == 0) {
			// Do nothing
		}
		else {
			if (this->m_coeff[i] < 0) {
				if (string == "")
					string += "-";
				else
					string += " - ";

				string += std::to_string(abs(this->m_coeff[i]));
			}

			else {
				if (string != "") {
					string += " + ";
				}
				string +=std::to_string (this->m_coeff[i]); 
			}

			if (this->m_degree - i > 1) {
				string += "*x^"; string += std::to_string(this->m_degree - i);
			}
			else {
				string += "*x";
			}
		}
	}

	return string;
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
	Polynomial r, q;
	q.m_coeff[0] = x;
	q.m_degree = 0;

	r = q % p;

	return r;
}

Polynomial operator%(const Polynomial & p, double x)
{
	Polynomial r;
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

bool operator==(const Polynomial & a, const Polynomial & b)
{
	if (a.m_degree == b.m_degree) {
		for (int i = 0; i < a.m_coeff.size(); i++) {
			if (a.m_coeff[i] != b.m_coeff[i]) {
				return false;
			}
		}
		return true;
	}
	return false;
}

bool operator==(const double & x, const Polynomial & p)
{
	return (p.m_degree == 0 && p.m_coeff[0] == x);
}

bool operator==(const Polynomial & p, const double & x)
{
	return (p.m_degree == 0 && p.m_coeff[0] == x);
}

bool operator!=(const Polynomial & a, const Polynomial & b)
{
	if (a.m_degree == b.m_degree) {
		for (int i = 0; i < a.m_coeff.size(); i++) {
			if (a.m_coeff[i] != b.m_coeff[i]) {
				return true;
			}
		}
		return false;
	}
	return true;
}

bool operator!=(const double & x, const Polynomial & p)
{
	return (p.m_degree != 0 || p.m_coeff[0] != x);
}

bool operator!=(const Polynomial & p, const double & x)
{
	return (p.m_degree != 0 || p.m_coeff[0] != x);
}

std::istream & operator>>(std::istream & stream, Polynomial & p)
{
	std::string string;
	char ch; stream.get(ch);
	while (ch != '\n' && !stream.eof()) {
		string += ch;
		stream.get(ch);
	}

	// remove the spaces from the input
	for (int i = 0; i < string.length(); i++) {
		if (string[i] == ' ') {
			string.erase(string.begin() + i);
		}
	}

	p.m_coeff.clear();
	p.m_degree = -1;


	int start_pos = 0, end_pos, start_coeff, end_coeff, curr_degree = 0, sign, substr_degree;
	double coeff;
	std::string substring;
	bool finished = false;
	while (!finished) {
		
		//find + or - or '\0'
		end_pos = 0;
		for (int i = start_pos + 1; i < string.length(); i++) {
			if (string[i] == '+' || string[i] == '-') {
				end_pos = i;
				break;
			}
		}
		if (end_pos == 0) {
			// end of string
			end_pos = string.length();
			finished = true;
		}

		substring = string.substr(start_pos, end_pos - start_pos);

		//get the sign
		if (substring[0] == '-') {
			sign = -1;
		}
		else {
			sign = 1;
		}

		//get the coefficient
		start_coeff = 0;
		if (substring[0] == '+' || substring[0] == '-') {
			start_coeff = 1;
		}

		end_coeff = substring.length();

		for (int j = start_coeff + 1; j < substring.length(); j++) {
			if (strchr("0123456789.", substring[j]) == NULL) {
				end_coeff = j;
				break;
			}
		}

		coeff = std::stod(substring.substr(start_coeff, end_coeff));


		//get degree
		if (end_coeff == substring.length()) {
			substr_degree = 0;
		}
		else {
			// look for '^'
			if (end_coeff + 2 >= substring.length()) {
				substr_degree = 1;
			}
			else {
				substr_degree = std::stoi(substring.substr(end_coeff + 3, substring.length() - (end_coeff + 3)));
			}
		}

		// if the difference of degree is not 1, then fill the coefficients with 0

		while (p.m_degree < substr_degree - 1) {
			p.m_coeff.push_back(0);
			p.m_degree++;
		}

		p.m_coeff.push_back(sign * coeff);
		p.m_degree++;


		start_pos = end_pos;
	}


	//reverse p
	std::reverse(p.m_coeff.begin(), p.m_coeff.end());
	return stream;
}

std::ostream & operator<<(std::ostream & stream, const Polynomial & p)
{
	Polynomial q = p;
	stream << q.ToString();
	return stream;
}


