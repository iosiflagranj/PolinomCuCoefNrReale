#pragma once
#include <vector>
#include <stdexcept>
#include <cmath>
#include <string>
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

	//Compund operators
	Polynomial& operator+=(const Polynomial &);
	Polynomial& operator+=(const double &);

	Polynomial& operator-=(const Polynomial &);
	Polynomial& operator-=(const double &);

	Polynomial& operator*=(const Polynomial &);
	Polynomial& operator*=(const double &);

	Polynomial& operator/=(const Polynomial &);
	Polynomial& operator/=(const double &);

	Polynomial& operator%=(const Polynomial &);
	Polynomial& operator%=(const double &);

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

	friend Polynomial operator%(const Polynomial &, const Polynomial &);
	friend Polynomial operator%(double, const Polynomial &);
	friend Polynomial operator%(const Polynomial &, double);

	friend Polynomial operator^(const Polynomial &, int);

	//Logical Operators
	friend bool operator==(const Polynomial &, const Polynomial &);
	friend bool operator==(const double &, const Polynomial &);
	friend bool operator==(const Polynomial &, const double &);

	friend bool operator!=(const Polynomial &, const Polynomial &);
	friend bool operator!=(const double &, const Polynomial &);
	friend bool operator!=(const Polynomial &, const double &);

	//Other Operators
	double operator()(const double &);
	double& operator[](const int &);

	//Conversion Operator
	operator std::string() const;

	// I/O Operators
	friend std::istream & operator>>(std::istream &, Polynomial &);
	friend std::ostream & operator<<(std::ostream &, const Polynomial &);

	//Methods
	std::string ToString();

protected:
private:
	std::vector<double> m_coeff;
	int m_degree;
};

