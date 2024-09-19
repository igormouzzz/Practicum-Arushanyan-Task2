#include "Header.h"
#include "Fredholm.h"

Fredholm Fredholm11()
{
	auto u = [](double x) { return x + exp(-x); };
	auto K = [](double x, double t) { return 0.5 * x * exp(t); };
	auto f = [](double x) { return exp(-x); };
	double a = 0.0, b = 1.0;
	return Fredholm(u, K, f, a, b);
}
Fredholm Fredholm12()
{
	auto u = [](double x) { return 1.0; };
	auto K = [](double x, double t) { return sin(x*t); };
	auto f = [](double x) { return std::fabs(x) < 1e-10 ? 0 : 1 + 1 / x * (cos(0.5*x) - 1); };
	double a = 0.0, b = 0.5;
	return Fredholm(u, K, f, a, b);
}
Fredholm Fredholm17()
{
	auto u = [](double x) { return 1.0 + 4.0/9.0 * x; };
	auto K = [](double x, double t) { return x*t*t; };
	auto f = [](double x) { return 1.0; };
	double a = 0.0, b = 1.0;
	return Fredholm(u, K, f, a, b);
}
Fredholm Fredholm18()
{
	auto u = [](double x) { return x; };
	auto K = [](double x, double t) { return 0.5*x*t; };
	auto f = [](double x) { return 5.0/6.0*x; };
	double a = 0.0, b = 1.0;
	return Fredholm(u, K, f, a, b);
}