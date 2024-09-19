#include "Task.h"
#include "Matrix.h"
#include "Fredholm.h"
#define N 1000000

void Task::Count(double (*K)(double x, double t), double (*f)(double x), const double alpha, const double beta, vector<double>& cn, vector<double>& Un, vector<double>& tn, const int k)
{
	vector<double> c = { (beta - alpha)/6.0, 4.0*(beta - alpha)/6.0, (beta - alpha)/6.0 };
	vector<double> U(3);
	vector<double> t = { alpha, (alpha + beta)/2.0, beta };
	//cout << "c = " << c; 
	cn[1+2*k] = c[1]; cn[1+2*k+1] = c[2];
	tn[1+2*k] = t[1]; tn[1+2*k+1] = t[2];

	vector<vector<double>> m(3, vector<double>(3));
	for (int i = 0; i < 3; i++)
	{
		m[i][i] = 1;
		for (int j = 0; j < 3; j++)
		{
			m[i][j] -= c[j] * K(t[i], t[j]);
		}
	}
	Matrix M(m);
	//cout << M;
	vector<double> F(3);
	for (int i = 0; i < 3; i++) F[i] = f(t[i]);
	U = M.Gauss(F);
	//cout << "U = " << U;
	if (k == 0) Un[0] = U[0];
	Un[1+2*k] = U[1]; Un[1+2*k+1] = U[2];
}

void Task::Count2(double (*K)(double x, double t), double (*f)(double x), const double alpha, const double beta, vector<double>& cn, vector<double>& Un, vector<double>& tn, const int k, Matrix& Mat_glob, vector<double>& F_glob)
{
	vector<double> c = { (beta - alpha) / 6.0, 4.0 * (beta - alpha) / 6.0, (beta - alpha) / 6.0 };
	vector<double> t = { alpha, (alpha + beta) / 2.0, beta };
	//cout << "c = " << c;
	cn[1+2*k] += c[1]; cn[1+2*k+1] += c[2];
	tn[1+2*k] = t[1]; tn[1+2*k+1] = t[2];

	vector<double> F(3);
	for (int i = 0; i < 3; i++) F[i] = f(t[i]);
	for (int i = 0; i < 3; i++) F_glob[2*k+i] += F[i];

	for (int i = 0; i < 3; i++)
	{
		Mat_glob.a[2*k+i][2*k+i] += 1;
		for (int j = 0; j < 3; j++)
		{
			Mat_glob.a[2*k+i][2*k+j] -= c[j] * K(t[i], t[j]);
		}
	}
}

double Task::Integral(std::function<double(double)> f, double a, double b, double n, double epsilon)
{
	double I1 = 0.0, I2 = 0.0;
	double c, d, h;

	h = (b-a)/n;
	c = a; d = c + h;

	for (int i = 0; i < n; i++)
	{
		I2 += (d - c)/6.0 * (f(c) + 4.0*f((c + d)/2) + f(d));
		c = d; d += h;
	}
	/*do
	{
		I1 = I2; I2 = 0.0;
		n *= 2;
		h = (b - a) / n; c = a; d = c + h;
		for (int i = 0; i < n; i++)
		{
			I2 += (d - c)/6.0 * (f(c) + 4.0*f((c+d)/2) + f(d));
			c = d; d += h;
		}
	} while (fabs(I1 - I2) > epsilon && n <= N);*/

	return I2;
}

double Task::L2_norm_sqr(std::function<double(double)> g, const double a, const double b, const double epsilon)
{
	auto f = [&](double x) { return pow(g(x), 2); };
	const int p = 4;
	int n_global = 10;
	double n;
	double h = (b-a)/n_global;
	double MODIFIED_epsilon = epsilon * h / (b-a);
	double rh;
	double I = 0.0;
	double Sh, Sh2;
	double alpha, beta = a;
	while (fabs(b-beta) > epsilon)
	{
		alpha = beta; beta = std::min(beta + 2*(b-a)/n_global, b);

		n = 1;
		h = (b-a)/n_global;
		Sh = Integral(f, alpha, beta, n, epsilon);		// h
		Sh2 = Integral(f, alpha, beta, 2*n, epsilon);	// h/2
		rh = std::fabs(Sh - Sh2) / (1 - pow(2, 1 - p));
		while (rh > MODIFIED_epsilon)
		{
			Sh = Sh2;
			n *= 2;
			h /= 2;
			Sh2 = Integral(f, alpha, beta, 2*n, epsilon);
			rh = std::fabs(Sh - Sh2) / (1 - pow(2, 1 - p));
			MODIFIED_epsilon = epsilon * h / (b - a);
		}
		I += Sh;
		//cout << "I = " << I << endl;
	}
	return I;
}

void Task::main_func()
{
	Fredholm Eq = Fredholm17();
	double a = Eq.a, b = Eq.b;
	//double (*u)(double x) = Eq.u;
	double (*K)(double x, double t) = Eq.K;
	double (*f)(double x) = Eq.f;
	const double epsilon = 1e-4;
	int n = 1 + 2*1;
	int num_of_subsegments = (n-1)/2;

	vector<double> cn1, cn2;
	vector<double> Un1, Un2;
	vector<double> tn1, tn2;

	cn2.resize(n); Un2.resize(n); tn2.resize(n);
	cn2[0] = (b-a)/6.0; tn2[0] = a;
	double alpha = a, beta = a + (b-a)/ num_of_subsegments;
	
	for (int i = 0; i < num_of_subsegments; i++)
	{
		Count(K, f, alpha, beta, cn2, Un2, tn2, i);
		alpha = beta; beta += (b-a)/ num_of_subsegments;
	}
	//cout << Un2 << endl;
	std::function<double(double)> un, u2n;
	u2n = [&](double x) {
		double S(0.0); 
		for (int j = 0; j < cn2.size(); j++) S += cn2[j] * K(x, tn2[j]) * Un2[j];
		return S + f(x);
	};
	std::function<double(double)> diff;
	do
	{
		cn1 = std::move(cn2); Un1 = std::move(Un2); tn1 = std::move(tn2);
		un = [&](double x) {
			double S(0.0);
			for (int j = 0; j < cn1.size(); j++) S += cn1[j] * K(x, tn1[j]) * Un1[j];
			return S + f(x);
		};
		n = 1 + (n - 1) * 2;
		num_of_subsegments = (n - 1) / 2;
		Matrix Mat_glob(n, n);
		vector<double> F_glob(n);
		cn2.resize(n); Un2.resize(n); tn2.resize(n);
		double alpha = a, beta = a + (b - a) / num_of_subsegments;
		cn2[0] = (beta-alpha)/6.0; //???
		cn2[0] = 1.0/6.0; tn2[0] = a;
		for (int k = 0; k < num_of_subsegments; k++)
		{
			//Count(K, f, alpha, beta, cn2, Un2, tn2, k);
			Count2(K, f, alpha, beta, cn2, Un2, tn2, k, Mat_glob, F_glob);
			alpha = beta; beta += (b - a) / num_of_subsegments;
		}
		cout << "F_glob = " << F_glob << endl;
		cout << Mat_glob;
		Un2 = Mat_glob.Gauss(F_glob);
		//cout << "Un2 = " << Un2;
		u2n = [&](double x) {
			double S(0.0);
			for (int j = 0; j < cn2.size(); j++) S += cn2[j] * K(x, tn2[j]) * Un2[j];
			return S + f(x);
		};
		cout << un(0.25) << " " << un(1.0) << endl;
		cout << u2n(0.25) << " " << u2n(1.0) << endl << endl;
		diff = [&](double x) { 
			return un(x) - u2n(x); 
		};
		//cout << Integral(un, a, b, n, epsilon) << "\t" << Integral(u2n, a, b, n, epsilon) << endl;
		//cout << diff(0.2) << " " << diff(0.4) << " " << diff(0.6) << "" << diff(0.8) << " " << diff(1.0) << endl;
		cout << L2_norm_sqr(diff, a, b, epsilon) << endl;
	} 
	while ( L2_norm_sqr(diff, a, b, epsilon) > epsilon && 0 > 1);

	ofstream out("out.txt");
	for (auto& tt : tn1) out << tt << " "; out << endl;
	for (auto& tt : tn1) out << un(tt) << " "; out << endl;
	for (auto& tt : tn2) out << tt << " "; out << endl;
	for (auto& tt : tn2) out << u2n(tt) << " "; out << endl;
	system("python script.py");
}