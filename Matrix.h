#include "Header.h"
#include "Task.h"
using row = vector<double>;

class Matrix
{
protected:
	vector<row> a;
	int N, M;
public:
	Matrix();
	Matrix(int N, int M);
	Matrix(vector<row> rows);
	Matrix(const Matrix& b);
	size_t GetN();
	size_t GetM();
	void SetValue(int i, int j, double value);
	Matrix& operator=(const Matrix& b);
	Matrix operator+(const Matrix& b);
	Matrix operator-(const Matrix& b);
	Matrix operator*(double k);
	Matrix operator*(const Matrix& b);
	vector<double> operator*(vector<double>& b);
	Matrix T();
	double Det2();
	Matrix Inv2();
	void Inversed(Matrix& Inv);

	void Move(vector<row> rows);
	
	vector<double> Gauss(vector<double> b);

	friend void Task::Count2(double (*K)(double x, double t), double (*f)(double x), const double a, const double b, vector<double>& cn, vector<double>& Un, vector<double>& tn, const int i, Matrix& Mat_glob, vector<double>& F_glob);
	friend ostream& operator<<(ostream& cout, const Matrix& b);
};
