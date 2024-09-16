#pragma once
#include <stdio.h>
#include <math.h>
#include <numeric>
#include <algorithm>
#include <cmath> 

#include <thread>
#include <omp.h>
#include <typeinfo>
#include <chrono>

#include <iostream>
#include <fstream>

#include <vector>
#include <algorithm>
#include <functional>

using namespace std;

ostream& operator<<(ostream& cout, const vector<double>& b);

class Matrix;
class Fredholm;

double l2_norm_square(vector<double>& x);

Fredholm Fredholm11();
Fredholm Fredholm12();
Fredholm Fredholm17();