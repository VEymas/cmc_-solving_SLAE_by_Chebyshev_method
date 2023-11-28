#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace {
    const double EPS = std::numeric_limits<double>::epsilon();
}

//res m*k A m*n B n*k
std::vector<double> MatrixMultiplication(const std::vector<double>& matrix1, 
                    const std::vector<double>& matrix2, int m, int n, int k);

std::vector<double> MulMatrixByNumber(const std::vector<double>& matrix, double num);

std::vector<double> MatrixSubtraction(const std::vector<double>& matrix1, 
                    const std::vector<double>& matrix2);

double VectorNorm(const std::vector<double>& vector);

double VectorNorm2(const std::vector<double>& vector);

int Sign(double a);

std::vector<double> HouseholderMethod(const std::vector<double>& matrix, int size);