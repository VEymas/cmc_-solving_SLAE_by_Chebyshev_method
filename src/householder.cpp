#include "householder.h"

std::vector<double> MatrixMultiplication (const std::vector<double>& matrix1, 
                    const std::vector<double>& matrix2, int m, int n, int k) {
    std::vector<double> res_matrix(m * k, 0);
    for (int i = 0; i < m; ++i) {
        for (int q = 0; q < n; ++q) {
            for (int j = 0; j < k; ++j) {
                res_matrix[i * k + j] += matrix1[i * n + q] * matrix2[q * k + j];
            }
        }
    }
    return res_matrix;
} 

std::vector<double> MulMatrixByNumber(const std::vector<double>& matrix, double num) {
    int size = matrix.size();
    std::vector<double> res_matrix(size);
    for (int i = 0; i < size; ++i) {
        res_matrix[i] = matrix[i] * num;
    }
    return res_matrix;
}

std::vector<double> MatrixSubtraction(const std::vector<double>& matrix1, 
                    const std::vector<double>& matrix2) {
    int size = matrix1.size();
    std::vector<double> res_matrix(size);
    for (int i = 0; i < size; ++i) {
        res_matrix[i] = matrix1[i] - matrix2[i];
    }
    return res_matrix;
}

double VectorNorm(const std::vector<double>& vector) {
    int n = vector.size();
    double res{};
    for (int i = 0; i < n; ++i) {
        res += vector[i] * vector[i];
    }
    return sqrt(res);
}

double VectorNorm2(const std::vector<double>& vector) {
    int n = vector.size();
    double res{};
    for (int i = 0; i < n; ++i) {
        res += vector[i] * vector[i];
    }
    return res;
}

int Sign(double a) {
    if (a > EPS) {
        return 1;
    } else if (a < -EPS) {
        return -1;
    }
    return 0;
}

std::vector<double> HouseholderMethod(const std::vector<double>& matrix, int size) {
    std::vector<double> matrix(size * (size + 1));
    for (int i = 0; i < size * (size + 1); ++i) {
        matrix[i] = matrix[i];
    }
    for (int i = 0; i < size - 1; ++i) {
        std::vector<double> vector_s(size);
        for (int j = 0; j < size; ++j) {
            if (j < i) {
                vector_s[j] = 0;
            } else {
                vector_s[j] = matrix[j * (size + 1) + i];
            }
        }
        double b = Sign(matrix[i * (size + 1) + i]) * VectorNorm(vector_s);
        double p = sqrt(2 * (VectorNorm2(vector_s) + b * matrix[i * (size + 1) + i]));

        std::vector<double> vector_w(size);
        for(int j = 0; j < size; ++j) {
            if (j < i) {
                vector_w[j] = 0;
            } else if (j == i) {
                vector_w[j] = (matrix[i * (size + 1) + i] + b) / p;
            } else {
                vector_w[j] = (matrix[j * (size + 1) + i]) / p;
            }
        }
        std::vector<double> tmp_matrix = MatrixMultiplication(vector_w, matrix,
                             1, size, size + 1);
        tmp_matrix = MatrixMultiplication(vector_w, tmp_matrix, size, 1, size + 1);
        tmp_matrix = MulMatrixByNumber(tmp_matrix, 2);
        matrix = MatrixSubtraction(matrix, tmp_matrix);
    }
    std::vector<double> res(size, 0);
    for (int i = size - 1; i >= 0; --i) {
        res[i] = matrix[i * (size + 1) + size];
        for (int j = i + 1; j < size; ++j) {
            res[i] -= res[j] * matrix[i * (size + 1) + j];
        }
        res[i] /= matrix[i * (size + 1) + i];
    }
    return res;
}