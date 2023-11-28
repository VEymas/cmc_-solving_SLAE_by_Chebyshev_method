#include "householder.h"
#include <string>
#include <fstream>
#include <iterator>
#include <regex>
#include <sstream>
#include <ctime>

const std::regex comma(",");

std::vector<double> ReadDataFromCSV(std::string fname) {
    std::ifstream ifs(fname.c_str());
    std::vector<double> matrix;
    std::string line;
    while (ifs.good()) {
        std::getline(ifs, line);
        std::replace(line.begin(), line.end(), ',',  ' ');
        std::stringstream ss(line);
        std::string current;
        while (ss >> current) {
            try {
                const double d = std::stod(current);
                matrix.push_back(d);
            } catch (const std::exception& e) {
                std::cerr << "Error" << std::endl;
            }
        }
    }
    ifs.close();
    return matrix;
}

std::vector<double> GenerateRandomDistributedVector(int size) {
    std::vector<double> x(size);
    for (int i = 0; i < size; ++i) {
        double random_num = (double)rand() / RAND_MAX; // generate random num from [0, 1]
        double scaled_num = random_num * 2 - 1; // scale this num to [-1, 1]
        x[i] = scaled_num;
    }
    return x;
}

std::vector<double> MatrixToExtendedMatrix(const std::vector<double>& matrix, 
                    const std::vector<double>& vector, int size) {
    std::vector<double> ExtendedMatrix(size * (size + 1));
    for (int i = 0; i < size; ++i) {
        for(int j = 0; j < size + 1; ++j) {
            if (j == size) {
                ExtendedMatrix[i * (size + 1) + j] = vector[i];
            } else {
                ExtendedMatrix[i * (size + 1) + j] = matrix[i * size + j];
            }
        }
    }
    return ExtendedMatrix;
}

std::vector<double> MulMatrixOnVector(const std::vector<double>& matrix, 
                    const std::vector<double>& vector, int size) {
    std::vector<double> res(size, 0);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            res[i] += matrix[i * size + j] * vector[j];
        }
    }
    return res;
}

std::vector<double> AddVectors(const std::vector<double>& vector1,
                    const std::vector<double>& vector2) {
    if (vector1.size() != vector2.size()) {
        std::cout << "Sizes don't match";
    }
    std::vector<double> res(vector1.size());
    for (int i = 0; i < vector1.size(); ++i) {
        res[i] = vector1[i] + vector2[i];
    }
    return res;
}

std::vector<double> FindMismatch(const std::vector<double>& vector1,
                    const std::vector<double>& vector2) {
    int size = vector1.size();
    std::vector<double> mismatch(size);
    for (int i = 0; i < size; ++i) {
        mismatch[i] = sqrt((vector1[i] - vector2[i]) * (vector1[i] - vector2[i]));
    }
    return mismatch;
}

double FindMaxNorm(const std::vector<double>& vector) {
    double res = abs(vector[0]);
    for (int i = 0; i < vector.size(); ++i) {
        if (abs(vector[i]) > res) {
            res = abs(vector[i]);
        }
    }
    return res;
}
 
int main()
{
    std::srand(std::time(nullptr));
    const std::string fname = "../SLAU_var_5.csv";
    std::vector<double> matrix{ReadDataFromCSV(fname)};

    int size = sqrt(matrix.size());

    std::vector<double> x{GenerateRandomDistributedVector(size)};

    std::vector<double> vector_f = MulMatrixOnVector(matrix, x, size);

    vector_f = AddVectors(x, vector_f);

    std::vector<double> extended_matrix{MatrixToExtendedMatrix(matrix, vector_f, size)};

    clock_t start = clock();
    std::vector<double> res = HouseholderMethod(extended_matrix, size);
    clock_t end = clock();
    double seconds = (double)(end - start) / CLOCKS_PER_SEC;

    std::vector<double> mismatch{FindMismatch(res, x)};

    double MaxNorm = FindMaxNorm(mismatch);

    std::cout << "TIME - " << seconds << "s" <<std::endl;
    std::cout << "Max-Norm - " << MaxNorm << std::endl;
}