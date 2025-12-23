#include <bits/stdc++.h>

using namespace std;

class Matrix {
public:
    vector<vector<double>> data;
    int rows, cols;

    Matrix(int r, int c) : rows(r), cols(c) {
        data.resize(r, vector<double>(c, 0.0));
    }

    Matrix operator*(const Matrix& other) const {
        Matrix result(rows, other.cols);
        for(int i = 0; i < rows; i++) {
            for(int j = 0; j < other.cols; j++) {
                for(int k = 0; k < cols; k++) {
                    result.data[i][j] += data[i][k] * other.data[k][j];
                }
            }
        }
        return result;
    }

    Matrix transpose() const {
        Matrix result(cols, rows);
        for(int i = 0; i < rows; i++) {
            for(int j = 0; j < cols; j++) {
                result.data[j][i] = data[i][j];
            }
        }
        return result;
    }

    Matrix operator+(const Matrix& other) const {
        Matrix result(rows, cols);
        for(int i = 0; i < rows; i++) {
            for(int j = 0; j < cols; j++) {
                result.data[i][j] = data[i][j] + other.data[i][j];
            }
        }
        return result;
    }

    Matrix operator-(const Matrix& other) const {
        Matrix result(rows, cols);
        for(int i = 0; i < rows; i++) {
            for(int j = 0; j < cols; j++) {
                result.data[i][j] = data[i][j] - other.data[i][j];
            }
        }
        return result;
    }

    Matrix inverse() const {
        // Упрощенная реализация для диагональных матриц
        Matrix result(rows, cols);
        for(int i = 0; i < min(rows, cols); i++) {
            if(fabs(data[i][i]) > 1e-10) {
                result.data[i][i] = 1.0 / data[i][i];
            }
        }
        return result;
    }

    static Matrix identity(int n) {
        Matrix result(n, n);
        for(int i = 0; i < n; i++) {
            result.data[i][i] = 1.0;
        }
        return result;
    }
};
