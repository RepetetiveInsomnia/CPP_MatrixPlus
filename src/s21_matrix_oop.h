#ifndef MATRIX_OOP_H_
#define MATRIX_OOP_H_
#include <cmath>
#include <iostream>
#define PRECISION 10E-7

typedef enum { BIN_OK = 0, BIN_CALCULATION_ERROR } bin_operation_result;
class S21Matrix {
 public:
  // constructors + destructors
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other);
  ~S21Matrix();
  // set + get
  int GetRows() const;
  int GetCols() const;
  void SetRows(int rows);
  void SetCols(int cols);
  // operation func
  bool EqMatrix(const S21Matrix& other) const;
  void MulNumber(const double& num);
  void MulMatrix(const S21Matrix& other);
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  S21Matrix Transpose();
  double Determinant();
  S21Matrix CalcComplements();
  S21Matrix InverseMatrix();
  // operators
  double& operator()(int i, int j);
  bool operator==(const S21Matrix& other) const;
  bool operator!=(const S21Matrix& other);
  void operator+=(const S21Matrix& other);
  void operator-=(const S21Matrix& other);
  void operator*=(const S21Matrix& other);
  void operator*=(const double& num);
  S21Matrix operator+(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix operator*(const S21Matrix& other);
  S21Matrix operator*(const double& num);
  S21Matrix& operator=(const S21Matrix& other);
  //  add func
  void s21_create_matrix();
  void s21_remove_matrix();
  bool s21_can_i_mult(const S21Matrix& other);
  bool s21_is_square();
  bool s21_is_eq(const S21Matrix& other) const;
  void s21_copy_m(double** array);
  S21Matrix s21_minor(int row, int col);
  void CheckIndexes(int i, int j);
  double s21_gauss();
  void s21_swapRows(double** matrix, int maxRow, int j);

 private:
  int rows_;
  int cols_;
  double** matrix_;
};
#endif  // MATRIX_OOP_H_