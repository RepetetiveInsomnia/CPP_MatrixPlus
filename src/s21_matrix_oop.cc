#include "s21_matrix_oop.h"

#include <unistd.h>
// constructors + destructors
S21Matrix::S21Matrix() : rows_(1), cols_(1), matrix_(nullptr) {
  s21_create_matrix();
}
S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows_ <= 0) {
    throw std::invalid_argument("Invalid number of rows " +
                                std::to_string(rows_));
  }
  if (cols_ <= 0) {
    throw std::invalid_argument("Invalid number of cols " +
                                std::to_string(cols_));
  }
  matrix_ = nullptr;
  s21_create_matrix();
}
S21Matrix::S21Matrix(const S21Matrix &other)
    : rows_(other.rows_), cols_(other.cols_) {
  s21_create_matrix();
  s21_copy_m(other.matrix_);
}
S21Matrix::S21Matrix(S21Matrix &&other) {
  matrix_ = other.matrix_;
  rows_ = other.rows_;
  cols_ = other.cols_;

  other.matrix_ = nullptr;
  other.cols_ = 0;
  other.rows_ = 0;
}
S21Matrix::~S21Matrix() {
  if (matrix_ != nullptr) {
    s21_remove_matrix();
  }
  rows_ = 0;
  cols_ = 0;
}
// set + get
int S21Matrix::GetRows() const { return rows_; }
int S21Matrix::GetCols() const { return cols_; }
void S21Matrix::SetRows(int rows) {
  if (rows != rows_) {
    S21Matrix tmpMatrix(rows, cols_);

    for (int i = 0; i < rows && i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        tmpMatrix.matrix_[i][j] = matrix_[i][j];
      }
    }

    s21_remove_matrix();
    rows_ = tmpMatrix.rows_;
    matrix_ = tmpMatrix.matrix_;

    tmpMatrix.matrix_ = nullptr;
  }
}
void S21Matrix::SetCols(int cols) {
  if (cols != cols_) {
    S21Matrix tmpMatrix(rows_, cols);

    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols && j < cols_; j++) {
        tmpMatrix.matrix_[i][j] = matrix_[i][j];
      }
    }

    s21_remove_matrix();
    cols_ = tmpMatrix.cols_;
    matrix_ = tmpMatrix.matrix_;

    tmpMatrix.matrix_ = nullptr;
  }
}
// operation func
S21Matrix S21Matrix::Transpose() {
  S21Matrix result(cols_, rows_);
  for (int i = 0; i < rows_; ++i) {
    for (int32_t j = 0; j < cols_; ++j) {
      result.matrix_[j][i] = matrix_[i][j];
    }
  }

  return result;
}
S21Matrix S21Matrix::s21_minor(int row, int col) {
  int mem_row = 0;
  int mem_col = 0;
  S21Matrix result(cols_ - 1, rows_ - 1);
  for (int i = 0; i < rows_; i++) {
    if (i == row) continue;
    mem_col = 0;
    for (int j = 0; j < cols_; j++) {
      if (j == col) continue;
      result.matrix_[mem_row][mem_col] = matrix_[i][j];
      mem_col++;
    }
    mem_row++;
  }
  return result;
}
double S21Matrix::Determinant() {
  double result = 0.0;
  if (s21_is_square()) {
    if (rows_ == 1) {
      result = matrix_[0][0];
    } else if (rows_ == 2) {
      result = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
    } else if (rows_ < 6) {
      for (int j = 0; j < cols_; j++) {
        result += (*this)(0, j) * s21_minor(0, j).Determinant() * pow(-1, j);
      }
    } else {
      result = s21_gauss();
    }
    if (isnan(result)) {
      result = 0;
    }
  } else {
    throw std::invalid_argument("matrix is not square");
  }
  return result;
}
double S21Matrix::s21_gauss() {
  S21Matrix gauss(*this);
  int n = rows_;
  double res = 1.0;
  for (int i = 0; i < n; i++) {
    int maxRow = i;
    for (int j = i + 1; j < n; j++) {
      if (fabs(gauss(j, i)) > fabs(gauss(maxRow, i))) {
        maxRow = j;
      }
    }
    if (i != maxRow) {
      s21_swapRows(gauss.matrix_, maxRow, i);
      res *= -1.0;
    }
    for (int j = i + 1; j < n; j++) {
      double factor = gauss(j, i) / gauss(i, i);
      for (int k = i + 1; k < n; k++) {
        gauss.matrix_[j][k] -= factor * gauss(i, k);
      }
      gauss.matrix_[j][i] = 0.0;
    }
    res *= gauss(i, i);
  }
  return res;
}
void S21Matrix::s21_swapRows(double **matrix, int maxRow, int j) {
  for (int c = 0; c < cols_; c++) {
    double temp = matrix[j][c];
    matrix[j][c] = matrix[maxRow][c];
    matrix[maxRow][c] = temp;
  }
}
S21Matrix S21Matrix::CalcComplements() {
  bool return_code = true;
  S21Matrix other(rows_, cols_);
  if (rows_ <= 1 || cols_ <= 1) {
    return_code = false;
    throw std::out_of_range("Invalid argument");
  } else {
    if (s21_is_square()) {
      for (int i = 0; i < rows_; i++)
        for (int j = 0; j < cols_; j++) {
          double determinant = 0.0;
          S21Matrix minor = s21_minor(i, j);
          determinant = minor.Determinant();
          other.matrix_[i][j] = pow(-1, (i + j)) * determinant;
          minor.s21_remove_matrix();
        }
    } else {
      return_code = false;
      throw std::invalid_argument("matrix is not square");
    }
    return other;
  }
}
S21Matrix S21Matrix::InverseMatrix() {
  bool return_code = true;
  S21Matrix result(rows_, cols_);
  if (rows_ == 1 && cols_ == 1 && matrix_[0][0] != 0) {
    result.SetCols(1);
    result.SetRows(1);
    result.matrix_[0][0] = 1.0 / matrix_[0][0];

  } else {
    double det = 0.0;
    det = Determinant();
    if (det != 0) {
      S21Matrix calc(rows_, cols_);
      calc = CalcComplements();

      result.SetCols(calc.cols_);
      result.SetRows(calc.rows_);

      result = calc.Transpose();

      result *= 1 / det;

    } else {
      return_code = false;
    }
  }

  return result;
}
bool S21Matrix::EqMatrix(const S21Matrix &other) const {
  bool return_code = true;

  if (s21_is_eq(other)) {
    for (int i = 0; i < rows_ && return_code == true; i++) {
      for (int j = 0; j < cols_; j++) {
        if (fabs(matrix_[i][j] - other.matrix_[i][j]) > PRECISION) {
          return_code = false;
        }
      }
    }
  } else {
    throw std::invalid_argument("Matrixes of different size");
  }
  return return_code;
}
void S21Matrix::SumMatrix(const S21Matrix &other) {
  if (s21_is_eq(other)) {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        matrix_[i][j] += other.matrix_[i][j];
      }
    }
  } else {
    throw std::invalid_argument("Matrixes of different size");
  }
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  if (s21_is_eq(other)) {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        matrix_[i][j] -= other.matrix_[i][j];
      }
    }
  } else {
    throw std::invalid_argument("Matrixes of different size");
  }
}
void S21Matrix::MulNumber(const double &num) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] *= num;
    }
  }
}
void S21Matrix::MulMatrix(const S21Matrix &other) {
  if (s21_can_i_mult(other)) {
    S21Matrix tmp(rows_, other.cols_);
    for (int i = 0; i < tmp.rows_; i++) {
      for (int j = 0; j < tmp.cols_; j++) {
        for (int k = 0; k < cols_; k++) {
          tmp.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
        }
      }
    }
    *this = tmp;
  } else {
    throw std::invalid_argument("cant't mult matrixes");
  }
}
// operators
double &S21Matrix::operator()(int i, int j) {
  CheckIndexes(i, j);
  return matrix_[i][j];
}
bool S21Matrix::operator==(const S21Matrix &other) const {
  return EqMatrix(other);
}
bool S21Matrix::operator!=(const S21Matrix &other) { return ~EqMatrix(other); }
S21Matrix S21Matrix::operator-(const S21Matrix &other) {
  S21Matrix res = *this;
  res.SubMatrix(other);
  return res;
}
S21Matrix S21Matrix::operator+(const S21Matrix &other) {
  S21Matrix res = *this;
  res.SumMatrix(other);
  return res;
}
S21Matrix S21Matrix::operator*(const S21Matrix &other) {
  S21Matrix res = *this;
  res.MulMatrix(other);
  return res;
}
S21Matrix S21Matrix::operator*(const double &num) {
  S21Matrix res = *this;
  res.MulNumber(num);
  return res;
}
void S21Matrix::operator+=(const S21Matrix &other) { SumMatrix(other); }
void S21Matrix::operator-=(const S21Matrix &other) { SubMatrix(other); }
void S21Matrix::operator*=(const S21Matrix &other) { MulMatrix(other); }
void S21Matrix::operator*=(const double &num) { MulNumber(num); }
S21Matrix &S21Matrix::operator=(const S21Matrix &other) {
  S21Matrix(other.rows_, other.cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
  return *this;
}
// add func
void S21Matrix::s21_create_matrix() {
  if (rows_ <= 0) {
    throw std::invalid_argument("Invalid number of rows " +
                                std::to_string(rows_));
  }
  if (cols_ <= 0) {
    throw std::invalid_argument("Invalid number of cols " +
                                std::to_string(cols_));
  }
  matrix_ = new double *[rows_]();

  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_]();
  }
}

void S21Matrix::s21_remove_matrix() {
  for (int i = 0; i < rows_; i++) {
    delete[] matrix_[i];
  }
  delete[] matrix_;
  matrix_ = nullptr;
}
bool S21Matrix::s21_can_i_mult(const S21Matrix &other) {
  bool return_code = true;
  if ((cols_ != other.rows_) || (rows_ != other.cols_)) {
    return_code = false;
    // throw std::invalid_argument("cant't mult matrixes");
  }
  return return_code;
}
bool S21Matrix::s21_is_square() {
  bool return_code = true;
  if (cols_ != rows_) {
    return_code = false;
    // throw std::invalid_argument("matrix is not square");
  }
  return return_code;
}
bool S21Matrix::s21_is_eq(const S21Matrix &other) const {
  bool return_code = true;
  if ((cols_ != other.cols_) || (rows_ != other.rows_)) {
    return_code = false;
  }
  return return_code;
}

void S21Matrix::s21_copy_m(double **array) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = array[i][j];
    }
  }
}
void S21Matrix::CheckIndexes(int i, int j) {
  if (i < 0 || i > rows_ - 1) {
    throw std::out_of_range(
        "Invalid argument i - number of rows out of range [0:" +
        std::to_string(rows_ - 1) + "]");
  }
  if (j < 0 || j > cols_ - 1) {
    throw std::out_of_range(
        "Invalid argument j - number of cols out of range [0:" +
        std::to_string(cols_ - 1) + "]");
  }
}