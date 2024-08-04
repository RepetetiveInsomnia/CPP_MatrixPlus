#include "gtest/gtest.h"
#include "s21_matrix_oop.h"

void s21_print_matrix(S21Matrix &other) {
  for (int i = 0; i < other.GetRows(); i++) {
    for (int j = 0; j < other.GetCols(); j++) {
      printf("%.2f  ", other(i, j));
    }
    printf("\n");
  }
}

void s21_fill_matrix(S21Matrix &other, double key) {
  double var = key;
  for (int i = 0; i < other.GetRows(); i++) {
    for (int j = 0; j < other.GetCols(); j++) {
      other(i, j) = var;
      var++;
    }
  }
}
void s21_rand_matrix(S21Matrix &other) {
  for (int i = 0; i < other.GetRows(); i++) {
    for (int j = 0; j < other.GetCols(); j++) {
      other(i, j) = rand() % 10;
    }
  }
}
void s21_put_matrix(S21Matrix &other, double *array) {
  int counter = 0;
  int size = sizeof(array);

  for (int i = 0; i < other.GetRows() && counter < size; i++) {
    for (int j = 0; j < other.GetCols() && counter < size; j++) {
      other(i, j) = array[counter];
      array++;
    }
  }
}

TEST(Matrix, Constructors) {
  S21Matrix first;
  first(0, 0) = 99.99;
  S21Matrix second(1, 1);
  second(0, 0) = 99.99;
  EXPECT_EQ(first, second);
}

TEST(Matrix, Copy) {
  S21Matrix first(rand() % 10, rand() % 10);
  s21_rand_matrix(first);
  S21Matrix second(first);
  EXPECT_EQ(first, second);
}

TEST(Matrix, Move) {
  S21Matrix first(rand() % 10, rand() % 10);
  s21_rand_matrix(first);
  S21Matrix second(first);
  S21Matrix move(std::move(first));
  EXPECT_EQ(move, second);
  EXPECT_THROW(move = first, std::invalid_argument);
}

TEST(Matrix, SetterGetter) {
  int row = rand() % 10;
  int col = rand() % 10;
  S21Matrix second(row, col);
  EXPECT_EQ(second.GetCols(), col);
  EXPECT_EQ(second.GetRows(), row);
  second.SetRows(10);
  EXPECT_EQ(second.GetRows(), 10);
  second.SetCols(5);
  EXPECT_EQ(second.GetCols(), 5);
}

TEST(Matrix, Operations) {
  S21Matrix first(3, 3);
  double array[] = {1, 2, 4, 3, 3, 5, 2, 4, 4};
  s21_put_matrix(first, array);
  S21Matrix one(3, 3);
  double array_one[] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  s21_put_matrix(one, array_one);
  S21Matrix calc1(3, 3);
  S21Matrix calc2(3, 3);
  calc1 = first.InverseMatrix();
  calc2 = first.InverseMatrix();
  EXPECT_TRUE(calc1 == calc2);
  calc1.MulMatrix(first);
  EXPECT_EQ(one, calc1);
}

TEST(Matrix, Determinant) {
  S21Matrix other = S21Matrix();
  other(0, 0) = 214214321.4325452;
  double result = 0;
  result = other.Determinant();
  EXPECT_EQ(result, other(0, 0));
  other.SetCols(4);
  other.SetRows(4);
  double array[] = {1, 2, 3, 5, 6, 1, 12, 3, 0, 11, 2, 33, 7, 68, 9, 71};
  s21_put_matrix(other, array);
  double det = -6984;
  result = 0;
  result = other.Determinant();
  EXPECT_EQ(result, det);
  other.SetCols(80);
  other.SetRows(80);
  result = 0;
  result = other.Determinant();
  EXPECT_EQ(result, 0.0);

  for (int i = 0; i < other.GetRows(); i++) {
    for (int j = 0; j < other.GetCols(); j++) {
      other(i, j) = (random() % 1000) * 0.01;
    }
  }
  det = 0;
  result = 0;
  det = other.Determinant();
  result = other.Determinant();
  EXPECT_EQ(det, result);
}

TEST(Matrix, Operators) {
  S21Matrix first(3, 3);
  double array[] = {11.01, 22.02, 33.03, 44.04, 55.05,
                    66.06, 77.07, 88.08, 99.09};
  s21_put_matrix(first, array);
  S21Matrix sum(3, 3);
  sum = first + first;
  first *= 2.0;
  EXPECT_EQ(sum, first);
  for (int i = 0; i < 3; i++) {
    sum += first;
    sum *= 2.0;
    EXPECT_FALSE(sum == first);
    sum *= 0.5;
    sum -= first;
  }
  EXPECT_EQ(sum, first);
  sum = sum * sum;
  first = first * first;
  EXPECT_EQ(first, sum);
  first = first - sum;
  sum = sum - first;
  EXPECT_FALSE(first == sum);
}

TEST(Matrix, Errors) {
  EXPECT_THROW(S21Matrix errors(-123, 0), std::invalid_argument);
  EXPECT_THROW(S21Matrix errors(1, 0), std::invalid_argument);
  S21Matrix errorfirst(1, 55);
  S21Matrix errorsecond(2, 55);
  EXPECT_THROW(errorfirst -= errorsecond, std::invalid_argument);
  EXPECT_THROW(errorfirst += errorsecond, std::invalid_argument);
  EXPECT_THROW(errorfirst *= errorsecond, std::invalid_argument);
  S21Matrix first(4, 4);
  EXPECT_THROW(first(-12, 1) = 1, std::out_of_range);
  EXPECT_THROW(first(1, 123) = 14, std::out_of_range);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
