GCC=gcc
SRC=s21_matrix_oop.cc
TEST_SRC=s21_matrix_test.cc
OBJ=s21_matrix_oop.o
CFLAGS=--std=c++17 -lstdc++
TESTFLAGS=-lgtest -lgcov -lm
GCOVFLAGS=--coverage
HTML=lcov -t test -o rep.info --no-external -c -d ./
OS = $(shell uname)

all: clean test

clean:
	rm -rf *.o *.a *.so *.out *.info *.gcda *.gcno *.gch rep.info *.html *.css test report *.txt *.dSYM

test:  clean
	$(GCC) $(CFLAGS) -g $(TEST_SRC) $(SRC) -lgtest -pthread -o test
	./test
s21_matrix_oop.a: clean
	$(GCC) $(CFLAGS) -c $(SRC)
	ar rcs s21_matrix_oop.a $(OBJ)
	ranlib s21_matrix_oop.a

main: s21_matrix_oop.a
	$(GCC) s21_matrix_oop.a main.cc $(CFLAGS) $(CFLAGS) -o main.out
	leaks --atExit -- ./main.out
	#./main.out

gcov_report: clean
	$(GCC) $(CFLAGS) -g $(TEST_SRC) $(GCOVFLAGS) $(SRC) -lgtest -pthread -o test
	./test
	rm -rf s21_matrix_test.gcda
	rm -rf s21_matrix_test.gcno
	$(HTML)
	genhtml -o report rep.info
	open report/index.html

check: test
	clang-format -style="{BasedOnStyle: Google}" -i *.cc *.h
	cppcheck --enable=all --suppress=missingIncludeSystem --inconclusive --check-config $(SRC) *.h
ifeq ($(OS), Darwin)
	leaks --atExit -- ./test
else
	CK_FORK=no valgrind --vgdb=no --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=RESULT_VALGRIND.txt ./test
endif