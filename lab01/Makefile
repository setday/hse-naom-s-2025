all:
	g++ main.cpp -lgslcblas -std=c++23 -lm -o main -lopenblas -llapacke -lgsl -fsanitize=address

debug:
	g++-12 -std=c++23 -pedantic-errors -Wall -Wextra -Werror -fsanitize=undefined -fno-sanitize-recover=all -fsanitize=address -o main -g main.cpp -lgsl -llapacke
