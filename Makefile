CXX = clang++
CC = clang++
CPPFLAGS = -W -Wextra -Wall -Wconversion -Werror -std=c++20 -g -O1
default: cma
cma: cma.o
cma.o: ost.h

