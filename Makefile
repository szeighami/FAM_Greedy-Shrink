all: 
	g++ arr.cpp qsort.cpp lp/*.cpp /csproject/kdd/share/szeighami/mmr/glpk/lib/libglpk.a -I. -I/csproject/kdd/share/szeighami/mmr/glpk/include -lm -o run
