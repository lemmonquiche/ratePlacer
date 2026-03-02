# to do: make into proper makefile
#
all: gfl ratePlacer

gfl: src/get_frac_like.c
	gcc -o gfl src/get_frac_like.c -lm -g

ratePlacer: src/RatePlacer.c
	#gcc -fsanitize=address -o ratePlacer_noError_concat_bfgs src/RatePlacer.c src/bfgs.c -lm -g
	gcc -o ratePlacer_error_longInt_noPrint src/RatePlacer.c src/minfunc.c -lm -ggdb3


