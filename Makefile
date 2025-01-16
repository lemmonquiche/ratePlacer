# to do: make into proper makefile
#
all: gfl ratePlacer

gfl: src/get_frac_like.c
	gcc -o gfl src/get_frac_like.c -lm -g

ratePlacer: src/RatePlacer.c
	gcc -o ratePlacer_noError_concat_restart src/RatePlacer.c -lm -g


