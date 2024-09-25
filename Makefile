all: alg1 alg2 alg_e1

alg1: alg1.c vbf.c
	gcc alg1.c vbf.c -lm -o check_lin_eq_2x_uniform_3to1.o

alg2: alg2.c vbf.c
	gcc alg2.c vbf.c -lm -o find_uniform_3to1_lin_eq_to_3to1.o

alg_e1: alg_e1.c vbf.c
	gcc alg_e1.c vbf.c -lm -o lin_eq_to_self.o
