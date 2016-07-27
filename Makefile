#compilar
compandlink:
	f2py -c -m lecture lecture.f90   
	gfortran -Ofast -fexternal-blas -fblas-matmul-limit=40 -faggressive-function-elimination -o menu menu.f90 -llapack -lblas

