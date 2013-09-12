# Compilador
CC=mpif90

# Opcoes de compilacao
CFLAG= -c -O2

# Opcoes de otimizacao
CFLAGOPT= -O2

mpca: estat.o ModuloRNA.o funcoes.o mpca.o
	$(CC) $(CFLAGOPT) -o mpca estat.o ModuloRNA.o funcoes.o mpca.o /home/juliana/nr2f90/mpif90/lib/libnr2f90.a

mpca.o: mpca.f90
	$(CC) $(CFLAG) mpca.f90 -I /home/juliana/nr2f90/mpif90/mod/

estat.o: estat.f90
	$(CC) $(CFLAG) estat.f90

funcoes.o: funcoes.f90
	$(CC) $(CFLAG) funcoes.f90

ModuloRNA.o: ModuloRNA.f90
	$(CC) $(CFLAG) ModuloRNA.f90

clean:
	rm -rf *.*~ Makefile~ *.o *.mod mpca
