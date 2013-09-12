!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RNA: PERCEPTRON DE MULTIPLAS CAMADAS (MLP)- AUTOCONFIGURAVEL 
! USANDO MULTIPLE PARTICLE COLLISION ALGORITHM (MPCA)
! 
! Desenvolvido por: Juliana A Anochi (CAP/INPE)
!					Sabrina Sambatti (CAP/INPE)
! Email: juliana.anochi@lac.inpe.br
!        sabrinabms@gmail.com          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


PROGRAM TREINAMENTO 

  USE ModuloRNA
  IMPLICIT NONE

  DOUBLE PRECISION, DIMENSION(16,1) :: parametros
  DOUBLE PRECISION a
  CHARACTER*40 caminho, arquivo, path
  INTEGER i

  caminho = './dados/'
  arquivo = 'entrada.txt'
  path = TRIM(caminho)//TRIM(arquivo)

  OPEN(UNIT = 30, FILE = path, STATUS = "old")
  DO i = 1,16
	READ(30, *) parametros(i, 1)
	WRITE(*,*) parametros(i,1)
  END DO
  CLOSE(30)

  a = Rede_Neural_BP (parametros)	

  WRITE(*,*)'Valor da funcao objetivo--> ', a

END PROGRAM 
