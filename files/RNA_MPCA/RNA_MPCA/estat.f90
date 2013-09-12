! Eduardo F. P. Luz
! CAP/INPE
! 
! Modulo de funcoes estatisticas (media e desvio padrao)
! 
! Versao: 30_11_12a

MODULE estat

CONTAINS

! Calcula a media e desvio padrao dos resultados
!   nl -> numero de rodadas do experimento (numero de linhas)
!   nc -> numero de variaveis (numero de colunas)
!   vetE -> matriz nl por nc com os resultados dos experimentos
!   vetM -> vetor que contera o resultado da media de vetE
!   vetDP -> vetor que contera o resultado do desvio padrao de vetE
  SUBROUTINE mediaDP (nl,nc,vetE,vetM,vetDP)
    IMPLICIT NONE
    INTEGER :: i, j
    INTEGER, INTENT(IN) :: nl, nc
    DOUBLE PRECISION :: soma
    DOUBLE PRECISION, INTENT(IN), DIMENSION(nl,nc) :: vetE
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(nc) :: vetM, vetDP
    ! Calculando a media
    DO i=1, nc
      DO j=1, nl
        vetM(i) = vetM(i) + vetE(j,i) / nl
      END DO
    END DO
    ! Calculando o desvio padrao (se houver mais de 1 experimento)
    IF (nl .GT. 1) THEN
      DO i=1, nc
        soma = 0.0D0
        DO j=1, nl
          soma = soma + (vetE(j,i) - vetM(i))**2
        END DO
        vetDP(i) = SQRT((1.0D0/(nl-1.0D0))*soma)
      END DO
    ELSE
      vetDP = 0.0D0
    END IF
  END SUBROUTINE

!Fim do modulo
END MODULE estat
