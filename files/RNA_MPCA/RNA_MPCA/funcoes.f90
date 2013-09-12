! Eduardo F. P. Luz
! CAP/INPE
! Versao: 30_11_12a

MODULE funcoes

CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Subrotina para definicao automatica dos limites
  ! Variaveis: nome da funcao, em character*16 (nomeFuncao)
  !            limite inferior, em double precision (linf)
  !            limite superior, em double precision (lsup)
  ! 
  ! Funcoes atualmente suportadas: Ackley(n); Beale(2); Easom(2); FourPeaks(2);
  !                                Griewank(n); Rastrigin(n); Rosenbrock(n);
  !                                Schwefel(n); Shubert(2); Sum Squares(n).
  !
  ! Ate a data da ultima modificacao, as funcoes acima podiam ser obtidas em:
  ! http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/
  ! TestGO_files/Page364.htm
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE limites (nomeFuncao,linf,lsup)

    !Nao permite declaracao implicita
    IMPLICIT NONE

    ! Variaveis de entrada da subrotina
    DOUBLE PRECISION, INTENT(OUT) :: linf, lsup
    CHARACTER, INTENT(IN) :: nomeFuncao*16

    ! Selecao e passagem dos valores de limites conforme "nomeFuncao" definido
    SELECT CASE (nomeFuncao)

      ! Limites para a funcao de Ackley (n-dimensional)
      CASE ("ackley")
        linf = -32.768D0
        lsup =  32.768D0

      ! Limites para a funcao de Beale (2-dimensional)
      CASE ("beale")
        linf = -4.5D0
        lsup =  4.5D0

      ! Limites para a funcao de Easom (2-dimensional)
      CASE ("easom")
        linf = -100.0D0
        lsup =  100.0D0

      ! Limites para a funcao Four Peaks (2-dimensional)
      CASE ("fourpeaks")
        linf = -5.0D0
        lsup =  5.0D0

      ! Limites para a funcao de Griewank (n-dimensional)
      CASE ("griewank")
        linf = -600.0D0
        lsup =  600.0D0

      ! Limites para a funcao de Rastrigin (n-dimensional)
      CASE ("rastrigin")
        linf = -5.12D0
        lsup =  5.12D0

      ! Limites para a funcao de Rosenbrock (n-dimensional)
      CASE ("rosenbrock")
        linf = -2.048D0
        lsup =  2.048D0

      ! Limites para a funcao de Schwefel (n-dimensional)
      CASE ("schwefel")
        linf = -500.0D0
        lsup =  500.0D0

      ! Limites para a funcao de Shubert (2-dimensional)
      CASE ("shubert")
        linf =  -10.0D0
        lsup =   10.0D0

      ! Limites para a funcao de Shekel (2-dimensional)
      CASE ("shekel")
        linf =  -65.536
        lsup =   65.536

      ! Limites padrões (funcao Sum Squares) (n-dimensional)
      CASE DEFAULT
        linf = -10.0D0
        lsup =  10.0D0
        WRITE(*,*) "Valores padroes adotados..."

    ! Fim da atribuicao de limites conforme "nomeFuncao"
    END SELECT

  ! Fim da funcao para definicao dos limites
  END SUBROUTINE


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Funcao para calcular e retornar o valor da funcao objetivo
  ! 
  ! Entrada: nome da funcao, em character*16 (nomeFuncao)
  !          dimensoes do problema, em integer (d)
  !          vetor de d-posicoes ou solucao candidata, em double precision (x)
  ! 
  ! Retorno: valor da funcao objetivo, em double precision
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION FUNCTION calcFuncao(nomeFuncao,d,x)

    !Nao permite declaracao implicita
    IMPLICIT NONE

    ! Variaveis de entrada da funcao
    CHARACTER, INTENT(IN) :: nomeFuncao*16
    INTEGER, INTENT(IN) :: d
    DOUBLE PRECISION, DIMENSION (d), INTENT(IN) :: x

    ! Variaveis locais:
    !                   i -> contador de loops
    !                   s1 -> somatorio 1
    !                   p1 -> produtorio 1
    !                   d1 -> divisao 1
    !                   ai, bi -> elementos adicionais das equacoes
    INTEGER :: i
    DOUBLE PRECISION :: s1, s2, p1, d1, ai, bi

    ! Constantes matematicas necessarias: PI
    DOUBLE PRECISION :: pi
    PARAMETER (pi = 3.1415926535897932384626433832795D0)

    ! Seleciona qual funcao sera calculada conforme o valor de "nomeFuncao"
    SELECT CASE (nomeFuncao)

      ! Calculo para a funcao de Ackley (n-dimensoes)
      ! -32.768 <= x(i) <= 32.786
      ! f(x*) = 0
      ! x*(i) = 0, i=1:n
      CASE ("ackley")
        s1 = 0.0D0
        s2 = 0.0D0
        DO i = 1, d
          s1 = s1 + x(i)**2
          s2 = s2 + DCOS(2.0D0*pi*x(i))
        ENDDO
        calcFuncao = -20.0D0*DEXP(-0.2D0*SQRT(1.0D0/d*s1))-DEXP(1.0D0/d*s2)+ &
                      20.0D0+DEXP(1.0D0)

      ! Calculo para a funcao de Beale (2-dimensoes)
      ! -4.5 <= x(i) <= 4.5
      ! f(x*) = 0
      ! x* = 3, 0.5
      CASE ("beale")
        calcFuncao = (1.5D0-x(1) * (1.0D0-x(2)))**2 + (2.25D0-x(1) * & 
                     (1.0D0-x(2)**2))**2 + (2.625D0-x(1) * (1.0D0-x(2)**3))**2

      ! Calculo para a funcao de Easom (2-dimensoes)
      ! -100 <= x(i) <= 100
      ! f(x*) = -1
      ! x*(i) = pi, i = 1,2
      CASE ("easom")
        calcFuncao = -DCOS(x(1))*DCOS(x(2))*DEXP(-(x(1)-pi)**2-(x(2)-pi)**2)

      ! Calculo para a funcao Four Peaks (2-dimensoes)
      ! -5 <= x(i) <= 5
      ! f(x*) = -2.0
      ! 2 minimos globais: (0,0) e (0,-4)
      CASE ("fourpeaks")
        calcFuncao = -(DEXP(-(x(1)-4.0D0)**2-(x(2)-4.0D0)**2)+ &
                     DEXP(-(x(1)+4.D0)**2-(x(2)-4.0D0)**2)+ &
                     2.0D0*DEXP(-x(1)**2-(x(2)+4.0D0)**2)+ &
                     2.0D0*DEXP(-x(1)**2-x(2)**2))

      ! Calculo para a funcao de Griewank (n-dimensoes)
      ! -600 <= x(i) <= 600
      ! f(x*) = 0
      ! x*(i) = 0, i=1:n
      CASE ("griewank")
        s1 = 0.0D0
        p1 = 1.0D0
        DO i = 1, d
          s1 = s1 + x(i)**2
          p1 = p1 * DCOS(x(i)/SQRT(REAL(i)))
        ENDDO
        calcFuncao = s1/4000.0D0 - p1 + 1.D0

      ! Calculo para a funcao de Rastrigin (n-dimensional)
      ! -5.12 <= x(i) <= 5.12
      ! f(x*) = 0
      ! x*(i) = 0, i=1:n
      CASE ("rastrigin")
        s1 = 0.0D0
        DO i = 1, d
          s1 = s1 + x(i)**2-10.0D0*DCOS(2.0D0*pi*x(i))
        END DO
        calcFuncao = 10.0D0 * d + s1

      ! Calculo para a funcao de Rosenbrock (n-dimensional)
      ! -2.048 <= x(i) <= 2.048
      ! f(x*) = 0
      ! x*(i) = 1, i=1:n
      CASE ("rosenbrock")
        s1 = 0.0D0
        DO i = 1, d-1
          s1 = s1 + 100.0D0*(x(i)**2 - x(i+1))**2 + (x(i) - 1.0D0)**2
        END DO
        calcFuncao = s1

      ! Calculo para a funcao de Schwefel (n-dimensional)
      ! -500 <= x(i) <= 500
      ! f(x*) = 0
      ! x*(i) = 420.9687, i=1:n
      CASE ("schwefel")
        s1 = 0.0D0
        DO i = 1, d
          s1 = s1 + (-x(i) * DSIN(SQRT(ABS(x(i)))))
        ENDDO
        calcFuncao = 418.9829D0 * d + s1

      ! Calculo para a funcao de Shubert (2-dimensional)
      ! -10 <= x(i) <= 10
      ! f(x*) = -186.7309
      ! 18 minimos globais
      CASE ("shubert")
        s1 = 0.0D0
        s2 = 0.0D0
        DO i = 1, 5
          s1 = s1 + i * DCOS((i+1)*x(1)+i)
          s2 = s2 + i * DCOS((i+1)*x(2)+i)
        ENDDO
        calcFuncao = s1 * s2

      ! Calculo para a funcao de Shekel (2-dimensional)
      ! -65.536 <= xi <= 65.536
      ! f(x*) = 499.002
      ! x* = (-32, -32)
      CASE ("shekel")
        s1 = 0.0D0
        DO i = 0, 24
          ai = 16.0D0 * ((MOD(i,5)) - 2.0D0)
          bi = 16.0D0 * ((i/5.0D0) - 2.0D0)
          s1 = s1+(1.0D0/(1.0D0+i+(x(1)-ai)**6+(x(2)-bi)**6))
        END DO
        d1 = 1.0D0 / (0.002D0 + s1)
        calcFuncao = -(500.0D0 - d1)

      ! Calculo da funcao "Sum Squares", do somatorio de "x" (n-dimensional)
      ! -10 <= x(i) <= 10
      ! f(x*) = 0
      ! x*(i) = 0, i=1:n
      CASE DEFAULT
        calcFuncao = 0.0D0
        DO i=1,d
          calcFuncao = calcFuncao + i*x(i)**2
        END DO

    ! Fim do SELECT para a funcao definida em "nomeFuncao"
    END SELECT

    ! Retorna o valor calculado (calcFuncao) pelo SELECT acima
    RETURN

  ! Fim da funcao para calcular o valor da funcao objetivo
  END FUNCTION

!Fim do modulo
END MODULE funcoes
