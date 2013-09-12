!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MULTIPLE PARTICLE COLLISION ALGORITHM (MPCA)
! 
! Desenvolvido por: Eduardo Favero Pacheco da Luz (CAP/INPE)
! Baseado no PCA por Wagner F. Sacco
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Inicio do programa
PROGRAM MPCA

  !Usa o modulo de funcoes
  USE nr
  USE nrtype
  USE ran_state, ONLY : ran_seed
  USE funcoes
  USE estat
  USE ModuloRNA

  !Evita uso implicito
  IMPLICIT NONE

  !Inclui Bibliotecas MPI
  INCLUDE 'mpif.h'

  !Variaveis de amplo aspecto e definicao de parametros
    ! l1a,b -> contador de loop generico
    ! lexp -> contador de loop de experimentos computacionais
    ! idum -> semente geradora específica
    ! i -> contador de loop generico
    ! j -> contador de loop generico
    ! vetseed -> vetor de sementes geradoras de numeros aleatorios
    ! dimensoes -> dimensoes do problema a ser resolvido
    ! iteracoes -> quantidade de iteracoes do PCA (modo global)
    ! iterPert -> quantidade de iteracoes do PCA interno (modo local)
    !             na maioria das vezes iterPert = iteracoes
    ! max_nafo -> numero maximo de avaliacoes da funcao objetivo
    ! nafo -> numero atual de avaliacoes da funcao objetivo
    ! ciclo -> ciclo de iteracoes para efetuar comunicacao via Blackboard
    ! tot_nafo -> total de AFO de todos os processadores
    ! p_scat -> probabilidade de espalhamento
    ! espaco_inf -> limite inferior do espaco de buscas (fixo)
    ! espaco_sup -> limite superior do espaco de buscas (fixo)
    ! alea -> numero aleatorio sorteado
    ! superior -> limite inferior (variavel)
    ! inferior -> limite superior (variavel)
    ! nomeF -> nome da funcao a ser otimizada
    ! texp -> numero de experimentos computacionais (semente geradora NR)
    ! tipo_p -> tipo de probabilidade usada na funcao Scattering
    ! vMed -> vetor para calcular a media dos resultados dos experimentos
    ! vDP -> vetor para calcular o desvio padrao dos resultados dos experimentos
    ! vFO -> vetor para armazenar o valor da funcao objetivo (para media)
    ! vNafo -> vetor para armazenar o NAFO (para media)
    ! vExp -> matriz para armazenar os resultados obtidos (para media)
    ! nPartNode -> numero de particulas por no
    ! lo_small -> limite inferior da perturbacao na busca local
    ! up_small -> limite superior da perturbacao na busca local
    ! paramRNA -> numero de parametros para a RNA (arquivo entrada.txt)
  INTEGER :: dimensoes, idum, texp, l1a, l1b, lexp, i, j, tipo_p, ciclo
  INTEGER :: nPartNode, paramRNA
  INTEGER*8 :: iteracoes, iterPert, max_nafo, nafo, tot_nafo
  INTEGER, ALLOCATABLE, DIMENSION(:) :: vetseed
  INTEGER*8, ALLOCATABLE, DIMENSION(:) :: vNafo
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: vMed, vDP, vFO
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: vExp
  REAL :: etime, total, elapsed(2)
  REAL*8 :: p_scat,superior,inferior
  REAL*8, ALLOCATABLE, DIMENSION(:) :: espaco_inf, espaco_sup
  REAL*8 :: lo_small, up_small
  CHARACTER*16 :: nomeF
  CHARACTER*30 :: caminho, path_ent, path_limEnt, path_result
  CHARACTER*20 :: arq_ent, arq_limEnt, arq_result

  REAL(SP), ALLOCATABLE, DIMENSION(:) :: harvest
!  PARAMETROS PASSADOS PELO EDUARDO
!  PARAMETER (iteracoes=100, iterPert=10, max_nafo=1000, dimensoes=7)
!  PARAMETER (texp = 10, nomeF = "griewank", nPartNode = 1)
!  PARAMETER (tipo_p = 1, ciclo = 10, lo_small = 0.8D0, up_small = 1.2D0)
  PARAMETER (iteracoes=20, iterPert=10, max_nafo=500, dimensoes=7)
  PARAMETER (texp = 10, nomeF = "griewank", nPartNode = 1)
  PARAMETER (tipo_p = 1, ciclo = 10, lo_small = 0.8D0, up_small = 1.2D0)

  PARAMETER (paramRNA = 16)

  !Variáveis para a integração com a RNA
  DOUBLE PRECISION, DIMENSION(paramRNA,1) :: parametros

  !Variáveis para a paralelizacao
    ! np -> numero do processador atual
    ! tp -> total de processadores envolvidos
    ! nump_np -> numero de particulas para cada processador
    ! total_p -> numero total de particulas que serao usadas em paralelo
    ! p_ini -> numero da particula inicial para cada processador
    ! p_fin -> numero da particula final para cada processador
    ! melhor_p -> numero do processador que obteve o melhor resultado
  INTEGER np, tp, ierr, status(MPI_STATUS_SIZE), t_part, p_ini, p_fin, part_np
  INTEGER melhor_p

  !Definicao da estrutura da solucao
  TYPE tipo_solucao
    REAL*8 :: custo
    REAL*8 :: solucao(dimensoes)
  END TYPE tipo_solucao

  !Construcao da solucao antiga (Old_Config)
  TYPE (tipo_solucao) Old_Config

  !Construcao da nova solucao (New_Config)
  TYPE (tipo_solucao) New_Config

  !Construcao da melhor solucao (Best_Config)
  TYPE (tipo_solucao) Best_Config

  !Inicializa o MPI
  CALL MPI_Init(ierr)

  !Determina o total de processadores em uso
  CALL MPI_Comm_Size(MPI_COMM_WORLD, tp, ierr)

  !Determina o numero deste processador
  CALL MPI_Comm_Rank(MPI_COMM_WORLD, np, ierr)

  !Construindo caminhos para os arquivos
  caminho = './dados/'
  arq_result = 'result.out'
  arq_limEnt = 'limEntradas.txt'
  arq_ent =	'entrada.txt'
  path_result = trim(caminho)//trim(arq_result)
  path_limEnt = trim(caminho)//trim(arq_limEnt)
  path_ent = trim(caminho)//trim(arq_ent)

  !Arquivo para armazenar evolucao da solucao e tempo de execucao
  IF (np .EQ. 0) THEN
    OPEN(UNIT=10, FILE = path_result, STATUS="replace")
    !OPEN(UNIT=20, FILE="decai.out", STATUS="replace")
  END IF

  !Definição do número de partículas por processador
  t_part = tp * nPartNode !Total de particulas que serao usadas (por enquanto uma por processador)
  part_np = t_part/tp !Calcula quantas particulas serao por processador
  p_ini = (np * part_np) + 1 !Calcula o numero da particula inicial do processador
  p_fin = p_ini + part_np - 1 !Calcula o numero da particula final do processador

  !Define os limites do espaco de buscas
  !CALL limites(nomeF,espaco_inf,espaco_sup)
  ALLOCATE(espaco_inf(dimensoes))
  ALLOCATE(espaco_sup(dimensoes))

  !Le os limites do arquivo limEntradas.txt
  OPEN(UNIT=40, FILE = path_limEnt, STATUS="old")
  DO i = 1, dimensoes
    READ(40,*) espaco_inf(i), espaco_sup(i)
  END DO
  CLOSE(40)

  !Le o restante das informações de entrada
  OPEN(UNIT=30, FILE = path_ent, STATUS="old")
  DO i = 1, paramRNA
    READ(30,*) parametros(i,1)
  END DO
  CLOSE(30)

  !Vetor de sementes geradoras de numeros aleatorios
  ALLOCATE(vetseed(25))
  vetseed = (/ 876, 354, 912, 87, 7234, 9299, 12, 9911, 3244, 10, &
               222, 88, 1111, 8876, 4334, 5138, 76, 1934, 45777, &
               5521, 231, 45552, 800, 5, 3332 /)

  !Alocacao de variaveis dinamicas
  ALLOCATE(vMed(dimensoes))
  ALLOCATE(vDP(dimensoes))
  ALLOCATE(vFO(texp))
  ALLOCATE(vNafo(texp))
  ALLOCATE(vExp(texp,dimensoes))
  ALLOCATE(harvest(dimensoes))

  !Loop para experimentos numericos
  DO lexp = 1, texp

    !Total de AFO
    tot_nafo = 0

    !Fixa a semente geradora de numeros aleatorios
    idum = vetseed(lexp) + ((tp*3)*np)
    CALL ran_seed(sequence=idum)

    !Gera uma solucao inicial para Old_Config
    CALL ran3(harvest)
    DO l1a = 1, dimensoes
      Old_Config%solucao(l1a) = (harvest(l1a)*(espaco_sup(l1a)-espaco_inf(l1a)))+espaco_inf(l1a)
      ! Monta a lista de parametros a ser avaliada
      parametros(l1a,1) = Old_Config%solucao(l1a)
    END DO
    ! Converte para inteiro os parametros da RNA
    Old_Config%solucao(1) = NINT(Old_Config%solucao(1))
    Old_Config%solucao(2) = NINT(Old_Config%solucao(2))
    Old_Config%solucao(3) = NINT(Old_Config%solucao(3))
    Old_Config%solucao(4) = NINT(Old_Config%solucao(4))
    Old_Config%solucao(5) = NINT(Old_Config%solucao(5))
    parametros(1,1) = NINT(parametros(1,1)) ! Numero de camadas escondidas
    parametros(2,1) = NINT(parametros(2,1)) ! Numero de neuronios na 1a camada
    parametros(3,1) = NINT(parametros(3,1)) ! Numero de neuronios na 2a camada
    parametros(4,1) = NINT(parametros(4,1)) ! Numero de neuronios na 3a camada
    parametros(5,1) = NINT(parametros(5,1)) ! Tipo da funcao de ativacao

    !Old_Config%custo = calcFuncao(nomeF,dimensoes,Old_Config%solucao)
    Old_Config%custo = Rede_Neural_BP (parametros)
    IF (np .EQ. 0) THEN
      write(*,*) "Experimento ", lexp,"-->",Old_Config%custo
    ENDIF
    !Uma AFO feita
    nafo = 1

    !No momento, a solucao inicial eh a melhor solucao
    Best_Config = Old_Config

    !Inicializa contador do loop
    l1a = 0

    !Loop principal
      !Criterios de parada:
        ! (l1a .LT. iteracoes) -> numero maximo de iteracoes do PCA (canonico)
        ! (tot_nafo .LT. max_nafo) -> numero maximo de avaliacoes da funcao objetivo
        ! (Best_Config%custo .GT. 1.0D-8) -> precisao da funcao objetivo
    DO WHILE (tot_nafo .LT. max_nafo)

      !Incrementa controle de laco
      l1a = l1a + 1

      !Retorno visual
      IF (np .EQ. 0) THEN
        WRITE(*,*) "Iteracao ", l1a, "-> NAFO = ", tot_nafo
      END IF

      !Chama e atualiza o Blackboard quando liberado pelo ciclo
      IF (MOD(l1a,ciclo) .EQ. 0) THEN
        CALL blackboard(tp,np,dimensoes,Best_Config,nafo,tot_nafo)
      END IF

      !Arquivo para armazenar evolucao da solucao
      !IF (np .EQ. 0) THEN
      !  WRITE(10,*) l1a, Best_Config%custo
      !END IF

      !Loop das particulas (para controle das particulas no MPCA)
      DO l1b = p_ini, p_fin

        !Chama a rotina de perturbacao da particula
        CALL Perturbation(nomeF,dimensoes,espaco_inf,espaco_sup,&
          Old_Config%solucao,Old_Config%custo,New_Config%solucao,&
          New_Config%custo,nafo,parametros)

        IF (New_Config%custo .LT. Old_Config%custo) THEN
          IF (New_Config%custo .LT. Best_Config%custo) THEN
            Best_Config = New_Config
          END IF
          Old_Config = New_Config
          !Chama a rotina de exploracao da particula
          CALL Exploration(nomeF,iterPert,dimensoes,espaco_inf,espaco_sup,&
               Old_Config%solucao,Old_Config%custo,New_Config%solucao,&
               New_Config%custo,Best_Config%solucao,Best_Config%custo,&
               nafo,lo_small,up_small,parametros)
        ELSE
          !Chama a rotina de espalhamento da particula
          CALL Scattering(nomeF,iterPert,dimensoes,espaco_inf,espaco_sup,&
               Old_Config%solucao,Old_Config%custo,New_Config%solucao,&
               New_Config%custo,Best_Config%solucao,Best_Config%custo,nafo,&
               tipo_p,lo_small,up_small,parametros)
        END IF
      !Fim do loop de controle das particulas do MPCA
      END DO

      !Plota decaimento da FO
      !IF (np .EQ. 0) THEN
      !  WRITE(20,*) l1a, Best_Config%custo
      !END IF

    !Fim do loop princpial
    END DO
    
    !Atualizacao final do Blackboard
    CALL blackboard(tp,np,dimensoes,Best_Config,nafo,tot_nafo)

    !Arquivo para armazenar evolucao da solucao
    !IF (np .EQ. 0) THEN
    !  WRITE(10,*) Best_Config%custo
    !END IF

    !Imprime a melhor solucao
    IF (np .EQ. 0) THEN
      WRITE(10,*) ''
      WRITE(10,*) 'Melhor solucao obtida:'
      WRITE(10,*) Best_Config%solucao
      WRITE(10,*) ''
      WRITE(10,*) 'Melhor fitness:'
      WRITE(10,*) Best_Config%custo
      WRITE(10,*) ''
      WRITE(10,*) 'Obtida no processador:'
      WRITE(10,*) melhor_p
      WRITE(10,*) ''
    END IF

    !Armazena os resultados para calcular a media e desvio padrao
    IF (np .EQ. 0) THEN
      WRITE(10,*) '-------------------------'
      WRITE(10,*) Best_Config%custo
      CALL FLUSH(10)
      vExp(lexp,:) = Best_Config%solucao
      vFO(lexp) = Best_Config%custo
      vNafo(lexp) = tot_nafo
    END IF

  !Fim do loop de experimentos computacionais
  END DO

  IF (np .EQ. 0) THEN

    ! Calcula e imprime a media e desvio padrao: Posicao
    Vmed = 0.0D0
    Vdp = 0.0D0
    CALL mediaDP (texp,dimensoes,vExp,vMed,vDP)
    WRITE(10,*) '======================='
    WRITE(10,*) 'Resultado Medio: ',vMed
    WRITE(10,*) 'Desvio Padrao:   ',vDP
    WRITE(10,*) ''

    ! Calcula e imprime a media e desvio padrao: valor da FO
    Vmed = 0.0D0
    Vdp = 0.0D0
    CALL mediaDP (texp,1,vFO,vMed,vDP)
    WRITE(10,*) '======================='
    WRITE(10,*) 'Valor FO Medio: ',vMed(1)
    WRITE(10,*) 'Desvio Padrao:  ',vDP(1)
    WRITE(10,*) ''

    ! Calcula e imprime a media e desvio padrao: NAFO
    Vmed = 0.0D0
    Vdp = 0.0D0
    CALL mediaDP (texp,1,DBLE(vNafo),vMed,vDP)
    WRITE(10,*) '======================='
    WRITE(10,*) 'Total NAFO:    ',INT(SUM(vNafo))
    WRITE(10,*) 'NAFO Medio:    ',vMed(1)
    WRITE(10,*) 'NAFO Medio por processador:    ',vMed(1)/tp
    WRITE(10,*) 'Desvio Padrao: ',vDP(1)
    WRITE(10,*) ''

    ! Tempo de execucao do programa
    total = etime(elapsed)
    WRITE(10,*) '======================='
    WRITE(10,*) 'Tempo de execucao do programa:'
    WRITE(10,*) 'Total   = ', total, 'segundos'
    WRITE(10,*) 'Usuario = ', elapsed(1), 'segundos'
    WRITE(10,*) 'Sistema = ', elapsed(2), 'segundos'
    WRITE(10,*) ''

   ! Fecha arquivos
   CLOSE(10)
   !CLOSE(20)

  END IF

  !Finaliza o MPI
  CALL MPI_Finalize(ierr)

  !Desaloca variaveis dinamicas
  DEALLOCATE(vetseed,vMed,vDP,vFO,vNafo,vExp)

  !Aqui vem as subrotinas
  CONTAINS

  !///////////////////////////////////////////////////////////////////////////////
  SUBROUTINE Perturbation(nomeF,dimensoes,espaco_inf,espaco_sup,Old_Config,&
             Old_Fitness,New_Config,New_Fitness,nafo,parametros)
  !///////////////////////////////////////////////////////////////////////////////
    USE nrtype
    USE funcoes
    IMPLICIT NONE
    INTEGER :: i, dimensoes
    INTEGER*8 :: nafo
    REAL*8 :: superior, inferior, espaco_inf(:), espaco_sup(:), Old_Config(:)
    REAL*8 :: Old_Fitness, New_Config(:), New_Fitness
    DOUBLE PRECISION :: parametros (:,:)
    REAL(SP) :: alea
    CHARACTER*16 :: nomeF

    DO i = 1, dimensoes
      !2 IF's para manter o valor se o numero de camadas nao for 3
      IF (i .EQ. 3) THEN
        IF (New_Config(1) .EQ. 1) THEN
          !Ao avaliar o 3o parametro, verifica se o numero de camadas eh 1
          !Se for, entao nao muda o valor
          New_Config(i) = Old_Config(i)
          CYCLE
        END IF
      END IF
      IF (i .EQ. 4) THEN
        IF ( (New_Config(1) .EQ. 1) .OR. (New_Config(1) .EQ. 2)) THEN
          !Ao avaliar o 4o parametro, verifica se o numero de camadas eh 1 ou 2
          !Se for, entao nao muda o valor
          New_Config(i) = Old_Config(i)
          CYCLE
        END IF
      END IF
      superior = espaco_sup(i)
      inferior = espaco_inf(i)
      CALL ran3(alea)
      New_Config(i) = Old_Config(i) + ((superior - Old_Config(i))*alea) - &
                      ((Old_Config(i) - inferior)*(1.0D0-alea))
      IF (New_Config(i) .GT. superior) THEN
        New_Config(i) = superior
      ELSE
        IF (New_Config(i) .LT. inferior) THEN
          New_Config(i) = inferior
        END IF
      END IF
      parametros(i,1) =  New_Config(i)
    END DO

    ! Converte para inteiro os parametros da RNA
    New_Config(1) = NINT(New_Config(1))
    New_Config(2) = NINT(New_Config(2))
    New_Config(3) = NINT(New_Config(3))
    New_Config(4) = NINT(New_Config(4))
    New_Config(5) = NINT(New_Config(5))
    parametros(1,1) = NINT(parametros(1,1)) ! Numero de camadas escondidas
    parametros(2,1) = NINT(parametros(2,1)) ! Numero de neuronios na 1a camada
    parametros(3,1) = NINT(parametros(3,1)) ! Numero de neuronios na 2a camada
    parametros(4,1) = NINT(parametros(4,1)) ! Numero de neuronios na 3a camada
    parametros(5,1) = NINT(parametros(5,1)) ! Tipo da funcao de ativacao

    !New_Fitness = calcFuncao(nomeF,dimensoes,New_Config)
    New_Fitness = Rede_Neural_BP (parametros)
    nafo = nafo + 1

  END SUBROUTINE ! Fim da subrotina

  !///////////////////////////////////////////////////////////////////////////////
  SUBROUTINE Exploration(nomeF,iteracoes,dimensoes,espaco_inf,espaco_sup,&
             Old_Config,Old_Fitness,New_Config,New_Fitness,Best_Config,&
             Best_Fitness,nafo,lo_small,up_small,parametros)
  !///////////////////////////////////////////////////////////////////////////////
    IMPLICIT NONE
    INTEGER :: l2, dimensoes
    INTEGER*8 :: nafo, iteracoes
    REAL*8 :: superior, inferior, espaco_inf(:), espaco_sup(:), lo_small, up_small
    REAL*8 :: Old_Config(:),Old_Fitness,New_Config(:)
    REAL*8 :: Best_Config(:), Best_Fitness, New_Fitness
    DOUBLE PRECISION :: parametros(:,:)
    CHARACTER*16 :: nomeF

    DO l2 = 1, iteracoes
      CALL Small_Perturbation (nomeF, dimensoes, espaco_inf, espaco_sup, &
           Old_Config, New_Config, New_Fitness, nafo, lo_small, up_small, &
           parametros)
      IF (New_Fitness .LT. Old_Fitness) THEN
        IF (New_Fitness .LT. Best_Fitness) THEN
          Best_Config = New_Config
          Best_Fitness = New_Fitness
        END IF
        Old_Config = New_Config
        Old_Fitness = New_Fitness
      END IF
    END DO

  END SUBROUTINE ! Fim da subrotina

  !///////////////////////////////////////////////////////////////////////////////
  SUBROUTINE Small_Perturbation(nomeF,dimensoes,espaco_inf,espaco_sup,&
             Old_Config,New_Config,New_Fitness,nafo,lo_small,up_small,&
             parametros)
  !///////////////////////////////////////////////////////////////////////////////
    USE nrtype
    USE funcoes
    IMPLICIT NONE
    INTEGER :: dimensoes, l3
    INTEGER*8 :: nafo
    REAL*8 :: espaco_inf(:), espaco_sup(:), inferior, superior, Old_Config(:)
    REAL*8 :: New_Config(:), New_Fitness, alea, lo_small, up_small
    DOUBLE PRECISION :: parametros(:,:)
    CHARACTER*16 :: nomeF
    REAL(SP), DIMENSION(3*dimensoes) :: harvest

    CALL ran3(harvest)
    DO l3 = 1, dimensoes
      !2 IF's para manter o valor se o numero de camadas nao for 3
      IF (l3 .EQ. 3) THEN
        IF (New_Config(1) .EQ. 1) THEN
          !Ao avaliar o 3o parametro, verifica se o numero de camadas eh 1
          !Se for, entao nao muda o valor
          New_Config(l3) = Old_Config(l3)
          CYCLE
        END IF
      END IF
      IF (l3 .EQ. 4) THEN
        IF ( (New_Config(1) .EQ. 1) .OR. (New_Config(1) .EQ. 2) ) THEN
          !Ao avaliar o 4o parametro, verifica se o numero de camadas eh 1 ou 2
          !Se for, entao nao muda o valor
          New_Config(l3) = Old_Config(l3)
          CYCLE
        END IF
      END IF
      superior = ((harvest(l3)*(up_small-1.0D0))+1.0D0) * Old_Config(l3)
      IF (superior .GT. espaco_sup(l3)) THEN
        superior = espaco_sup(l3)
      END IF
      inferior = ((harvest(l3+dimensoes)*(1.0D0-lo_small))+lo_small) * Old_Config(l3)
      IF (inferior .LT. espaco_inf(l3)) THEN
        inferior = espaco_inf(l3)
      END IF
      alea = harvest(l3+2*dimensoes)
      New_Config(l3) = Old_Config(l3) + ((superior - Old_Config(l3))*alea) -&
                       ((Old_Config(l3) - inferior)*(1.0D0-alea))
      parametros(l3,1) = New_Config(l3)
    END DO

    ! Converte para inteiro os parametros da RNA
    New_Config(1) = NINT(New_Config(1))
    New_Config(2) = NINT(New_Config(2))
    New_Config(3) = NINT(New_Config(3))
    New_Config(4) = NINT(New_Config(4))
    New_Config(5) = NINT(New_Config(5))
    parametros(1,1) = NINT(parametros(1,1)) ! Numero de camadas escondidas
    parametros(2,1) = NINT(parametros(2,1)) ! Numero de neuronios na 1a camada
    parametros(3,1) = NINT(parametros(3,1)) ! Numero de neuronios na 2a camada
    parametros(4,1) = NINT(parametros(4,1)) ! Numero de neuronios na 3a camada
    parametros(5,1) = NINT(parametros(5,1)) ! Tipo da funcao de ativacao

    !New_Fitness = calcFuncao(nomeF,dimensoes,New_Config)
    New_Fitness = Rede_Neural_BP (parametros)
    nafo = nafo + 1

  END SUBROUTINE ! Fim da subrotina

  !///////////////////////////////////////////////////////////////////////////////
  SUBROUTINE Scattering(nomeF,iteracoes,dimensoes,espaco_inf,espaco_sup,&
             Old_Config,Old_Fitness,New_Config,New_Fitness,Best_Config,&
             Best_Fitness,nafo,tipo_p,lo_small,up_small,parametros)
  !///////////////////////////////////////////////////////////////////////////////
    USE nrtype
    USE funcoes
    IMPLICIT NONE
    INTEGER :: l2, dimensoes, tipo_p
    INTEGER*8 :: nafo, iteracoes
    REAL*8 :: espaco_inf(:), espaco_sup(:), Old_Config(:), Old_Fitness
    REAL*8 :: New_Config(:), New_Fitness, Best_Config(:)
    REAL*8 :: Best_Fitness, p_scat, lo_small, up_small
    DOUBLE PRECISION :: parametros(:,:)
    REAL(SP) :: alea
    REAL(SP), DIMENSION(dimensoes) :: harvest
    CHARACTER*16 :: nomeF

    SELECT CASE (tipo_p)
      CASE (1)
        !Probabilidade associada a uma função exponencial truncada
        p_scat = 1.0D0-(Best_Fitness/New_Fitness)
      CASE (2)
        !Probabilidade associada à FDP de Cauchy
        p_scat=1.0D0-(1.0D0/(pi*1.0D0*(1.0D0+((New_Fitness-Best_Fitness) / &
               1.0D0)**2)))
      CASE (3)
        !Probabilidade associada à FDP de Cauchy (complemento)
        p_scat=1.0D0/(pi*1.0D0*(1.0D0+((New_Fitness-Best_Fitness)/1.0D0)**2))
    END SELECT

    CALL ran3(alea)
    IF (p_scat .GT. alea) THEN
      CALL ran3(harvest)
      DO l2 = 1, dimensoes
        !2 IF's para manter o valor se o numero de camadas nao for 3
        IF (l2 .EQ. 3) THEN
          IF (Old_Config(1) .EQ. 1) THEN
            !Ao avaliar o 3o parametro, verifica se o numero de camadas eh 1
            !Se for, entao nao muda o valor
            CYCLE
          END IF
        END IF
        IF (l2 .EQ. 4) THEN
          IF ( (Old_Config(1) .EQ. 1) .OR. (Old_Config(1) .EQ. 2) ) THEN
            !Ao avaliar o 4o parametro, verifica se o numero de camadas eh 1 ou 2
            !Se for, entao nao muda o valor
            CYCLE
          END IF
        END IF
        !Atribui novos valores aleatorios
        Old_Config(l2) = (harvest(l2)*(espaco_sup(l2)-espaco_inf(l2)))+espaco_inf(l2)
      END DO
      Old_Fitness = calcFuncao(nomeF,dimensoes,Old_Config)
      nafo = nafo + 1
    ELSE
      CALL Exploration(nomeF,iteracoes, dimensoes, espaco_inf, espaco_sup, &
           Old_Config,Old_Fitness, New_Config, New_Fitness, Best_Config, &
           Best_Fitness, nafo, lo_small, up_small, parametros)
    END IF

  END SUBROUTINE! Fim da subrotina

  !///////////////////////////////////////////////////////////////////////////////
  SUBROUTINE blackboard(tp,np,dimensoes,Best_Config,nafo,tot_nafo)
  !///////////////////////////////////////////////////////////////////////////////
    INCLUDE 'mpif.h'
    INTEGER :: i,ierr,status(MPI_STATUS_SIZE)
    INTEGER*8, INTENT(INOUT) :: nafo, tot_nafo
    INTEGER, INTENT(IN) :: tp, np, dimensoes
    REAL*8, DIMENSION(dimensoes+2) :: Send_Config
    TYPE (tipo_solucao), INTENT(INOUT) :: Best_Config

    !Atualiza a melhor solucao
    IF (np .EQ. 0) THEN ! Se for o processador zero
      !Soma o NAFO do Zero e reinicia seu valor
      tot_nafo = tot_nafo + nafo
      nafo = 0
      !Recebe os Config's dos outros processadores
      DO i=1, tp-1
        !Recebe
        CALL MPI_Recv(Send_Config, dimensoes+2, MPI_DOUBLE_PRECISION, &
             MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, status, ierr)
        !Soma NAFO
        tot_nafo = tot_nafo + INT(Send_Config(dimensoes+2))
        !Verifica se o recebido é melhor
        IF (Send_Config(dimensoes+1) .LT. Best_Config%custo) THEN
          DO j=1, dimensoes
            Best_Config%solucao(j) = Send_Config(j)
          ENDDO
          Best_Config%custo = Send_Config(dimensoes+1)
        ENDIF
      ENDDO
      !Depois de termos o melhor de todos monta o Send_Config
      DO i=1,dimensoes
        Send_Config(i) = Best_Config%solucao(i)
      ENDDO
      Send_Config(dimensoes+1) = Best_Config%custo
      Send_Config(dimensoes+2) = DBLE(tot_nafo)
      !Envia para todos os processadores
!     DO i=1, tp-1
!       CALL MPI_Send(Send_Config, dimensoes+2, MPI_DOUBLE_PRECISION, i, 0, &
!            MPI_COMM_WORLD, ierr)
!     ENDDO
      CALL MPI_Bcast(Send_Config, dimensoes+2, MPI_DOUBLE_PRECISION, 0, &
                     MPI_COMM_WORLD, ierr)

    ELSE ! Os outros processadores
      !Monta o Send_Config
      DO i=1,dimensoes
        Send_Config(i) = Best_Config%solucao(i)
      ENDDO
      Send_Config(dimensoes+1) = Best_Config%custo
      Send_Config(dimensoes+2) = DBLE(nafo)
      !Envia para o processador Zero
      CALL MPI_Send(Send_Config, dimensoes+2, MPI_DOUBLE_PRECISION, 0, 0, &
           MPI_COMM_WORLD, ierr)
      !Depois do processador zero identificar o melhor, este retorna
!     CALL MPI_Recv(Send_Config, dimensoes+2, MPI_DOUBLE_PRECISION, 0, 0, &
!          MPI_COMM_WORLD, status, ierr)
      CALL MPI_Bcast(Send_Config, dimensoes+2, MPI_DOUBLE_PRECISION, 0, &
                     MPI_COMM_WORLD, ierr)
      DO i=1, dimensoes
        Best_Config%solucao(i) = Send_Config(i)
      ENDDO
      Best_Config%custo = Send_Config(dimensoes+1)
      tot_nafo = INT(Send_Config(dimensoes+2))
      nafo = 0
    ENDIF

  END SUBROUTINE

!Fim do programa
END PROGRAM

