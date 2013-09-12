!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Instituto Nacional de Pesquisas Espaciais - INPE
! Programa de Doutorado em Computação Aplicada - CAP
!
! Rede neural Perceptron de Multiplas Camadas (MLP) auto-configurada
! utilizando o Multiple Particle Collision Algorithm (MPCA)
! 
! Desenvolvedores:  Juliana A Anochi (CAP/INPE)
!					Sabrina Sambatti (CAP/INPE)
!
! Email: [juliana.anochi@lac.inpe.br, sabrinabms@gmail.com]
!
! Consulte o LEIA-ME para mais informações sobre a RNA&MPCA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


MODULE ModuloRNA

	CONTAINS

	!Funcao para treinar a rede neural com os parametros calculados pela heuristica
	DOUBLE PRECISION FUNCTION  Rede_Neural_BP(parametros)
	
		IMPLICIT NONE
		DOUBLE PRECISION :: parametros(16,1)

		CHARACTER*24, caminho
		CHARACTER*24, arq_ent, arq_sai, arq_result, arq_ent_valid, arq_sai_valid
		CHARACTER*24, arq_wh1, arq_wh2, arq_wh3, arq_ws
		CHARACTER*24, arq_bh1, arq_bh2, arq_bh3, arq_bs
		CHARACTER*52, path_ent, path_sai
		CHARACTER*52, path_wh1, path_wh2, path_wh3, path_ws
		CHARACTER*52, path_bh1, path_bh2, path_bh3, path_bs
		CHARACTER*52, path_sai_valid, path_ent_valid, path_result
		INTEGER :: hide_camada

		hide_camada = parametros(1,1)

		!Definindo os nomes dos arquivos de entrada, saida e resultados
		arq_ent ='x.txt'
		arq_sai ='yd.txt'
		arq_sai_valid = 'yd_valid.txt'
		arq_ent_valid = 'x_valid.txt'
		arq_wh1 = 'wh1.txt'
		arq_wh2 = 'wh2.txt'
		arq_wh3 = 'wh3.txt'
		arq_ws  = 'ws.txt'
		arq_bh1 = 'bh1.txt'
		arq_bh2 = 'bh2.txt'
		arq_bh3 = 'bh3.txt'
		arq_bs  = 'bs.txt'
		arq_result ='compara.txt'

		!Definindo o caminho (path)
		caminho = './dados/'
		path_ent = TRIM(caminho)//TRIM(arq_ent)
		path_sai = TRIM(caminho)//TRIM(arq_sai)
		path_sai_valid = TRIM(caminho)//TRIM(arq_sai_valid)
		path_ent_valid = TRIM(caminho)//TRIM(arq_ent_valid)
		path_wh1 = TRIM(caminho)//TRIM(arq_wh1)
		path_wh2 = TRIM(caminho)//TRIM(arq_wh2)
		path_wh3 = TRIM(caminho)//TRIM(arq_wh3)
		path_ws = TRIM(caminho)//TRIM(arq_ws)
		path_bh1 = TRIM(caminho)//TRIM(arq_bh1)
		path_bh2 = TRIM(caminho)//TRIM(arq_bh2)
		path_bh3 = TRIM(caminho)//TRIM(arq_bh3)
		path_bs = TRIM(caminho)//TRIM(arq_bs)
		path_result =TRIM(caminho)//TRIM(arq_result)

		IF (hide_camada .EQ. 1) THEN
			Rede_Neural_BP = Rede_Neural_C1(parametros, path_ent, path_sai, path_ent_valid, path_sai_valid,&
			path_wh1, path_bh1, path_ws, path_bs, path_result)
		ELSE IF (hide_camada .EQ. 2) THEN
			Rede_Neural_BP = Rede_Neural_C2(parametros, path_ent, path_sai, path_ent_valid, path_sai_valid,&
			path_wh1, path_wh2, path_ws, path_bh1, path_bh2, path_bs, path_result)
		ELSE IF (hide_camada .EQ. 3) THEN
			Rede_Neural_BP = Rede_Neural_C3(parametros, path_ent, path_sai, path_ent_valid, path_sai_valid,&
			path_wh1, path_wh2, path_wh3, path_ws, path_bh1, path_bh2, path_bh3, path_bs, path_result)
		ENDIF

	END FUNCTION Rede_Neural_BP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	DOUBLE PRECISION FUNCTION Rede_Neural_C1(parametros, path_ent, path_sai, path_ent_valid, path_sai_valid,&
		path_wh1, path_bh1, path_ws, path_bs, path_result)

		IMPLICIT NONE
	
		!Parametros recebidos pela funcao (passados pelo algoritmo de otimizacao)
		double precision, dimension(:,:) :: parametros
		character*52, path_ent,path_sai
		character*52, path_wh1, path_bh1, path_ws, path_bs
		character*52, path_sai_valid, path_ent_valid, path_result

		!Variaveis usadas para o calculo da funcao objetivo
		double precision :: penaltyObj, alphaObj, betaObj
	
		double precision, allocatable,dimension(:,:) :: x
		double precision, allocatable,dimension(:,:) :: yd
		double precision, allocatable,dimension(:,:) :: x_valid
		double precision, allocatable,dimension(:,:) :: yd_valid
		double precision, allocatable,dimension(:,:) :: vs, ys, ys_melhor
		double precision, allocatable,dimension(:,:) :: ws, ws_velho, ws_menor
		double precision, allocatable,dimension(:,:) :: bs, bs_velho, bs_menor
		double precision, allocatable,dimension(:,:) :: eqm_menor, eqm_valid_menor
		double precision, allocatable,dimension(:,:) :: erro
		double precision, allocatable,dimension(:,:) :: erro_pad
		double precision, allocatable,dimension(:,:) :: eqm
		double precision, allocatable,dimension(:,:) :: eqm_valid
		double precision, allocatable,dimension(:,:) :: grad_saida, grad_hide1
		double precision, allocatable,dimension(:,:) :: deltaw_saida
		double precision, allocatable,dimension(:,:) :: deltabs, deltabh1
		double precision, allocatable,dimension(:,:) :: vh1, vh2
		double precision, allocatable,dimension(:,:) :: wh1, wh1_velho, wh1_menor
		double precision, allocatable,dimension(:,:) :: bh1, bh1_velho, bh1_menor
		double precision, allocatable,dimension(:,:) :: yh1
		double precision, allocatable,dimension(:,:) :: deltaw_h1

		double precision, parameter :: a = 1

		double precision	:: entrada(16,1)
		double precision	:: dv,soma, erro_menor,alpha
		double precision   	:: eta		! taxa de aprendizado
		double precision   	:: ed		! erro desejado
		double precision   	:: aleat    ! gerar pesos e bias aleatoriamente

		integer		:: hide_neuron(3,1) ! numero de neuronios na camada escondida
		integer		:: cont, l, unit
		integer		:: max_it	     	! numero maximo de iteracoes
		integer		:: vetor_entrada	! dimensao do vetor de entrada
		integer		:: num_pad		    ! numero de padroes a serem treinados
		integer		:: num_pad_valid
		integer		:: vetor_saida	 	! dimensao do vetor de saida
		integer		:: i, j, k, hn		! variavel para controle do programa
		integer		:: f_ativa		    ! controle de qual funcao ativacao
		integer		:: f_deriva		    ! controle de qual funcao derivada
		integer		:: hide_camada		! numero de camadas escondidas
		integer		:: mpca			    ! variavel para ativar o mpca (mpca = 0 --> sequencial, mpca = 1 --> paralelo)
		integer		:: peso		      	! variavel para controlar como os pesos serao inicializados (0=fixo, 1=aleat)
		integer		:: validacao		! variavel para controlar se terah ou nao validacao (validacao = 1 --> com validacao cruzada)

		!Atribuicao empirica dos pesos para a funcao objetivo
		penaltyObj = 1.0D0
		alphaObj = 1.0D0
		betaObj = 1.0D0
	
		DO I = 1, 16
			entrada(I,1) = parametros(I,1)
		END DO

		hide_camada = entrada(1,1)
		hide_neuron(1,1) = entrada(2,1)	!numero de neuronios na camada escondida
		hide_neuron(2,1) = entrada(3,1)	!numero de neuronios na camada escondida
		hide_neuron(3,1) = entrada(4,1)	!numero de neuronios na camada escondida
		f_ativa = entrada(5,1)		    !parametro que corresponde ao tipo de funcao que sera utilizada
		alpha = entrada(6,1)  		    !termo momentum
		eta = entrada(7,1)
		num_pad = entrada(8,1)
		num_pad_valid = entrada(9,1)
		vetor_entrada = entrada(10,1)	!numero de entradas
		vetor_saida = entrada(11,1)	    !numero de saidas		
		ed = entrada(12,1)	        	!erro desejado
		max_it = entrada(13,1)	    	!numero maximo de iteracao
		mpca = entrada(14,1)		    !mpca = 0 --> serial, mpca = 1 --> paralelo
		peso = entrada(15,1)	        !peso = 0 --> valor fixo, peso = 1 --> valor aleatorio
		validacao = entrada(16,1)	    ! = 0 --> sem validacao, = 1 --> com validacao
		f_deriva = f_ativa

		!inicializando o menor erro
		erro_menor = 10**9
		
		!Bias 1. camada oculta
		allocate(bh1(hide_neuron(1,1),1))
		allocate(bh1_menor(hide_neuron(1,1),1))
		allocate(bh1_velho(hide_neuron(1,1),1))
		!Bias camada de saida
		allocate(bs(vetor_saida,1))
		allocate(bs_menor(vetor_saida,1))
		allocate(bs_velho(vetor_saida,1))
		!DeltaBias
		allocate(deltabh1(hide_neuron(1,1),1))
		allocate(deltabs(vetor_saida,1))
		!DeltaW's
		allocate(deltaw_h1(vetor_entrada,hide_neuron(1,1)))
		allocate(deltaw_saida(hide_neuron(hide_camada,1),vetor_saida))
		!Erros
		allocate(eqm(1,max_it))
		allocate(eqm_valid(1,max_it))
		allocate(eqm_valid_menor(1,max_it))
		allocate(eqm_menor(1,max_it))
		allocate(erro(vetor_saida,num_pad))
		allocate(erro_pad(num_pad,1))
		!Gradientes
		allocate(grad_hide1(hide_neuron(1,1),1))
		allocate(grad_saida(vetor_saida,1))
		!Campo local induzido
		allocate(vh1(hide_neuron(1,1),1))
		allocate(vs(vetor_saida,1))
		!Pesos 1. camada oculta
		allocate(wh1(vetor_entrada,hide_neuron(1,1)))
		allocate(wh1_menor(vetor_entrada,hide_neuron(1,1)))
		allocate(wh1_velho(vetor_entrada,hide_neuron(1,1)))
		!Pesos camada de saida
		allocate(ws(hide_neuron(hide_camada,1),vetor_saida))
		allocate(ws_menor(hide_neuron(hide_camada,1),vetor_saida))
		allocate(ws_velho(hide_neuron(hide_camada,1),vetor_saida))
		!Dados de entrada e saida
		allocate(x(vetor_entrada,num_pad))
		allocate(x_valid(vetor_entrada,num_pad_valid))
		allocate(yd(vetor_saida,num_pad))
		allocate(yd_valid(vetor_saida,num_pad_valid))
		!Saidas obtidas
		allocate(yh1(hide_neuron(1,1),1))
		allocate(ys(vetor_saida,num_pad))
		allocate(ys_melhor(vetor_saida,num_pad)) 
		
		!------------------------------------------------------------!
		!LENDO OS PARAMETROS DO ARQUIVO DE ENTRADA
		!------------------------------------------------------------!
		OPEN (1, file = path_sai)
		DO I = 1, vetor_saida
			READ(1,*) (yd(I,J), J = 1, num_pad)
		END DO
		CLOSE (1)

		OPEN (2, file = path_ent)
		DO I = 1, vetor_entrada
			READ(2,*) (x(I,J), J = 1, num_pad)
		END DO
		CLOSE (2)

		IF (validacao .EQ. 1) THEN
			OPEN (1, file = path_sai_valid)
			DO I = 1, vetor_saida
				READ(1,*) (yd_valid(I,J), J = 1, num_pad_valid)
			END DO
			CLOSE (1)

			OPEN (2, file = path_ent_valid)
			DO I = 1, vetor_entrada
				READ(2,*) (x_valid(I,J), J=1, num_pad_valid)
			END DO
			CLOSE (2)
		ENDIF
		!--------------------------------------------------------------------!
		!INICIALIZANDO OS PARAMETROS: wh1, bh1, ws e bs
		!--------------------------------------------------------------------!

		!PRIMEIRA CAMADA OCULTA
        	DO l = 1,vetor_entrada
			DO k = 1, hide_neuron(1,1)
				if (peso .EQ. 0) then
					wh1(l,k) = 0.5
				else 
					call random_number(aleat)
					wh1(l,k) = aleat
				endif
            		ENDDO
        	ENDDO
		DO k = 1, hide_neuron(1,1)
			if (peso .EQ. 0) then
				bh1(k,1) = 0.5
			else
				call random_number(aleat)
				bh1(k,1) = aleat
			endif
		ENDDO

		!CAMADA DE SAIDA
		DO l = 1,hide_neuron(hide_camada,1)
			DO k = 1, vetor_saida
				if (peso .EQ. 0) then
					ws(l,k) = 0.5
				else
					call random_number(aleat)
					ws(l,k) = aleat
				endif
			ENDDO
		ENDDO
		DO k = 1,vetor_saida
			if (peso .EQ. 0) then
				bs(k,1) = 0.5
			else
				call random_number(aleat)
				bs(k,1) = aleat
			endif
		ENDDO

		!----------------------------------------------------------------------!
		! INICIO DA REDE: FEEDFORWARD
		!----------------------------------------------------------------------!
		cont = 1	!Contador de epocas, usado para alterar o valor do ETA
		l = 0
	
		!DO WHILE ((erro_menor .GE. ed) .AND. (l .LT. max_it))
		DO WHILE (l .LT. max_it)

			l = l + 1

			DO i = 1, num_pad	!DO NUMERO TOTAL DE PADROES

				! ATIVACAO: 1. CAMADA OCULTA
				vh1 = 0.d0
 				vh1(:,1) = matmul(x(:,i), wh1(:,:))
				vh1(:,1) = vh1(:,1) - bh1(:,1);
				select case(f_ativa)
					case (1) !LOGISTICA
						yh1(:,1) = 1.d0/(1.d0+DEXP(-a * vh1(:,1)))
					case (2) !TANGENTE
						yh1(:,1) = (1.d0-DEXP(-vh1(:,1)))/(1.d0+DEXP(-vh1(:,1)))
					case (3) !GAUSS
						yh1(:,1) = DEXP(-vh1(:,1)**2)
				end select

				! ATIVACAO: CAMADA DE SAIDA
				vs = 0.d0
				vs(:,1) = matmul(yh1(:,1),ws(:,:))
				vs(:,1) = vs(:,1) - bs(:,1)

				select case(f_ativa)
					case (1) !LOGISTICA
						ys(:,i) = 1.d0/(1.d0+DEXP(-a*vs(:,1)))
					case (2) !TANGENTE
						ys(:,i) = (1.d0-DEXP(-vs(:,1)))/(1.d0+DEXP(-vs(:,1)))
					case (3) !GAUSS
						ys(:,i) = DEXP(-vs(:,1)**2)
				end select

				! PARA O CALCULO DO NOVO PESO
				wh1_velho = wh1
				bh1_velho = bh1
				ws_velho = ws
				bs_velho = bs
				
				! CALCULO ERRO TREINAMENTO
				erro(:,i) = yd(:,i) - ys(:,i)

		!-------------------------------------------------------------------------!
		!                        BACKPROPAGATION
		!-------------------------------------------------------------------------!
				! TREINAMENTO: CAMADA DE SAIDA
				DO j = 1,vetor_saida

					select case(f_deriva)
						case (1)!LOGISTICA: e^-x / (1+e^-x)^2
							dv = ((a*DEXP(-a*vs(j,1)))/((1.d0+DEXP(-a*vs(j,1)))**2.d0))
						case (2)!TANGENTE: 2e^-x / (1+e^-x)^2
							dv = (2*DEXP(-vs(j,1)))/ ((1+DEXP(-vs(j,1)))**2)
						case (3)!GAUSS
							dv = -ys(j,i)/a	
					end select

					grad_saida(j,1) = erro(j,i) * dv
					deltaw_saida(:,j) = eta * grad_saida(j,1) * yh1(:,1)
					ws(:,j) = ws(:,j) + alpha * (ws(:,j) - ws_velho(:,j)) + deltaw_saida(:,j)
		   			deltabs(j,1) = eta * grad_saida(j,1) * (-1.d0)
		   			bs(j,1) = bs(j,1) + deltabs(j,1)

				enddo !CAMADA SAIDA

				! TREINAMENTO: 1. CAMADA OCULTA
				DO j = 1,hide_neuron(1,1)

					soma = 0.d0
					DO k = 1,vetor_saida
						soma = soma + ( grad_saida(k,1) * ws_velho(j,k) )
					enddo

					select case(f_deriva)
						case (1)!LOGISTICA: e^-x / (1+e^-x)^2
							dv = ((a*DEXP(-a*vh1(j,1)))/((1.d0+DEXP(-a*vh1(j,1)))**2.d0))
						case (2)!TANGENTE
							dv = (2*DEXP(-vh1(j,1)))/ ((1+DEXP(-vh1(j,1)))**2)
						case (3)!GAUSS
							dv = -yh1(j,1)/a	
					end select

					grad_hide1(j,1) = dv * soma
   					deltaw_h1(:,j) = eta * grad_hide1(j,1) * x(:,i)
					wh1(:,j) = wh1(:,j) + alpha * (wh1(:,j) - wh1_velho(:,j)) + deltaw_h1(:,j)      
					deltabh1(j,1) = eta * grad_hide1(j,1) * (-1.d0)
  					bh1(j,1) = bh1(j,1) + deltabh1(j,1)

				ENDDO ! 1. CAMADA OCULTA

	
				!CALCULO PADRAO DO ERRO!!!!
					erro_pad(i,1) = sum(erro(:,i),dim=1)
					erro_pad(i,1) = 0.5d0*(erro_pad(i,1)**2.d0)
				
			ENDDO !NUMERO PADROES
			eqm(1,l) = sum(erro_pad(:,1))
			eqm(1,l) = (1.d0/(num_pad))*eqm(1,l)
			
			
			!***********************************************************************
			! VALIDACAO CRUZADA
			!***********************************************************************

			DO i = 1,num_pad_valid

				! ATIVACAO: 1. CAMADA OCULTA
				vh1 = 0.d0
 				vh1(:,1) = matmul(x_valid(:,i),wh1(:,:))
				vh1(:,1) = vh1(:,1) - bh1(:,1);
			
				select case(f_ativa)
					case (1) !LOGISTICA
						yh1(:,1) = 1.d0/(1.d0+DEXP(-a * vh1(:,1)))
					case (2) !TANGENTE
						yh1(:,1) = (1.d0-DEXP(-vh1(:,1)))/(1.d0+DEXP(-vh1(:,1)))
					case (3) !GAUSS
						yh1(:,1) = DEXP(-vh1(:,1)**2)
				end select

				! ATIVACAO: CAMADA DE SAIDA
				vs = 0.d0
				vs(:,1) = matmul(yh1(:,1),ws(:,:))
				vs(:,1) = vs(:,1) - bs(:,1)

				select case(f_ativa)
					case (1) !LOGISTICA
						ys(:,i) = 1.d0/(1.d0+DEXP(-a*vs(:,1)))
					case (2) !TANGENTE
						ys(:,i) = (1.d0-DEXP(-vs(:,1)))/(1.d0+DEXP(-vs(:,1)))
					case (3) !GAUSS
						ys(:,i) = DEXP(-vs(:,1)**2)
				end select

				!CALCULO DO ERRO RNA
				erro(:,i) = yd_valid(:,i) - ys(:,i)
				
				!CALCULO PADRAO DO ERRO!!!!
					erro_pad(i,1) = sum(erro(:,i),dim=1)
					erro_pad(i,1) = 0.5d0*(erro_pad(i,1)**2.d0)

			ENDDO !DO VALIDACAO

			!CALCULO DOS ERROS
			eqm_valid(1,l) = sum(erro_pad(:,1))
			eqm_valid(1,l) = (1.d0/(num_pad_valid))*eqm_valid(1,l)

			if (eqm_valid(1,l) < erro_menor) then
				erro_menor = eqm_valid(1,l)
				wh1_menor = wh1
				bh1_menor = bh1
				ys_melhor = ys
				ws_menor = ws
				bs_menor = bs
				eqm_menor = eqm
				eqm_valid_menor = eqm_valid
			endif

			if (cont >= 100) then
				if (mpca .EQ. 0) then
					WRITE(*,*)l,eqm(1,l), eqm_valid(1,l)
				endif
				cont = 1
				eta = eta * 0.99
			else
				cont = cont+1
			endif
				
		ENDDO !DO MAXIMO ITERAÇÃO

		IF ( mpca .EQ. 0 ) THEN

			Open(12, FILE = path_result)!,status="REPLACE", ACCESS="SEQUENTIAL",ACTION="write")
				do i=1,vetor_saida
					write(12,110)(ys(i,j),j=1,num_pad_valid)
				enddo
			close(12)

			Open(12, FILE = path_wh1)!,status="REPLACE", ACCESS="SEQUENTIAL",ACTION="write")
				do i = 1, vetor_entrada
					write(12,110)(wh1_menor(i,j), j = 1, hide_neuron(1,1))
				enddo
			close(12)

			Open(12, FILE = path_bh1)!,status="REPLACE", ACCESS="SEQUENTIAL",ACTION="write")
				do i = 1, 1
					write(12,110)(bh1_menor(j,i), j = 1, hide_neuron(1,1))
				enddo
			close(12)

			Open(12, FILE = path_ws)!,status="REPLACE", ACCESS="SEQUENTIAL",ACTION="write")
				do i = 1, hide_neuron(hide_camada,1)
					write(12,110)(ws_menor(i,j), j = 1, vetor_saida)
				enddo
			close(12)

			Open(12, FILE = path_bs)!,status="REPLACE", ACCESS="SEQUENTIAL",ACTION="write")
				do i = 1, 1
					write(12,110)(bs_menor(j,i), j = 1, vetor_saida)
				enddo
			close(12)
110			format (424F32.16)
		ENDIF

		!CALCULO DA FUNCAO OBJETIVO DEFINIDA POR: Adenilson
		IF (validacao .EQ. 1) THEN
			Rede_Neural_C1 = penaltyObj * ( (alphaObj * eqm(1,l) + betaObj * &
                	eqm_valid(1,l)) / (alphaObj + betaObj) )
		ELSE
			 Rede_Neural_C1 = penaltyObj * eqm(1,l)
		ENDIF

999		deallocate(x,yd,x_valid,yd_valid,vh1,wh1,yh1)
		deallocate(wh1_velho, wh1_menor,bh1,bh1_velho,bh1_menor)
		deallocate(vs,ys,ys_melhor,ws,ws_velho,bs)
		deallocate(bs_velho,bs_menor,erro,erro_pad,grad_saida,grad_hide1)
		deallocate(deltabs,deltabh1,deltaw_saida,deltaw_h1,eqm)
		deallocate(eqm_valid,eqm_menor,eqm_valid_menor)
		

	!FIM DA FUNCAO
	END FUNCTION Rede_Neural_C1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	DOUBLE PRECISION FUNCTION Rede_Neural_C2(parametros, path_ent, path_sai, path_ent_valid, path_sai_valid,&
		path_wh1, path_wh2, path_ws, path_bh1, path_bh2, path_bs, path_result)

		IMPLICIT NONE
	
		!Parametros recebidos pela funcao (passados pelo algoritmo de otimizacao)
		double precision, dimension(:,:) :: parametros
		character*52, path_ent,path_sai
		character*52, path_wh1, path_wh2, path_bh1, path_bh2, path_ws, path_bs
		character*52, path_sai_valid, path_ent_valid, path_result	
	
		!Variaveis usadas para o calculo da funcao objetivo
		double precision :: penaltyObj, alphaObj, betaObj

		double precision, allocatable,dimension(:,:) :: x
		double precision, allocatable,dimension(:,:) :: yd
		double precision, allocatable,dimension(:,:) :: x_valid
		double precision, allocatable,dimension(:,:) :: yd_valid
		double precision, allocatable,dimension(:,:) :: vs, ys, ys_melhor
		double precision, allocatable,dimension(:,:) :: ws, ws_velho, ws_menor
		double precision, allocatable,dimension(:,:) :: bs, bs_velho, bs_menor
		double precision, allocatable,dimension(:,:) :: eqm_menor, eqm_valid_menor
		double precision, allocatable,dimension(:,:) :: erro
		double precision, allocatable,dimension(:,:) :: erro_pad
		double precision, allocatable,dimension(:,:) :: eqm
		double precision, allocatable,dimension(:,:) :: eqm_valid
		double precision, allocatable,dimension(:,:) :: grad_saida, grad_hide1, grad_hide2
		double precision, allocatable,dimension(:,:) :: deltaw_saida
		double precision, allocatable,dimension(:,:) :: deltabs, deltabh1, deltabh2
		double precision, allocatable,dimension(:,:) :: vh1, vh2
		double precision, allocatable,dimension(:,:) :: wh1, wh1_velho, wh1_menor
		double precision, allocatable,dimension(:,:) :: wh2, wh2_velho, wh2_menor
		double precision, allocatable,dimension(:,:) :: bh1, bh1_velho, bh1_menor
		double precision, allocatable,dimension(:,:) :: bh2, bh2_velho, bh2_menor
		double precision, allocatable,dimension(:,:) :: yh1, yh2
		double precision, allocatable,dimension(:,:) :: deltaw_h1, deltaw_h2

		double precision, parameter :: a = 1

		!Parametros recebidos pela funcao (passados pelo algoritmo de otimizacao)
		double precision	:: entrada(16,1)
		double precision	:: dv,soma, erro_menor,alpha
		double precision    :: eta			! taxa de aprendizado
		double precision    :: ed			! erro desejado
		double precision    :: aleat        ! gerar pesos e bias aleatoriamente

		integer		:: hide_neuron(3,1) ! número de neurônios na camada escondida
		integer		:: cont, l, unit
		integer		:: max_it	     	! número máximo de iterações
		integer		:: vetor_entrada	! dimensão do vetor de entrada
		integer		:: num_pad			! número de padrões a serem treinados
		integer		:: num_pad_valid
		integer		:: vetor_saida	 	! dimensão do vetor de saída
		integer		:: i, j, k, hn		! variável para controle do programa
		integer		:: f_ativa			! controle de qual funcao ativacao
		integer		:: f_deriva			! controle de qual funcao derivada
		integer		:: hide_camada		! numero de camadas escondidas
		integer		:: mpca				! variavel para ativar o mpca (mpca = 0 --> sequencial, mpca = 1 --> paralelo)
		integer		:: peso				! variavel para controlar como os pesos serao inicializados (0=fixo, 1=aleat)
		integer		:: validacao		! variavel para controlar se terah ou nao validacao cruzada (validacao = 1 --> com validacao cruzada)


		!Atribuicao empirica dos pesos para a funcao objetivo
		penaltyObj = 1.0D0
		alphaObj = 1.0D0
		betaObj = 1.0D0

		DO I = 1, 16
			entrada(I,1) = parametros(I,1)
		END DO

		hide_camada = entrada(1,1)
		hide_neuron(1,1) = entrada(2,1)	!numero de neuronios na camada escondida
		hide_neuron(2,1) = entrada(3,1)	!numero de neuronios na camada escondida
		hide_neuron(3,1) = entrada(4,1)	!numero de neuronios na camada escondida
		f_ativa = entrada(5,1)			!parametro que corresponde ao tipo de funcao que sera utilizada
		alpha = entrada(6,1)			!termo momentum
		eta = entrada(7,1)		
		num_pad = entrada(8,1)
		num_pad_valid = entrada(9,1)
		vetor_entrada = entrada(10,1)	!numero de entradas
		vetor_saida = entrada(11,1)		!numero de saidas		
		ed = entrada(12,1)				!erro desejado
		max_it = entrada(13,1)			!numero maximo de iteração
		mpca = entrada(14,1)			!mpca = 0 --> serial, mpca = 1 --> paralelo
		peso = entrada(15,1)			!peso = 0 --> valor fixo, peso = 1 --> valor aleatorio
		validacao = entrada(16,1)		!validacao = 0 --> sem validacao, validacao = 1 --> com validacao

		f_deriva = f_ativa

		erro_menor = 10**9		!inicializando o menor erro

		!Bias 1. camada escondida
		allocate(bh1(hide_neuron(1,1),1))
		allocate(bh1_menor(hide_neuron(1,1),1))
		allocate(bh1_velho(hide_neuron(1,1),1))
		!Bias 2. camada escondida
		allocate(bh2(hide_neuron(2,1),1))
		allocate(bh2_menor(hide_neuron(2,1),1))
		allocate(bh2_velho(hide_neuron(2,1),1))
		!Bias camada de saida
		allocate(bs(vetor_saida,1))
		allocate(bs_menor(vetor_saida,1))
		allocate(bs_velho(vetor_saida,1))
		!DeltaBias
		allocate(deltabh1(hide_neuron(1,1),1))
		allocate(deltabh2(hide_neuron(2,1),1))
		allocate(deltabs(vetor_saida,1))
		!DeltaW's
		allocate(deltaw_h1(vetor_entrada,hide_neuron(1,1)))
		allocate(deltaw_h2(hide_neuron(1,1),hide_neuron(2,1)))
		allocate(deltaw_saida(hide_neuron(hide_camada,1),vetor_saida))
		!Erros
		allocate(eqm(1,max_it))
		allocate(eqm_valid(1,max_it))
		allocate(eqm_valid_menor(1,max_it))
		allocate(eqm_menor(1,max_it))
		allocate(erro(vetor_saida,num_pad))
		allocate(erro_pad(num_pad,1))
		!Gradientes
		allocate(grad_hide1(hide_neuron(1,1),1))
		allocate(grad_hide2(hide_neuron(2,1),1))
		allocate(grad_saida(vetor_saida,1))
		!Campo local induzido
		allocate(vh1(hide_neuron(1,1),1))
		allocate(vh2(hide_neuron(2,1),1))
		allocate(vs(vetor_saida,1))
		!Pesos 1. camada escondida
		allocate(wh1(vetor_entrada,hide_neuron(1,1)))
		allocate(wh1_menor(vetor_entrada,hide_neuron(1,1)))
		allocate(wh1_velho(vetor_entrada,hide_neuron(1,1)))
		!Pesos 2. camada escondida
		allocate(wh2(hide_neuron(1,1),hide_neuron(2,1)))
		allocate(wh2_menor(hide_neuron(1,1),hide_neuron(2,1)))
		allocate(wh2_velho(hide_neuron(1,1),hide_neuron(2,1)))
		!Pesos camada de saida
		allocate(ws(hide_neuron(hide_camada,1),vetor_saida))
		allocate(ws_menor(hide_neuron(hide_camada,1),vetor_saida))
		allocate(ws_velho(hide_neuron(hide_camada,1),vetor_saida))
		!Dados de entrada e saida
		allocate(x(vetor_entrada,num_pad))
		allocate(x_valid(vetor_entrada,num_pad_valid))
		allocate(yd(vetor_saida,num_pad))
		allocate(yd_valid(vetor_saida,num_pad_valid))
		!Saidas obtidas
		allocate(yh1(hide_neuron(1,1),1))
		allocate(yh2(hide_neuron(2,1),1))
		allocate(ys(vetor_saida,num_pad))
		allocate(ys_melhor(vetor_saida,num_pad)) 

		!------------------------------------------------------------!
		!LENDO OS PARAMETROS DO ARQUIVO DE ENTRADA
		!------------------------------------------------------------!
		OPEN (1, file = path_sai)
		DO I = 1, vetor_saida
			READ(1,*) (yd(I,J), J = 1, num_pad)
		END DO
		CLOSE (1)

		OPEN (2, file = path_ent)
		DO I = 1, vetor_entrada
			READ(2,*) (x(I,J), J = 1, num_pad)
		END DO
		CLOSE (2)

		IF (validacao .EQ. 1) THEN
			OPEN (1, file = path_sai_valid)
			DO I = 1, vetor_saida
				READ(1,*) (yd_valid(I,J), J = 1, num_pad_valid)
			END DO
			CLOSE (1)

			OPEN (2, file = path_ent_valid)
			DO I = 1, vetor_entrada
				READ(2,*) (x_valid(I,J), J=1, num_pad_valid)
			END DO
			CLOSE (2)
		ENDIF

		!--------------------------------------------------------------------!
		!INICIALIZANDO OS PARAMETROS: wh1, bh1, ws e bs
		!--------------------------------------------------------------------!

		!PRIMEIRA CAMADA OCULTA
	        DO l = 1,vetor_entrada
			DO k = 1, hide_neuron(1,1)
				if (peso .EQ. 0) then
					wh1(l,k) = 0.5
				else 
					call random_number(aleat)
					wh1(l,k) = aleat
				endif
        	    ENDDO
        	ENDDO
		DO k = 1, hide_neuron(1,1)
			if (peso .EQ. 0) then
				bh1(k,1) = 0.5
			else
				call random_number(aleat)
				bh1(k,1) = aleat
			endif
		ENDDO

		!SEGUNDA CAMADA OCULTA
		DO k = 1, hide_neuron(1,1)
			DO i = 1,hide_neuron(2,1)
				if (peso .EQ. 0) then
					wh2(k,i) = 0.5
				else
					call random_number(aleat)
					wh2(k,i) = aleat
				endif
			ENDDO
		ENDDO
		DO i = 1,hide_neuron(2,1)
			if (peso .EQ. 0) then
				bh2(i,1) = 0.5
			else
                call random_number(aleat)
                bh2(i,1) = aleat
			endif
		ENDDO

		!CAMADA DE SAIDA
		DO l = 1,hide_neuron(hide_camada,1)
			DO k = 1, vetor_saida
				if (peso .EQ. 0) then
					ws(l,k) = 0.5
				else
	                call random_number(aleat)
		            ws(l,k) = aleat
				endif
			ENDDO
		ENDDO
		DO k = 1,vetor_saida
			if (peso .EQ. 0) then
				bs(k,1) = 0.5
			else
				call random_number(aleat)
				bs(k,1) = aleat
			endif
		ENDDO

		!----------------------------------------------------------------------!
		! INICIO DA REDE: FEEDFORWARD
		!----------------------------------------------------------------------!
		cont = 1	!Contador de epocas
		l = 0 

		!DO WHILE ((erro_menor .GE. ed) .AND. (l .LT. max_it))
		DO WHILE (l .LT. max_it)

			l = l + 1

			DO i = 1, num_pad	!DO NUMERO TOTAL DE PADROES
				
				! ATIVACAO: 1. CAMADA OCULTA
				vh1 = 0.d0
 				vh1(:,1) = matmul(x(:,i),wh1(:,:))
				vh1(:,1) = vh1(:,1) - bh1(:,1);
			
				SELECT CASE(f_ativa)
					CASE (1) !LOGISTICA
						yh1(:,1) = 1.d0/(1.d0+DEXP(-a * vh1(:,1)))
					CASE (2) !TANGENTE
						yh1(:,1) = (1.d0-DEXP(-vh1(:,1)))/(1.d0+DEXP(-vh1(:,1)))
					CASE (3) !GAUSS
						yh1(:,1) = DEXP(-vh1(:,1)**2)
				END SELECT

				! ATIVACAO: 2. CAMADA OCULTA
				vh2(:,1) = 0.0
				vh2(:,1) = MATMUL(yh1(:,1),wh2(:,:))
				vh2(:,1) = vh2(:,1) - bh2(:,1)
				SELECT CASE(f_ativa)
					CASE (1) !LOGISTICA
						yh2(:,1) = 1.d0/(1.d0+DEXP(-a*vh2(:,1)))
					CASE (2) !TANGENTE
						yh2(:,1) = (1.d0-DEXP(-vh2(:,1)))/(1.d0+DEXP(-vh2(:,1)))
					CASE (3) !GAUSS
						yh2(:,1) = DEXP(-vh2(:,1)**2)
				END SELECT
					
				! ATIVACAO: CAMADA DE SAIDA
				vs = 0.d0
				vs(:,1) = MATMUL(yh2(:,1), ws(:,:))
				vs(:,1) = vs(:,1) - bs(:,1)
				SELECT CASE(f_ativa)
					CASE (1) !LOGISTICA
						ys(:,i) = 1.d0 / (1.d0 + DEXP(-a * vs(:,1)))
					CASE (2) !TANGENTE
						ys(:,i) = (1.d0 - DEXP(-vs(:,1))) / (1.d0 + DEXP(-vs(:,1)))
					CASE (3) !GAUSS
						ys(:,i) = DEXP(-vs(:,1)**2)
				END SELECT

				! PARA CALCULO DO NOVO PESO
				wh1_velho = wh1
				bh1_velho = bh1
				wh2_velho = wh2
				bh2_velho = bh2
				ws_velho = ws
				bs_velho = bs
					
				! CALCULO ERRO TREINAMENTO
				erro(:,i) = yd(:,i) - ys(:,i)

		!-------------------------------------------------------------------------!
		!                        BACKPROPAGATION
		!-------------------------------------------------------------------------!
				! TREINAMENTO: CAMADA DE SAIDA
				DO j = 1,vetor_saida

					SELECT CASE(f_deriva)
						CASE (1)!LOGISTICA: e^-x / (1+e^-x)^2
							dv = ((a * DEXP(-a * vs(j,1))) / ((1.d0 + DEXP(-a * vs(j,1))) ** 2.d0))
						CASE (2)!TANGENTE: 2e^-x / (1+e^-x)^2
							dv = (2 * DEXP(-vs(j,1))) / ((1 + DEXP(-vs(j,1))) ** 2)
						CASE (3)!GAUSS
							dv = -ys(j,i) / a	
					END SELECT

					grad_saida(j,1) = erro(j,i) * dv
					deltaw_saida(:,j) = eta * grad_saida(j,1) * yh2(:,1)
					ws(:,j) = ws(:,j) + alpha * (ws(:,j) - ws_velho(:,j)) + deltaw_saida(:,j)
		   			deltabs(j,1) = eta * grad_saida(j,1) * (-1.d0)
		   			bs(j,1) = bs(j,1) + deltabs(j,1)

				enddo !CAMADA SAIDA

				!TREINAMENTO: 2. CAMADA OCULTA
				DO j = 1,hide_neuron(2,1)

					soma = 0.d0
					do k = 1, vetor_saida
						soma = soma + ( grad_saida(k,1) * ws_velho(j,k) )
					enddo

					select case(f_deriva)
						case (1)!LOGISTICA: e^-x / (1+e^-x)^2
							dv = ((a*DEXP(-a*vh2(j,1)))/((1.d0+DEXP(-a*vh2(j,1)))**2.d0))
						case (2)!TANGENTE
							dv = (2*DEXP(-vh2(j,1)))/ ((1+DEXP(-vh2(j,1)))**2)
						case (3)!GAUSS
							dv = -yh2(j,1)/a
					end select
			
					grad_hide2(j,1) = dv * soma
					deltaw_h2(:,j) = eta * grad_hide2(j,1) * yh1(:,1)
					wh2(:,j) = wh2(:,j) + alpha * (wh2(:,j) - wh2_velho(:,j)) + deltaw_h2(:,j)      
					deltabh2(j,1) = eta * grad_hide2(j,1) * (-1.d0)
					bh2(j,1) = bh2(j,1) + deltabh2(j,1)

				ENDDO !2. CAMADA OCULTA

				!TREINAMENTO: 1. CAMADA OCULTA
				DO j = 1,hide_neuron(1,1)

					soma = 0.d0
					do k = 1,hide_neuron(2,1)
						soma = soma + ( grad_hide2(k,1) * wh2_velho(j,k) )
					enddo

					select case(f_deriva)
						case (1)!LOGISTICA: e^-x / (1+e^-x)^2
							dv = ((a*DEXP(-a * vh1(j,1))) / ((1.d0 + DEXP(-a * vh1(j,1))) ** 2.d0))
						case (2)!TANGENTE
							dv = (2 * DEXP(-vh1(j,1))) / ((1 + DEXP(-vh1(j,1))) ** 2)
						case (3)!GAUSS
							dv = -yh1(j,1) / a
					end select

					grad_hide1(j,1) = dv * soma
   					deltaw_h1(:,j) = eta * grad_hide1(j,1) * x(:,i)
					wh1(:,j) = wh1(:,j) + alpha * (wh1(:,j) - wh1_velho(:,j)) + deltaw_h1(:,j)      
					deltabh1(j,1) = eta * grad_hide1(j,1) * (-1.d0)
  					bh1(j,1) = bh1(j,1) + deltabh1(j,1)

				ENDDO !1. CAMADA OCULTA

               	!CALCULO PADRAO DO ERRO!!!!
					erro_pad(i,1) = sum(erro(:,i),dim=1)
					erro_pad(i,1) = 0.5d0*(erro_pad(i,1)**2.d0)
							
			ENDDO !NUMERO PADROES

			eqm(1,l) = SUM(erro_pad(:,1))
			eqm(1,l) = (1.d0 / (num_pad)) * eqm(1,l)
			
			
			!***********************************************************************
			! VALIDACAO CRUZADA	
			!***********************************************************************

			do i = 1,num_pad_valid

				! ATIVACAO: 1. CAMADA OCULTA
				vh1 = 0.d0
 				vh1(:,1) = matmul(x_valid(:,i), wh1(:,:))
				vh1(:,1) = vh1(:,1) - bh1(:,1);
				select case(f_ativa)
					case (1) !LOGISTICA
						yh1(:,1) = 1.d0 / (1.d0 + DEXP(-a * vh1(:,1)))
					case (2) !TANGENTE
						yh1(:,1) = (1.d0 - DEXP(-vh1(:,1))) / (1.d0 + DEXP(-vh1(:,1)))
					case (3) !GAUSS
						yh1(:,1) = DEXP(-vh1(:,1)**2)
				end select

				! ATIVACAO: 2. CAMADA OCULTA
				vh2(:,1) = 0.0
				vh2(:,1) = MATMUL(yh1(:,1), wh2(:,:))
				vh2(:,1) = vh2(:,1) - bh2(:,1)
					
				SELECT CASE(f_ativa)
					CASE (1) !LOGISTICA
						yh2(:,1) = 1.d0/(1.d0+DEXP(-a*vh2(:,1)))
					CASE (2) !TANGENTE
						yh2(:,1) = (1.d0-DEXP(-vh2(:,1)))/(1.d0+DEXP(-vh2(:,1)))
					CASE (3) !GAUSS
						yh2(:,1) = DEXP(-vh2(:,1)**2)
				END SELECT
					
				! ATIVACAO: CAMADA DE SAIDA
				vs = 0.d0
				vs(:,1) = MATMUL(yh2(:,1), ws(:,:))
				vs(:,1) = vs(:,1) - bs(:,1)

				select case(f_ativa)
					CASE (1) !LOGISTICA
						ys(:,i) = 1.d0 / (1.d0 + DEXP(-a * vs(:,1)))
					CASE (2) !TANGENTE
						ys(:,i) = (1.d0 - DEXP(-vs(:,1))) / (1.d0 + DEXP(-vs(:,1)))
					CASE (3) !GAUSS
						ys(:,i) = DEXP(-vs(:,1)**2)
				END SELECT

				!CALCULO DO ERRO RNA
				erro(:,i) = yd_valid(:,i) - ys(:,i)
				
				!CALCULO PADRAO DO ERRO!!!!
				erro_pad(i,1) = sum(erro(:,i),dim=1)
				erro_pad(i,1) = 0.5d0*(erro_pad(i,1)**2.d0)

			ENDDO !DO VALIDACAO

			!CALCULO DOS ERROS
			eqm_valid(1,l) = sum(erro_pad(:,1))
			eqm_valid(1,l) = (1.d0/(num_pad_valid))*eqm_valid(1,l)

			if (eqm_valid(1,l) < erro_menor) then
				erro_menor = eqm_valid(1,l)
				wh1_menor = wh1
				wh2_menor = wh2
				bh1_menor = bh1
				bh2_menor = bh2
				ys_melhor = ys
				ws_menor = ws
				bs_menor = bs
				eqm_menor = eqm
				eqm_valid_menor = eqm_valid
			endif

			if (cont >= 100) then
				if (mpca .EQ. 0) then
					WRITE(*,*)l,eqm(1,l), eqm_valid(1,l)
				endif
				cont = 1
				eta = eta * 0.99
			else
				cont = cont+1
			endif
				
		ENDDO !DO MAXIMO ITERAÇÃO

		IF ( mpca .EQ. 0 ) THEN
			Open(12, FILE = path_result)!,status="REPLACE", ACCESS="SEQUENTIAL",ACTION="write")
				do i=1,vetor_saida
					write(12,110)(ys(i,j),j=1,num_pad_valid)
				enddo
			close(12)

			Open(12, FILE = path_wh1)!,status="REPLACE", ACCESS="SEQUENTIAL",ACTION="write")
				do i = 1, vetor_entrada
					write(12,110)(wh1_menor(i,j), j = 1, hide_neuron(1,1))
				enddo
			close(12)

			Open(12, FILE = path_bh1)!,status="REPLACE", ACCESS="SEQUENTIAL",ACTION="write")
				write(12,110)(bh1_menor(j,1), j = 1, hide_neuron(1,1))
			close(12)

			Open(12, FILE = path_wh2)!,status="REPLACE", ACCESS="SEQUENTIAL",ACTION="write")
				do i = 1, hide_neuron(1,1)
					write(12,110)(wh2_menor(i,j), j = 1, hide_neuron(2,1))
				enddo
			close(12)

			Open(12, FILE = path_bh2)!,status="REPLACE", ACCESS="SEQUENTIAL",ACTION="write")
				write(12,110)(bh2_menor(j,1), j = 1, hide_neuron(2,1))
			close(12)

			Open(12, FILE = path_ws)!,status="REPLACE", ACCESS="SEQUENTIAL",ACTION="write")
				do i = 1, hide_neuron(hide_camada,1)
					write(12,110)(ws_menor(i,j), j = 1, vetor_saida)
				enddo
			close(12)

			Open(12, FILE = path_bs)!,status="REPLACE", ACCESS="SEQUENTIAL",ACTION="write")
				write(12,110)(bs_menor(j,1), j = 1, vetor_saida)
			close(12)
110			format (424F32.16)
		ENDIF

		!CALCULO DO VALOR DA FUNCAO OBJETIVO DEFINIDA POR: Adenilson
                IF (validacao .EQ. 1) THEN
                        Rede_Neural_C2 = penaltyObj * ( (alphaObj * eqm(1,l) + betaObj * &
                        eqm_valid(1,l)) / (alphaObj + betaObj) )
                ELSE
                         Rede_Neural_C2 = penaltyObj * eqm(1,l)
                ENDIF

		deallocate(x,yd,x_valid,yd_valid,vh1,wh1,yh1)
		deallocate(wh1_velho, wh1_menor,bh1,bh1_velho,bh1_menor)
		deallocate(vs,ys,ys_melhor,ws,ws_velho,bs)
		deallocate(bs_velho,bs_menor,erro,erro_pad,grad_saida,grad_hide1)
		deallocate(deltabs,deltabh1,deltaw_saida,deltaw_h1,eqm)
		deallocate(eqm_valid,eqm_menor,eqm_valid_menor,grad_hide2)
		deallocate(vh2,wh2,wh2_velho,wh2_menor,bh2,bh2_velho,bh2_menor)
		deallocate(deltabh2,deltaw_h2,yh2)
		


	END FUNCTION Rede_Neural_C2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	DOUBLE PRECISION FUNCTION Rede_Neural_C3(parametros, path_ent, path_sai, path_ent_valid, path_sai_valid,&
		path_wh1, path_wh2, path_wh3, path_ws, path_bh1, path_bh2, path_bh3, path_bs, path_result)

		IMPLICIT NONE
	
		!Parametros recebidos pela funcao (passados pelo algoritmo de otimizacao)
		double precision, dimension(:,:) :: parametros
		character*52, path_ent,path_sai
		character*52, path_wh1, path_wh2, path_wh3, path_bh1, path_bh2, path_bh3, path_ws, path_bs
		character*52, path_sai_valid, path_ent_valid, path_result

		!Variaveis usadas para o calculo da funcao objetivo
		double precision :: penaltyObj, alphaObj, betaObj

		double precision, allocatable,dimension(:,:) :: bh1, bh1_velho, bh1_menor
		double precision, allocatable,dimension(:,:) :: bh2, bh2_velho, bh2_menor
		double precision, allocatable,dimension(:,:) :: bh3, bh3_velho, bh3_menor
		double precision, allocatable,dimension(:,:) :: bs, bs_velho, bs_menor
		double precision, allocatable,dimension(:,:) :: deltabs
		double precision, allocatable,dimension(:,:) :: deltabh1, deltabh2, deltabh3
		double precision, allocatable,dimension(:,:) :: deltaw_h1, deltaw_h2, deltaw_h3
		double precision, allocatable,dimension(:,:) :: deltaw_saida
		double precision, allocatable,dimension(:,:) :: eqm
		double precision, allocatable,dimension(:,:) :: eqm_menor, eqm_valid_menor
		double precision, allocatable,dimension(:,:) :: eqm_valid
		double precision, allocatable,dimension(:,:) :: erro
		double precision, allocatable,dimension(:,:) :: erro_pad
		double precision, allocatable,dimension(:,:) :: grad_saida
		double precision, allocatable,dimension(:,:) :: grad_hide1, grad_hide2, grad_hide3
		double precision, allocatable,dimension(:,:) :: vh1, vh2, vh3
		double precision, allocatable,dimension(:,:) :: vs, ys, ys_melhor
		double precision, allocatable,dimension(:,:) :: wh1, wh1_velho, wh1_menor
		double precision, allocatable,dimension(:,:) :: wh2, wh2_velho, wh2_menor
		double precision, allocatable,dimension(:,:) :: wh3, wh3_velho, wh3_menor
		double precision, allocatable,dimension(:,:) :: ws, ws_velho, ws_menor
		double precision, allocatable,dimension(:,:) :: x
		double precision, allocatable,dimension(:,:) :: x_valid
		double precision, allocatable,dimension(:,:) :: yd
		double precision, allocatable,dimension(:,:) :: yd_valid
		double precision, allocatable,dimension(:,:) :: yh1, yh2, yh3

		double precision, parameter :: a = 1

		!Parametros recebidos pela funcao (passados pelo algoritmo de otimizacao)
		double precision	:: entrada(16,1)
		double precision	:: dv,soma, erro_menor,alpha
		double precision    :: eta			! taxa de aprendizado
		double precision    :: ed			! erro desejado
		double precision    :: aleat        ! gerar pesos e bias aleatoriamente

		integer		:: hide_neuron(3,1) ! número de neurônios na camada escondida
		integer		:: cont, l, unit
		integer		:: max_it	     	! número máximo de iterações
		integer		:: vetor_entrada	! dimensão do vetor de entrada
		integer		:: num_pad			! número de padrões a serem treinados
		integer		:: num_pad_valid
		integer		:: vetor_saida	 	! dimensão do vetor de saída
		integer		:: i, j, k, hn		! variável para controle do programa
		integer		:: f_ativa			! controle de qual funcao ativacao
		integer		:: f_deriva			! controle de qual funcao derivada
		integer		:: hide_camada		! numero de camadas escondidas
		integer		:: mpca				! variavel para ativar o mpca (mpca = 0 --> sequencial, mpca = 1 --> paralelo)
		integer		:: peso				! variavel para controlar como os pesos serao inicializados (0=fixo, 1=aleat)
		integer		:: validacao		! variavel para controlar se terah ou nao validacao (validacao = 1 --> com validacao cruzada)

		!Atribuicao empirica dos pesos para a funcao objetivo
		penaltyObj = 1.0D0
		alphaObj = 1.0D0
		betaObj = 1.0D0

		DO I = 1, 16
			entrada(I,1) = parametros(I,1)
		END DO

		hide_camada = entrada(1,1)
		hide_neuron(1,1) = entrada(2,1)	!numero de neuronios na camada escondida
		hide_neuron(2,1) = entrada(3,1)	!numero de neuronios na camada escondida
		hide_neuron(3,1) = entrada(4,1)	!numero de neuronios na camada escondida
		f_ativa = entrada(5,1)			!parametro que corresponde ao tipo de funcao que sera utilizada
		alpha = entrada(6,1)			!termo momentum
		eta = entrada(7,1)		
		num_pad = entrada(8,1)
		num_pad_valid = entrada(9,1)
		vetor_entrada = entrada(10,1)	!numero de entradas
		vetor_saida = entrada(11,1)		!numero de saidas		
		ed = entrada(12,1)				!erro desejado
		max_it = entrada(13,1)			!numero máximo de iteração
		mpca = entrada(14,1)			!mpca = 0 --> serial, mpca = 1 --> paralelo
		peso = entrada(15,1)			!peso = 0 --> valor fixo, peso = 1 --> valor aleatorio
		validacao = entrada(16,1)		!validacao = 0 --> sem validacao cruzada, validacao = 1 --> com validacao cruzada
		
		f_deriva = f_ativa

		!inicializando o menor erro
		erro_menor = 10**9

		!BIAS 1. CAMADA ESCONDIDA
		allocate(bh1(hide_neuron(1,1),1))
		allocate(bh1_menor(hide_neuron(1,1),1))
		allocate(bh1_velho(hide_neuron(1,1),1))
		!BIAS 2. CAMADA ESCONDIDA
		allocate(bh2(hide_neuron(2,1),1))
		allocate(bh2_menor(hide_neuron(2,1),1))
		allocate(bh2_velho(hide_neuron(2,1),1))
		!BIAS 3. CAMADA ESCONDIDA
		allocate(bh3(hide_neuron(3,1),1))
		allocate(bh3_menor(hide_neuron(3,1),1))
		allocate(bh3_velho(hide_neuron(3,1),1))
		!BIAS CAMADA DE SAIDA
		allocate(bs(vetor_saida,1))
		allocate(bs_menor(vetor_saida,1))
		allocate(bs_velho(vetor_saida,1))
		!DELTABIAS
		allocate(deltabh1(hide_neuron(1,1),1))
		allocate(deltabh2(hide_neuron(2,1),1))
		allocate(deltabh3(hide_neuron(3,1),1))
		allocate(deltabs(vetor_saida,1))
		!DELTAW'S
		allocate(deltaw_h1(vetor_entrada,hide_neuron(1,1)))
		allocate(deltaw_h2(hide_neuron(1,1),hide_neuron(2,1)))
		allocate(deltaw_h3(hide_neuron(2,1),hide_neuron(3,1)))
		allocate(deltaw_saida(hide_neuron(hide_camada,1),vetor_saida))
		!ERROS
		allocate(eqm(1,max_it))
		allocate(eqm_valid(1,max_it))
		allocate(eqm_valid_menor(1,max_it))
		allocate(eqm_menor(1,max_it))
		allocate(erro(vetor_saida,num_pad))
		allocate(erro_pad(num_pad,1))
		!GRADIENTES
		allocate(grad_hide1(hide_neuron(1,1),1))
		allocate(grad_hide2(hide_neuron(2,1),1))
		allocate(grad_hide3(hide_neuron(3,1),1))
		allocate(grad_saida(vetor_saida,1))
		!CAMPO LOCAL INDUZIDO
		allocate(vh1(hide_neuron(1,1),1))
		allocate(vh2(hide_neuron(2,1),1))
		allocate(vh3(hide_neuron(3,1),1))
		allocate(vs(vetor_saida,1))
		!PESOS 1. CAMADA ESCONDIDA
		allocate(wh1(vetor_entrada,hide_neuron(1,1)))
		allocate(wh1_menor(vetor_entrada,hide_neuron(1,1)))
		allocate(wh1_velho(vetor_entrada,hide_neuron(1,1)))
		!PESOS 2. CAMADA ESCONDIDA
		allocate(wh2(hide_neuron(1,1),hide_neuron(2,1)))
		allocate(wh2_menor(hide_neuron(1,1),hide_neuron(2,1)))
		allocate(wh2_velho(hide_neuron(1,1),hide_neuron(2,1)))
		!PESOS 3. CAMADA ESCONDIDA
		allocate(wh3(hide_neuron(2,1),hide_neuron(3,1)))
		allocate(wh3_menor(hide_neuron(2,1),hide_neuron(3,1)))
		allocate(wh3_velho(hide_neuron(2,1),hide_neuron(3,1)))
		!PESOS CAMADA DE SAIDA
		allocate(ws(hide_neuron(hide_camada,1),vetor_saida))
		allocate(ws_menor(hide_neuron(hide_camada,1),vetor_saida))
		allocate(ws_velho(hide_neuron(hide_camada,1),vetor_saida))
		!DADOS DE ENTRADA E SAIDA
		allocate(x(vetor_entrada,num_pad))
		allocate(x_valid(vetor_entrada,num_pad_valid))
		allocate(yd(vetor_saida,num_pad))
		allocate(yd_valid(vetor_saida,num_pad_valid))
		!SAIDAS OBTIDAS
		allocate(yh1(hide_neuron(1,1),1))
		allocate(yh2(hide_neuron(2,1),1))
		allocate(yh3(hide_neuron(3,1),1))
		allocate(ys(vetor_saida,num_pad))
		allocate(ys_melhor(vetor_saida,num_pad)) 
	
		!------------------------------------------------------------!
		!LENDO OS PARAMETROS DO ARQUIVO DE ENTRADA
		!------------------------------------------------------------!
		OPEN (1, file = path_sai)
		DO I = 1, vetor_saida
			READ(1,*) (yd(I,J), J = 1, num_pad)
		END DO
		CLOSE (1)

		OPEN (2, file = path_ent)
		DO I = 1, vetor_entrada
			READ(2,*) (x(I,J), J = 1, num_pad)
		END DO
		CLOSE (2)

		IF (validacao .EQ. 1) THEN
			OPEN (1, file = path_sai_valid)
			DO I = 1, vetor_saida
				READ(1,*) (yd_valid(I,J), J = 1, num_pad_valid)
			END DO
			CLOSE (1)

			OPEN (2, file = path_ent_valid)
			DO I = 1, vetor_entrada
				READ(2,*) (x_valid(I,J), J=1, num_pad_valid)
			END DO
			CLOSE (2)
		ENDIF

		!--------------------------------------------------------------------!
		!INICIALIZANDO OS PARAMETROS: wh1, bh1, ws e bs
		!--------------------------------------------------------------------!
	
		!PRIMEIRA CAMADA OCULTA
        	DO l = 1,vetor_entrada
			DO k = 1, hide_neuron(1,1)
				if (peso .EQ. 0) then
					wh1(l,k) = 0.5
				else 
					call random_number(aleat)
					wh1(l,k) = aleat
				endif
            		ENDDO
        	ENDDO
		DO k = 1, hide_neuron(1,1)
			if (peso .EQ. 0) then
				bh1(k,1) = 0.5
			else
				call random_number(aleat)
				bh1(k,1) = aleat
			endif
		ENDDO

		!SEGUNDA CAMADA OCULTA
		DO k = 1, hide_neuron(1,1)
			DO i = 1,hide_neuron(2,1)
				if (peso .EQ. 0) then
					wh2(k,i) = 0.5
				else
					call random_number(aleat)
					wh2(k,i) = aleat
				endif
			ENDDO
		ENDDO
		DO i = 1,hide_neuron(2,1)
			if (peso .EQ. 0) then
				bh2(i,1) = 0.5
			else
                call random_number(aleat)
                bh2(i,1) = aleat
			endif
		ENDDO

		!INICIALIZACAO: 3.CAMADA OCULTA
		DO k = 1, hide_neuron(2,1)
			DO i = 1,hide_neuron(3,1)
				if (peso .EQ. 0) then
					wh3(k,i) = 0.5
				else
					call random_number(aleat)
					wh3(k,i) = aleat
				endif
			ENDDO
		ENDDO
		DO i = 1,hide_neuron(3,1)
			if (peso .EQ. 0) then
				bh3(i,1) = 0.5
			else
                	call random_number(aleat)
                	bh3(i,1) = aleat
			endif
		ENDDO

		!CAMADA DE SAIDA
		DO l = 1,hide_neuron(hide_camada,1)
			DO k = 1, vetor_saida
				if (peso .EQ. 0) then
					ws(l,k) = 0.5
				else
	               	call random_number(aleat)
		        	ws(l,k) = aleat
				endif
			ENDDO
		ENDDO
		DO k = 1,vetor_saida
			if (peso .EQ. 0) then
				bs(k,1) = 0.5
			else
				call random_number(aleat)
				bs(k,1) = aleat
			endif
		ENDDO

		!----------------------------------------------------------------------!
		! INICIO DA REDE: FEEDFORWARD
		!----------------------------------------------------------------------!
		cont = 1 !contador do numero de epocas usado para alterar o valor do ETA
		l = 0 
	
		DO WHILE ((erro_menor .GE. ed) .AND. (l .LT. max_it))

			l = l + 1
	
			DO i = 1, num_pad	! DO NUMERO TOTAL DE PADROES
		
				IF (hide_camada .EQ. 3) THEN
					! ATIVACAO: 1. CAMADA OCULTA
					vh1 = 0.d0
					vh1(:,1) = matmul(x(:,i),wh1(:,:))
					vh1(:,1) = vh1(:,1) - bh1(:,1);
					SELECT CASE(f_ativa)
						CASE (1) ! LOGISTICA
							yh1(:,1) = 1.d0 / (1.d0 + DEXP(-a * vh1(:,1)))
						CASE (2) ! TANGENTE
							yh1(:,1) = (1.d0 - DEXP(-vh1(:,1))) / (1.d0 + DEXP(-vh1(:,1)))
						CASE (3) ! GAUSS
							yh1(:,1) = DEXP(-vh1(:,1)**2)
					END SELECT
			
					! ATIVACAO: 2. CAMADA OCULTA
					vh2(:,1) = 0.0
					vh2(:,1) = matmul(yh1(:,1),wh2(:,:))
					vh2(:,1) = vh2(:,1) - bh2(:,1)
					SELECT CASE(f_ativa)
						CASE (1) ! LOGISTICA
							yh2(:,1) = 1.d0 / (1.d0 + DEXP(-a * vh2(:,1)))
						CASE (2) ! TANGENTE
							yh2(:,1) = (1.d0 - DEXP(-vh2(:,1))) / (1.d0 + DEXP(-vh2(:,1)))
						CASE (3) ! GAUSS
							yh2(:,1) = DEXP(-vh2(:,1)**2)
					END SELECT

					! ATIVACAO: 3. CAMADA OCULTA
					vh3(:,1) = 0.0
					vh3(:,1) = matmul(yh2(:,1),wh3(:,:))
					vh3(:,1) = vh3(:,1) - bh3(:,1)
					SELECT CASE(f_ativa)
						CASE (1) ! LOGISTICA
							yh3(:,1) = 1.d0 / (1.d0 + DEXP(-a * vh3(:,1)))
						CASE (2) ! TANGENTE
							yh3(:,1) = (1.d0 - DEXP(-vh3(:,1))) / (1.d0 + DEXP(-vh3(:,1)))
						CASE (3) ! GAUSS
							yh3(:,1) = DEXP(-vh3(:,1))
					END SELECT
			
					! ATIVACAO: CAMADA DE SAIDA
					vs = 0.d0
					vs(:,1) = matmul(yh3(:,1),ws(:,:))
					vs(:,1) = vs(:,1) - bs(:,1)
					SELECT CASE(f_ativa)
						CASE (1) ! LOGISTICA
							ys(:,i) = 1.d0 / (1.d0 + DEXP(-a * vs(:,1)))
						CASE (2) ! TANGENTE
							ys(:,i) = (1.d0 - DEXP(-vs(:,1))) / (1.d0 + DEXP(-vs(:,1)))
						CASE (3) ! GAUSS
							ys(:,i) = DEXP(-vs(:,1)**2)
					END SELECT

				ENDIF ! HIDE_CAMADA .EQ. 3
	
				! PARA CALCULO DO NOVO PESO
				wh1_velho = wh1
				bh1_velho = bh1
				wh2_velho = wh2
				bh2_velho = bh2
				wh3_velho = wh3
				bh3_velho = bh3
				ws_velho = ws
				bs_velho = bs
			
				! CALCULO ERRO TREINAMENTO
				erro(:,i) = yd(:,i) - ys(:,i)

		!-------------------------------------------------------------------------!
		!                        BACKPROPAGATION
		!-------------------------------------------------------------------------!

				! TREINAMENTO: CAMADA DE SAIDA
				DO j = 1,vetor_saida
	
					SELECT CASE(f_deriva)
						CASE (1)! LOGISTICA: e^-x / (1+e^-x)^2
							dv = ((a * DEXP(-a * vs(j,1))) / ((1.d0 + DEXP(-a * vs(j,1))) ** 2.d0))
						CASE (2)! TANGENTE: 2e^-x / (1+e^-x)^2
							dv = (2 * DEXP(-vs(j,1))) / ((1 + DEXP(-vs(j,1))) ** 2)
						CASE (3)! GAUSS
							dv = -ys(j,i) / a	
					END SELECT
	
					grad_saida(j,1) = erro(j,i) * dv
					deltaw_saida(:,j) = eta * grad_saida(j,1) * yh3(:,1)
					ws(:,j) = ws(:,j) + alpha * (ws(:,j) - ws_velho(:,j)) + deltaw_saida(:,j)
			   		deltabs(j,1) = eta * grad_saida(j,1) * (-1.d0)
			   		bs(j,1) = bs(j,1) + deltabs(j,1)

				ENDDO ! CAMADA SAIDA

				! TREINAMENTO: 3. CAMADA OCULTA
				DO j = 1,hide_neuron(3,1)
						soma = 0.d0
					DO k = 1, vetor_saida
						soma = soma + ( grad_saida(k,1) * ws_velho(j,k) )
					ENDDO

					SELECT CASE(f_deriva)
						CASE (1)! LOGISTICA: e^-x / (1+e^-x)^2
							dv = ((a * DEXP(-a * vh3(j,1))) / ((1.d0 + DEXP(-a * vh3(j,1))) ** 2.d0))
						CASE (2)! TANGENTE
							dv = (2 * DEXP(-vh3(j,1))) / ((1 + DEXP(-vh3(j,1))) ** 2)
						CASE (3)! GAUSS
							dv = -yh3(j,1) / a
					END SELECT
	
					grad_hide3(j,1) = dv * soma
					deltaw_h3(:,j) = eta * grad_hide3(j,1) * yh2(:,1)
					wh3(:,j) = wh3(:,j) + alpha * (wh3(:,j) - wh3_velho(:,j)) + deltaw_h3(:,j)      
					deltabh3(j,1) = eta * grad_hide3(j,1) * (-1.d0)
					bh3(j,1) = bh3(j,1) + deltabh3(j,1)

				ENDDO ! 3. CAMADA OCULTA

				! TREINAMENTO: 2. CAMADA OCULTA
				DO j = 1,hide_neuron(2,1)
						soma = 0.d0
					DO k = 1, hide_neuron(3,1)
						soma = soma + ( grad_hide3(k,1) * wh3_velho(j,k) )
					ENDDO

					SELECT CASE(f_deriva)
						CASE (1)!LOGISTICA: e^-x / (1+e^-x)^2
							dv = ((a * DEXP(-a * vh2(j,1))) / ((1.d0 + DEXP(-a * vh2(j,1))) ** 2.d0))
						CASE (2)!TANGENTE
							dv = (2 * DEXP(-vh2(j,1)))/ ((1 + DEXP(-vh2(j,1))) ** 2)
						CASE (3)!GAUSS
							dv = -yh2(j,1) / a
					END SELECT
	
					grad_hide2(j,1) = dv * soma
					deltaw_h2(:,j) = eta * grad_hide2(j,1) * yh1(:,1)
					wh2(:,j) = wh2(:,j) + alpha * (wh2(:,j) - wh2_velho(:,j)) + deltaw_h2(:,j)      
					deltabh2(j,1) = eta * grad_hide2(j,1) * (-1.d0)
					bh2(j,1) = bh2(j,1) + deltabh2(j,1)

				ENDDO ! 2. CAMADA OCULTA

				! TREINAMENTO: 1. CAMADA OCULTA
				DO j = 1,hide_neuron(1,1)
					soma = 0.d0
					DO k = 1,hide_neuron(2,1)
						soma = soma + ( grad_hide2(k,1) * wh2_velho(j,k) )
					ENDDO

					SELECT CASE(f_deriva)
						CASE (1)!LOGISTICA: e^-x / (1+e^-x)^2
							dv = ((a * DEXP(-a * vh1(j,1))) / ((1.d0 + DEXP(-a * vh1(j,1))) ** 2.d0))
						CASE (2)!TANGENTE
							dv = (2 * DEXP(-vh1(j,1))) / ((1 + DEXP(-vh1(j,1))) ** 2)
						CASE (3)!GAUSS
							dv = -yh1(j,1) / a	
					END SELECT

					grad_hide1(j,1) = dv * soma
		   			deltaw_h1(:,j) = eta * grad_hide1(j,1) * x(:,i)
					wh1(:,j) = wh1(:,j) + alpha * (wh1(:,j) - wh1_velho(:,j)) + deltaw_h1(:,j)      
					deltabh1(j,1) = eta * grad_hide1(j,1) * (-1.d0)
		  			bh1(j,1) = bh1(j,1) + deltabh1(j,1)
				
				ENDDO ! 1. CAMADA OCULTA

               	!CALCULO PADRAO DO ERRO!!!!
					erro_pad(i,1) = sum(erro(:,i),dim=1)
					erro_pad(i,1) = 0.5d0*(erro_pad(i,1)**2.d0)
			
			ENDDO ! NUMERO PADROES

			eqm(1,l) = sum(erro_pad(:,1))
			eqm(1,l) = (1.d0/(num_pad))*eqm(1,l)
	
			!***********************************************************************
			! VALIDACAO CRUZADA	
			!***********************************************************************

			DO i = 1,num_pad_valid

				IF (hide_camada .EQ. 3) THEN
	
					! ATIVACAO: 1. CAMADA OCULTA
					vh1 = 0.d0
			 		vh1(:,1) = matmul(x(:,i),wh1(:,:))
					vh1(:,1) = vh1(:,1) - bh1(:,1);
					SELECT CASE(f_ativa)
						CASE (1) !LOGISTICA
							yh1(:,1) = 1.d0 / (1.d0 + DEXP(-a * vh1(:,1)))
						CASE (2) !TANGENTE
							yh1(:,1) = (1.d0 - DEXP(-vh1(:,1))) / (1.d0 + DEXP(-vh1(:,1)))
						CASE (3) !GAUSS
							yh1(:,1) = DEXP(-vh1(:,1)**2)
					END SELECT

					! ATIVACAO: 2. CAMADA OCULTA
					vh2(:,1) = 0.0
					vh2(:,1) = matmul(yh1(:,1),wh2(:,:))
					vh2(:,1) = vh2(:,1) - bh2(:,1)
					SELECT CASE(f_ativa)
						CASE (1) !LOGISTICA
							yh2(:,1) = 1.d0 / (1.d0 + DEXP(-a * vh2(:,1)))
						CASE (2) !TANGENTE
							yh2(:,1) = (1.d0 - DEXP(-vh2(:,1))) / (1.d0 + DEXP(-vh2(:,1)))
						CASE (3) !GAUSS
							yh2(:,1) = DEXP(-vh2(:,1)**2)
					END SELECT

					! ATIVACAO: 3. CAMADA OCULTA
					vh3(:,1) = 0.0
					vh3(:,1) = matmul(yh2(:,1),wh3(:,:))
					vh3(:,1) = vh3(:,1) - bh3(:,1)
					SELECT CASE(f_ativa)
						CASE (1) !LOGISTICA
							yh3(:,1) = 1.d0 / (1.d0 + DEXP(-a * vh3(:,1)))
						CASE (2) !TANGENTE
							yh3(:,1) = (1.d0 - DEXP(-vh3(:,1))) / (1.d0 + DEXP(-vh3(:,1)))
						CASE (3) !GAUSS
							yh3(:,1) = DEXP(-vh3(:,1)**2)
					END SELECT
		
					! ATIVACAO: CAMADA DE SAIDA
					vs = 0.d0
					vs(:,1) = MATMUL(yh3(:,1),ws(:,:))
					vs(:,1) = vs(:,1) - bs(:,1)
					SELECT CASE(f_ativa)
						CASE (1) !LOGISTICA
							ys(:,i) = 1.d0 / (1.d0 + DEXP(-a * vs(:,1)))
						CASE (2) !TANGENTE
							ys(:,i) = (1.d0 - DEXP(-vs(:,1))) / (1.d0 + DEXP(-vs(:,1)))
						CASE (3) !GAUSS
							ys(:,i) = DEXP(-vs(:,1)**2)
					END SELECT

				ENDIF !IF hide_camada .EQ. 3

				!CALCULO DO ERRO RNA
				erro(:,i) = yd_valid(:,i) - ys(:,i)
				
				!CALCULO PADRAO DO ERRO!!!!
				erro_pad(i,1) = sum(erro(:,i),dim=1)
				erro_pad(i,1) = 0.5d0*(erro_pad(i,1)**2.d0)
				
			ENDDO !DO VALIDACAO

			!CALCULO DOS ERROS
			eqm_valid(1,l) = sum(erro_pad(:,1))
			eqm_valid(1,l) = (1.d0/(num_pad_valid))*eqm_valid(1,l)

			if (eqm_valid(1,l) < erro_menor) then
				erro_menor = eqm_valid(1,l)
				wh1_menor = wh1
				bh1_menor = bh1
				wh2_menor = wh2
				bh2_menor = bh2
				wh3_menor = wh3
				bh3_menor = bh3
				ys_melhor = ys
				ws_menor = ws
				bs_menor = bs
				eqm_menor = eqm
				eqm_valid_menor = eqm_valid
			endif

			if (cont >= 100) then
				if (mpca .EQ. 0) then
					WRITE(*,*)l,eqm(1,l), eqm_valid(1,l)
				endif
				cont = 1
				eta = eta * 0.99
			else
				cont = cont+1
			endif
				
		ENDDO !DO MAXIMO ITERAÇÃO

		IF ( mpca .EQ. 0 ) THEN
			Open(12, FILE = path_result)!,status="REPLACE", ACCESS="SEQUENTIAL",ACTION="write")
				do i=1,vetor_saida
					write(12,110)(ys(i,j),j=1,num_pad_valid)
				enddo
			close(12)
			Open(12, FILE = path_wh1)!,status="REPLACE", ACCESS="SEQUENTIAL",ACTION="write")
				do i = 1, vetor_entrada
					write(12,110)(wh1_menor(i,j), j = 1, hide_neuron(1,1))
				enddo
			close(12)

			Open(12, FILE = path_bh1)!,status="REPLACE", ACCESS="SEQUENTIAL",ACTION="write")
				write(12,110)(bh1_menor(j,1), j = 1, hide_neuron(1,1))
			close(12)

			Open(12, FILE = path_wh2)!,status="REPLACE", ACCESS="SEQUENTIAL",ACTION="write")
				do i = 1, hide_neuron(1,1)
					write(12,110)(wh2_menor(i,j), j = 1, hide_neuron(2,1))
				enddo
			close(12)

			Open(12, FILE = path_bh2)!,status="REPLACE", ACCESS="SEQUENTIAL",ACTION="write")
				write(12,110)(bh2_menor(j,1), j = 1, hide_neuron(2,1))
			close(12)

			Open(12, FILE = path_wh3)!,status="REPLACE", ACCESS="SEQUENTIAL",ACTION="write")
				do i = 1, hide_neuron(2,1)
					write(12,110)(wh3_menor(i,j), j = 1, hide_neuron(3,1))
				enddo
			close(12)

			Open(12, FILE = path_bh3)!,status="REPLACE", ACCESS="SEQUENTIAL",ACTION="write")
				write(12,110)(bh3_menor(j,1), j = 1, hide_neuron(3,1))
			close(12)

			Open(12, FILE = path_ws)!,status="REPLACE", ACCESS="SEQUENTIAL",ACTION="write")
				do i = 1, hide_neuron(hide_camada,1)
					write(12,110)(ws_menor(i,j), j = 1, vetor_saida)
				enddo
			close(12)

			Open(12, FILE = path_bs)!,status="REPLACE", ACCESS="SEQUENTIAL",ACTION="write")
				write(12,110)(bs_menor(j,1), j = 1, vetor_saida)
			close(12)
110			format (424F32.16)
		ENDIF

		!CALCULO DO VALOR DA FUNCAO OBJETIVO DEFINIDA POR: Adenilson
                IF (validacao .EQ. 1) THEN
                        Rede_Neural_C3 = penaltyObj * ( (alphaObj * eqm(1,l) + betaObj * &
                        eqm_valid(1,l)) / (alphaObj + betaObj) )
                ELSE
                         Rede_Neural_C3 = penaltyObj * eqm(1,l)
                ENDIF

		deallocate(x,yd,x_valid,yd_valid,vh1,wh1,yh1)
		deallocate(wh1_velho, wh1_menor,bh1,bh1_velho,bh1_menor)
		deallocate(vs,ys,ys_melhor,ws,ws_velho,bs)
		deallocate(bs_velho,bs_menor,erro,erro_pad,grad_saida,grad_hide1)
		deallocate(deltabs,deltabh1,deltaw_saida,deltaw_h1,eqm)
		deallocate(eqm_valid,eqm_menor,eqm_valid_menor,grad_hide2)
		deallocate(vh2,wh2,wh2_velho,wh2_menor,bh2,bh2_velho,bh2_menor)
		deallocate(deltabh2,deltaw_h2,yh2,yh3)
		deallocate(vh3,wh3,wh3_velho,wh3_menor,bh3,bh3_velho,bh3_menor)
		deallocate(grad_hide3,deltabh3,deltaw_h3)	

	!Fim da funcao
	END FUNCTION Rede_Neural_C3

!Fim do modulo
END MODULE ModuloRNA
