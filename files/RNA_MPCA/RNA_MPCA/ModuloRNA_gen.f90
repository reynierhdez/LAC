!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RNA: PERCEPTRON DE MULTIPLAS CAMADAS (MLP)- AUTOCONFIGURAVEL 
! USANDO MULTIPLE PARTICLE COLLISION ALGORITHM (MPCA)
! 
! Desenvolvido por: Juliana A Anochi (CAP/INPE)
!					Sabrina Sambatti (CAP/INPE)
! Email: juliana.anochi@lac.inpe.br
!        sabrinabms@gmail.com          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE ModuloRNA_gen
	CONTAINS

	!Funcao para treinar a rede neural com os parametros calculados pela heuristica
	DOUBLE PRECISION FUNCTION  Rede_Neural_BP_gen(parametros)
		
		IMPLICIT NONE

		DOUBLE PRECISION :: parametros(8,1)
		CHARACTER*24, caminho
		CHARACTER*24, arq_ent, arq_sai, arq_result_ys
		CHARACTER*24, arq_wh1, arq_wh2, arq_wh3, arq_ws
		CHARACTER*24, arq_bh1, arq_bh2, arq_bh3, arq_bs
		CHARACTER*52, path_ent, path_sai
		CHARACTER*52, path_wh1, path_wh2, path_wh3, path_ws
		CHARACTER*52, path_bh1, path_bh2, path_bh3, path_bs
		CHARACTER*52, path_result_ys
		INTEGER :: hide_camada
		INTEGER :: i

		hide_camada = parametros(1,1)

		!Definindo os nomes dos arquivos de entrada, saida e resultados
		arq_ent ='x_gen.txt'
		arq_sai ='yd_gen.txt'
		arq_wh1 = 'wh1.txt'
		arq_wh2 = 'wh2.txt'
		arq_wh3 = 'wh3.txt'
		arq_ws = 'ws.txt'
		arq_bh1 = 'bh1.txt'
		arq_bh2 = 'bh2.txt'
		arq_bh3 = 'bh3.txt'
		arq_bs = 'bs.txt'
		arq_result_ys =	'sgen_fixo_p0208.txt'
		
		!Definindo o caminho (path)
		caminho = './dados/ComRuido/'
		path_ent = TRIM(caminho)//TRIM(arq_ent)
		path_sai = TRIM(caminho)//TRIM(arq_sai)
		path_wh1 = TRIM(caminho)//TRIM(arq_wh1)
		path_wh2 = TRIM(caminho)//TRIM(arq_wh2)
		path_wh3 = TRIM(caminho)//TRIM(arq_wh3)
		path_ws = TRIM(caminho)//TRIM(arq_ws)
		path_bh1 = TRIM(caminho)//TRIM(arq_bh1)
		path_bh2 = TRIM(caminho)//TRIM(arq_bh2)
		path_bh3 = TRIM(caminho)//TRIM(arq_bh3)
		path_bs = TRIM(caminho)//TRIM(arq_bs)
		path_result_ys =TRIM(caminho)//TRIM(arq_result_ys)

		IF (hide_camada .EQ. 1) THEN
			Rede_Neural_BP_gen = Rede_Neural_C1_GEN(parametros, path_ent, path_sai, &
			path_wh1, path_ws, path_bh1, path_bs, path_result_ys)
		ELSE IF (hide_camada .EQ. 2) THEN
			Rede_Neural_BP_gen = Rede_Neural_C2_GEN(parametros, path_ent, path_sai, &
			path_wh1, path_wh2, path_ws, path_bh1, path_bh2, path_bs, path_result_ys)
		ELSE IF (hide_camada .EQ. 3) THEN
			Rede_Neural_BP_gen = Rede_Neural_C3_GEN(parametros, path_ent, path_sai, &
			path_bh1, path_bh2, path_bh3, path_bs, path_wh1, path_wh2, path_wh3, path_ws, path_result_ys)
		ENDIF

	END FUNCTION Rede_Neural_BP_gen

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	DOUBLE PRECISION FUNCTION Rede_Neural_C1_GEN(parametros,path_ent,path_sai, path_wh1, path_ws, &
		path_bh1, path_bs, path_result_ys)
	
		IMPLICIT NONE
	
		!Parametros recebidos pela funcao (passados pelo algoritmo de otimizacao)
		double precision, dimension(:,:) :: parametros
		CHARACTER*52, path_ent,path_sai
		CHARACTER*52, path_wh1, path_ws
		CHARACTER*52, path_bh1, path_bs
		CHARACTER*52, path_result_ys

		double precision, allocatable,dimension(:,:) :: x_gen
		double precision, allocatable,dimension(:,:) :: yd_gen
		double precision, allocatable,dimension(:,:) :: vs, vh1
		double precision, allocatable,dimension(:,:) :: ws, wh1
		double precision, allocatable,dimension(:,:) :: bs, bh1
		double precision, allocatable,dimension(:,:) :: ys, yh1
		double precision, allocatable,dimension(:,:) :: erro
		double precision, allocatable,dimension(:,:) :: emq
			
		double precision, parameter :: a = 1

		double precision	:: entrada(8,1)
		integer		:: hide_neuron(3,1) ! número de neurônios na camada escondida
		integer		:: vetor_entrada	! dimensão do vetor de entrada
		integer		:: num_pad_gen		! número de padrões a serem treinados
		integer		:: vetor_saida	 	! dimensão do vetor de saída
		integer		:: I, J			! variável para controle do programa
		integer		:: f_ativa			! controle de qual funcao ativacao
		integer		:: hide_camada		! numero de camadas escondidas

		DO I = 1, 8
			entrada(I,1) = parametros(I,1)
		END DO

		hide_camada = entrada(1,1)
		hide_neuron(1,1) = entrada(2,1)	!numero de neuronios na camada escondida
		hide_neuron(2,1) = entrada(3,1)	!numero de neuronios na camada escondida
		hide_neuron(3,1) = entrada(4,1)	!numero de neuronios na camada escondida
		f_ativa = entrada(5,1)			!parâmetro que corresponde ao tipo de funcao que será utilizada
		num_pad_gen = entrada(6,1)
		vetor_entrada = entrada(7,1)	!numero de entradas
		vetor_saida = entrada(8,1)		!numero de saidas	

		!Bias 1. camada escondida
		allocate(bh1(hide_neuron(1,1),1))
		!Bias camada de saida
		allocate(bs(vetor_saida,1))
		!Erros
		allocate(emq(1,1))
		allocate(erro(vetor_saida,num_pad_gen))
		!Campo local induzido
		allocate(vh1(hide_neuron(1,1),1))
		allocate(vs(vetor_saida,1))
		!Pesos 1. camada escondida
		allocate(wh1(vetor_entrada,hide_neuron(1,1)))
		!Pesos camada de saida
		allocate(ws(hide_neuron(hide_camada,1),vetor_saida))
		!Dados de entrada e saida
		allocate(x_gen(vetor_entrada,num_pad_gen))
		allocate(yd_gen(vetor_saida,num_pad_gen))
		!Saidas obtidas
		allocate(yh1(hide_neuron(1,1),1))
		allocate(ys(vetor_saida,num_pad_gen))

		!------------------------------------------------------------!
		!LENDO OS PARAMETROS DO ARQUIVO DE ENTRADA
		!------------------------------------------------------------!
		OPEN (1, file = path_sai)
			DO I = 1, vetor_saida
				READ(1,*) (yd_gen(I,J), J = 1, num_pad_gen)
			END DO
		CLOSE (1)

		OPEN (2, file = path_ent)
			DO I = 1, vetor_entrada
				READ(2,*) (x_gen(I,J), J=1, num_pad_gen)
			END DO
		CLOSE (2)

		!--------------------------------------------------------------------!
		!LENDO OS PARAMETROS: wh1, bh1, ws e bs
		!--------------------------------------------------------------------!
		OPEN (1, file = path_wh1)
			DO I = 1, vetor_entrada
				READ(1,*) (wh1(I,J),J=1,hide_neuron(1,1))
			END DO
		CLOSE (1)

		OPEN (1, file = path_bh1)
			READ(1,*) (bh1(J,1),J=1,hide_neuron(1,1))
		CLOSE (1)

		OPEN (1, file = path_ws)
			DO I = 1, hide_neuron(1,1)
				READ(1,*) (ws(I,J),J=1,vetor_saida)
			END DO
		CLOSE (1)

		OPEN (1, file = path_bs)
			READ(1,*) (bs(J,1),J=1,vetor_saida)
		CLOSE (1)

		!----------------------------------------------------------------------!
		! INICIO DA REDE: FEEDFORWARD
		!----------------------------------------------------------------------!
		DO I = 1, num_pad_gen	!DO NUMERO TOTAL DE PADROES GENERALIZACAO
				
			! Ativação das Camadas Ocultas
			vh1 = 0.d0
			vh1(:,1) = matmul(x_gen(:,I), wh1(:,:))
			vh1(:,1) = vh1(:,1) - bh1(:,1);
			select case(f_ativa)
				case (1) !LOGISTICA
					yh1(:,1) = 1.d0/(1.d0+DEXP(-a * vh1(:,1)))
				case (2) !TANGENTE
					yh1(:,1) = (1.d0-DEXP(-vh1(:,1)))/(1.d0+DEXP(-vh1(:,1)))
				case (3) !GAUSS
					yh1(:,1) = DEXP(-vh1(:,1))
			end select

			! Ativação da Camada de Saída
			!--------------------------------------------------------------!
			vs = 0.d0
			vs(:,1) = matmul(yh1(:,1),ws(:,:))
			vs(:,1) = vs(:,1) - bs(:,1)

			select case(f_ativa)
				case (1) !LOGISTICA
					ys(:,I) = 1.d0/(1.d0+DEXP(-a*vs(:,1)))
				case (2) !TANGENTE
					ys(:,I) = (1.d0-DEXP(-vs(:,1)))/(1.d0+DEXP(-vs(:,1)))
				case (3) !GAUSS
					ys(:,I) = DEXP(-vs(:,1))
			end select

			erro(:,I) = yd_gen(:,I) - ys(:,I)
					
		ENDDO !NUMERO PADROES

		OPEN(12, FILE = path_result_ys)
			DO I = 1, vetor_saida
				WRITE(12,110)(ys(I,J), J = 1, num_pad_gen)
			ENDDO
		CLOSE(12)
110		FORMAT (424F32.16)

		!CALCULO DO ERRO QUADRATICO MEDIO
		emq(1,1) = sum(erro(:,:))
		emq(1,1) = (1.d0/(num_pad_gen))*emq(1,1)
			
		Rede_Neural_C1_GEN = emq(1,1)
				
		deallocate(x_gen)
		deallocate(yd_gen)
		deallocate(vh1)
		deallocate(wh1)
		deallocate(bh1)
		deallocate(vs)
		deallocate(ys)
		deallocate(ws)
		deallocate(bs)
		deallocate(erro)
		deallocate(emq)
		deallocate(yh1)

	!FIM DA FUNCAO
	END FUNCTION Rede_Neural_C1_GEN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	DOUBLE PRECISION FUNCTION Rede_Neural_C2_GEN(parametros,path_ent,path_sai, path_wh1, path_wh2, &
		path_ws, path_bh1, path_bh2, path_bs, path_result_ys)

		IMPLICIT NONE

		!Parametros recebidos pela funcao (passados pelo algoritmo de otimizacao)
		double precision, dimension(:,:) :: parametros
		CHARACTER*52, path_ent,path_sai
		CHARACTER*52, path_wh1, path_wh2, path_ws
		CHARACTER*52, path_bh1, path_bh2, path_bs
		CHARACTER*52, path_result_ys

		double precision, allocatable,dimension(:,:) :: x_gen
		double precision, allocatable,dimension(:,:) :: yd_gen
		double precision, allocatable,dimension(:,:) :: vs, vh1, vh2
		double precision, allocatable,dimension(:,:) :: ws, wh1, wh2
		double precision, allocatable,dimension(:,:) :: bs, bh1, bh2
		double precision, allocatable,dimension(:,:) :: ys, yh1, yh2
		double precision, allocatable,dimension(:,:) :: erro
		double precision, allocatable,dimension(:,:) :: emq
			
		double precision, parameter :: a = 1

		double precision	:: entrada(8,1)

		integer		:: hide_neuron(3,1) ! número de neurônios na camada escondida
		integer		:: vetor_entrada	! dimensão do vetor de entrada
		integer		:: num_pad_gen		! número de padrões a serem treinados
		integer		:: vetor_saida	 	! dimensão do vetor de saída
		integer		:: I, J				! variável para controle do programa
		integer		:: f_ativa			! controle de qual funcao ativacao
		integer		:: hide_camada		! numero de camadas escondidas

		DO I = 1, 8
			entrada(I,1) = parametros(I,1)
		END DO

		hide_camada = entrada(1,1)
		hide_neuron(1,1) = entrada(2,1)!numero de neuronios na camada escondida
		hide_neuron(2,1) = entrada(3,1)!numero de neuronios na camada escondida
		hide_neuron(3,1) = entrada(4,1)!numero de neuronios na camada escondida
		f_ativa = entrada(5,1)!parâmetro que corresponde ao tipo de funcao que será utilizada
		num_pad_gen = entrada(6,1)
		vetor_entrada = entrada(7,1)	!numero de entradas
		vetor_saida = entrada(8,1)		!numero de saidas	

		!Bias 1. camada escondida
		allocate(bh1(hide_neuron(1,1),1))
		!Bias 2. camada escondida
		allocate(bh2(hide_neuron(2,1),1))
		!Bias camada de saida
		allocate(bs(vetor_saida,1))
		!Erros
		allocate(emq(1,1))
		allocate(erro(vetor_saida,num_pad_gen))
		!Campo local induzido
		allocate(vh1(hide_neuron(1,1),1))
		allocate(vh2(hide_neuron(2,1),1))
		allocate(vs(vetor_saida,1))
		!Pesos 1. camada escondida
		allocate(wh1(vetor_entrada,hide_neuron(1,1)))
		!Pesos 2. camada escondida
		allocate(wh2(hide_neuron(1,1),hide_neuron(2,1)))
		!Pesos camada de saida
		allocate(ws(hide_neuron(hide_camada,1),vetor_saida))
		!Dados de entrada e saida
		allocate(x_gen(vetor_entrada,num_pad_gen))
		allocate(yd_gen(vetor_saida,num_pad_gen))
		!Saidas obtidas
		allocate(yh1(hide_neuron(1,1),1))
		allocate(yh2(hide_neuron(2,1),1))
		allocate(ys(vetor_saida,num_pad_gen))

		!------------------------------------------------------------!
		!LENDO OS PARAMETROS DO ARQUIVO DE ENTRADA
		!------------------------------------------------------------!

		OPEN (1, file = path_sai)
			DO I = 1, vetor_saida
				READ(1,*) (yd_gen(I,J), J = 1, num_pad_gen)
			END DO
		CLOSE (1)

		OPEN (2, file = path_ent)
			DO I = 1, vetor_entrada
				READ(2,*) (x_gen(I,J), J=1, num_pad_gen)
			END DO
		CLOSE (2)

		!--------------------------------------------------------------------!
		!LENDO OS PARAMETROS: wh1, bh1, ws e bs
		!--------------------------------------------------------------------!
		OPEN (1, file = path_wh1)
			DO I = 1, vetor_entrada
				READ(1,*) (wh1(I,J),J=1,hide_neuron(1,1))
			END DO
		CLOSE (1)

		OPEN (1, file = path_wh2)
			DO I = 1, hide_neuron(1,1)
				READ(1,*) (wh2(I,J),J=1,hide_neuron(2,1))
			END DO
		CLOSE (1)

		OPEN (1, file = path_bh1)
			READ(1,*) (bh1(J,1),J=1,hide_neuron(1,1))
		CLOSE (1)

		OPEN (1, file = path_bh2)
			READ(1,*) (bh2(J,1),J=1,hide_neuron(2,1))
		CLOSE (1)

		OPEN (1, file = path_ws)
			DO I = 1, hide_neuron(2,1)
				READ(1,*) (ws(I,J),J=1,vetor_saida)
			END DO
		CLOSE (1)

		OPEN (1, file = path_bs)
			READ(1,*) (bs(J,1),J=1,vetor_saida)
		CLOSE (1)

		!----------------------------------------------------------------------!
		! INICIO DA REDE: FEEDFORWARD
		!----------------------------------------------------------------------!
		DO I = 1, num_pad_gen	!DO NUMERO TOTAL DE PADROES GENERALIZACAO

			! ATIVACAO: 1. CAMADA OCULTA
			vh1 = 0.d0
			vh1(:,1) = matmul(x_gen(:,I),wh1(:,:))
			vh1(:,1) = vh1(:,1) - bh1(:,1);
			SELECT CASE(f_ativa)
				CASE (1) !LOGISTICA
					yh1(:,1) = 1.d0/(1.d0+DEXP(-a * vh1(:,1)))
				CASE (2) !TANGENTE
					yh1(:,1) = (1.d0-DEXP(-vh1(:,1)))/(1.d0+DEXP(-vh1(:,1)))
				CASE (3) !GAUSS
					yh1(:,1) = DEXP(-vh1(:,1))
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
					yh2(:,1) = DEXP(-vh2(:,1))
			END SELECT
						
			! ATIVACAO: CAMADA DE SAIDA
			vs = 0.d0
			vs(:,1) = MATMUL(yh2(:,1), ws(:,:))
			vs(:,1) = vs(:,1) - bs(:,1)
			SELECT CASE(f_ativa)
				CASE (1) !LOGISTICA
					ys(:,I) = 1.d0 / (1.d0 + DEXP(-a * vs(:,1)))
				CASE (2) !TANGENTE
					ys(:,I) = (1.d0 - DEXP(-vs(:,1))) / (1.d0 + DEXP(-vs(:,1)))
				CASE (3) !GAUSS
					ys(:,I) = DEXP(-vs(:,1))
			END SELECT

			! CALCULO ERRO GENERALIZACAO
			erro(:,I) = yd_gen(:,I) - ys(:,I)

		ENDDO !NUMERO PADROES

		OPEN(12, FILE = path_result_ys)
			DO i=1,vetor_saida
				WRITE(12,110)(ys(i,j),j=1,num_pad_gen)
			ENDDO
		CLOSE(12)
110		FORMAT (424F32.16)

		!CALCULO DO ERRO QUADRATICO MEDIO
		emq(1,1) = sum(erro(:,:))
		emq(1,1) = (1.d0/(num_pad_gen))*emq(1,1)
			
		Rede_Neural_C2_GEN = emq(1,1)
				
		deallocate(x_gen)
		deallocate(yd_gen)
		deallocate(vh1)
		deallocate(wh1)
		deallocate(bh1)
		deallocate(vs)
		deallocate(ys)
		deallocate(ws)
		deallocate(bs)
		deallocate(vh2)
		deallocate(wh2)
		deallocate(bh2)
		deallocate(erro)
		deallocate(emq)
		deallocate(yh1)
		deallocate(yh2)

	!FIM DA FUNCAO
	END FUNCTION Rede_Neural_C2_GEN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	DOUBLE PRECISION FUNCTION Rede_Neural_C3_GEN(parametros, path_ent, path_sai, path_bh1, path_bh2, &
		path_bh3, path_bs, path_wh1, path_wh2, path_wh3, path_ws, path_result_ys)

		IMPLICIT NONE

		!Parametros recebidos pela funcao (passados pelo algoritmo de otimizacao)
		double precision, dimension(:,:) :: parametros
		CHARACTER*52, path_ent,path_sai
		CHARACTER*52, path_wh1, path_wh2, path_wh3, path_ws
		CHARACTER*52, path_bh1, path_bh2, path_bh3, path_bs
		CHARACTER*52, path_result_ys
		
		double precision, allocatable,dimension(:,:) :: x_gen
		double precision, allocatable,dimension(:,:) :: yd_gen
		double precision, allocatable,dimension(:,:) :: vs, vh1, vh2, vh3
		double precision, allocatable,dimension(:,:) :: ws, wh1, wh2, wh3
		double precision, allocatable,dimension(:,:) :: bs, bh1, bh2, bh3
		double precision, allocatable,dimension(:,:) :: ys, yh1, yh2, yh3
		double precision, allocatable,dimension(:,:) :: erro
		double precision, allocatable,dimension(:,:) :: emq
			
		double precision, parameter :: a = 1

		double precision	:: entrada(8,1)

		integer		:: hide_neuron(3,1) ! número de neurônios na camada escondida
		integer		:: vetor_entrada	! dimensão do vetor de entrada
		integer		:: num_pad_gen		! número de padrões a serem treinados
		integer		:: vetor_saida	 	! dimensão do vetor de saída
		integer		:: I,J				! variável para controle do programa
		integer		:: f_ativa			! controle de qual funcao ativacao
		integer		:: hide_camada		! numero de camadas escondidas

		DO I = 1, 8
			entrada(I,1) = parametros(I,1)
		END DO

		hide_camada = entrada(1,1)
		hide_neuron(1,1) = entrada(2,1)!numero de neuronios na camada escondida
		hide_neuron(2,1) = entrada(3,1)!numero de neuronios na camada escondida
		hide_neuron(3,1) = entrada(4,1)!numero de neuronios na camada escondida
		f_ativa = entrada(5,1)!parâmetro que corresponde ao tipo de funcao que será utilizada
		num_pad_gen = entrada(6,1)
		vetor_entrada = entrada(7,1)	!numero de entradas
		vetor_saida = entrada(8,1)		!numero de saidas	

		!BIAS 1. CAMADA ESCONDIDA
		allocate(bh1(hide_neuron(1,1),1))
		!BIAS 2. CAMADA ESCONDIDA
		allocate(bh2(hide_neuron(2,1),1))
		!BIAS 3. CAMADA ESCONDIDA
		allocate(bh3(hide_neuron(3,1),1))
		!BIAS CAMADA DE SAIDA
		allocate(bs(vetor_saida,1))
		!ERROS
		allocate(erro(vetor_saida,num_pad_gen))
		allocate(erro(num_pad_gen,1))
		!CAMPO LOCAL INDUZIDO
		allocate(vh1(hide_neuron(1,1),1))
		allocate(vh2(hide_neuron(2,1),1))
		allocate(vh3(hide_neuron(3,1),1))
		allocate(vs(vetor_saida,1))
		!PESOS 1. CAMADA ESCONDIDA
		allocate(wh1(vetor_entrada,hide_neuron(1,1)))
		!PESOS 2. CAMADA ESCONDIDA
		allocate(wh2(hide_neuron(1,1),hide_neuron(2,1)))
		!PESOS 3. CAMADA ESCONDIDA
		allocate(wh3(hide_neuron(2,1),hide_neuron(3,1)))
		!PESOS CAMADA DE SAIDA
		allocate(ws(hide_neuron(hide_camada,1),vetor_saida))
		!DADOS DE ENTRADA E SAIDA
		allocate(x_gen(vetor_entrada,num_pad_gen))
		allocate(yd_gen(vetor_saida,num_pad_gen))
		!SAIDAS OBTIDAS
		allocate(yh1(hide_neuron(1,1),1))
		allocate(yh2(hide_neuron(2,1),1))
		allocate(yh3(hide_neuron(3,1),1))
		allocate(ys(vetor_saida,num_pad_gen))

		!------------------------------------------------------------!
		!LENDO OS PARAMETROS DO ARQUIVO DE ENTRADA
		!------------------------------------------------------------!

		OPEN (1, file = path_sai)
			DO I = 1, vetor_saida
				READ(1,*) (yd_gen(I,J), J = 1, num_pad_gen)
			END DO
		CLOSE (1)

		OPEN (2, file = path_ent)
			DO I = 1, vetor_entrada
				READ(2,*) (x_gen(I,J), J=1, num_pad_gen)
			END DO
		CLOSE (2)

		!--------------------------------------------------------------------!
		!LENDO OS PARAMETROS: wh1, bh1, ws e bs
		!--------------------------------------------------------------------!
		OPEN (1, file = path_wh1)
			DO I = 1, vetor_entrada
				READ(1,*) (wh1(I,J),J=1,hide_neuron(1,1))
			END DO
		CLOSE (1)

		OPEN (1, file = path_wh2)
			DO I = 1, hide_neuron(1,1)
				READ(1,*) (wh2(I,J),J=1,hide_neuron(2,1))
			END DO
		CLOSE (1)

			OPEN (1, file = path_wh3)
			DO I = 1, hide_neuron(2,1)
				READ(1,*) (wh3(I,J),J=1,hide_neuron(3,1))
			END DO
		CLOSE (1)


		OPEN (1, file = path_bh1)
			READ(1,*) (bh1(J,1),J=1,hide_neuron(1,1))
		CLOSE (1)

		OPEN (1, file = path_bh2)
			READ(1,*) (bh2(J,1),J=1,hide_neuron(2,1))
		CLOSE (1)

			OPEN (1, file = path_bh3)
			READ(1,*) (bh3(J,1),J=1,hide_neuron(3,1))
		CLOSE (1)

		OPEN (1, file = path_ws)
			DO I = 1, hide_neuron(3,1)
				READ(1,*) (ws(I,J),J=1,vetor_saida)
			END DO
		CLOSE (1)

		OPEN (1, file = path_bs)
			READ(1,*) (bs(J,1),J=1,vetor_saida)
		CLOSE (1)

		!----------------------------------------------------------------------!
		! INICIO DA REDE: FEEDFORWARD
		!----------------------------------------------------------------------!
		DO I = 1, num_pad_gen	!DO NUMERO TOTAL DE PADROES GENERALIZACAO

			! ATIVACAO: 1. CAMADA OCULTA
			vh1 = 0.d0
			vh1(:,1) = matmul(x_gen(:,I),wh1(:,:))
			vh1(:,1) = vh1(:,1) - bh1(:,1);
			SELECT CASE(f_ativa)
				CASE (1) ! LOGISTICA
					yh1(:,1) = 1.d0 / (1.d0 + DEXP(-a * vh1(:,1)))
				CASE (2) ! TANGENTE
					yh1(:,1) = (1.d0 - DEXP(-vh1(:,1))) / (1.d0 + DEXP(-vh1(:,1)))
				CASE (3) ! GAUSS
					yh1(:,1) = DEXP(-vh1(:,1))
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
					yh2(:,1) = DEXP(-vh2(:,1))
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
					ys(:,I) = 1.d0 / (1.d0 + DEXP(-a * vs(:,1)))
				CASE (2) ! TANGENTE
					ys(:,I) = (1.d0 - DEXP(-vs(:,1))) / (1.d0 + DEXP(-vs(:,1)))
				CASE (3) ! GAUSS
					ys(:,I) = DEXP(-vs(:,1))
			END SELECT

		
			! CALCULO ERRO GENERALIZACAO
			erro(:,I) = yd_gen(:,I) - ys(:,I)

		ENDDO ! NUMERO PADROES

		OPEN(12, FILE = path_result_ys)
			DO i=1,vetor_saida
				WRITE(12,110)(ys(i,j),j=1,num_pad_gen)
			ENDDO
		CLOSE(12)
110		FORMAT (424F32.16)

		!CALCULO DO ERRO QUADRATICO MEDIO
		emq(1,1) = sum(erro(:,:))
		emq(1,1) = (1.d0/(num_pad_gen))*emq(1,1)

		Rede_Neural_C3_GEN = emq(1,1)

		deallocate(x_gen)
		deallocate(yd_gen)
		deallocate(vh1)
		deallocate(wh1)
		deallocate(bh1)
		deallocate(vs)
		deallocate(ys)
		deallocate(ws)
		deallocate(bs)
		deallocate(vh2)
		deallocate(wh2)
		deallocate(bh2)
		deallocate(vh3)
		deallocate(wh3)
		deallocate(bh3)
		deallocate(erro)
		deallocate(emq)
		deallocate(yh1)
		deallocate(yh2)
		deallocate(yh3)

		!Fim da funcao
	END FUNCTION Rede_Neural_C3_GEN

!Fim do modulo
END MODULE ModuloRNA_gen
