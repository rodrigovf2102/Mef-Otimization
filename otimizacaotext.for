	program otimizacao
	implicit none


	!declaracao das variaveis do programa
	real*8,allocatable :: fe(:) !frequencias experimentais
	real*8,allocatable :: fn(:,:) !frequencias numericas e parametros
	real*8,allocatable :: res(:)  !vetor de residuos
	integer mv !modos de vibração
	integer i,j,ii,k,kk,z,zz,iii,jjj,vv,jj,zzz,xx !contadores e auxiliares
	real x !aux residuos
	character*30 entrada,saida,entrad	!nome arquivos de entrada e saida
	integer ient,isaida, ios,ent,isolve	 !variaveis auxiliares
	integer nprop !numero de propriedades variadas
	integer nvar !numero total de variações
	integer noti
	integer nlin,ncol !numero de linhas e colunas da matriz fn
	integer ne !numero de elementos do modelo numerico
	integer nce !numero de nós com molas
	integer ncc !numero de nos com forças
	real*8, allocatable :: mim(:,:)	!residuo caso1
	real*8, allocatable :: auxa(:,:)	!auxiliar
	integer, allocatable :: auxb(:)	!auxiliar
	real*8, allocatable :: minimo(:) !residuo caso2
	real*8, allocatable :: aux(:) !residuo caso2
	integer, allocatable :: var(:) !auxiliar
	integer, allocatable :: auxv(:) !auxiliar
		

	!abertura de arquivos de entrada e saida
	ient=11
	isaida=12
	ent=13
      saida='Saidaotimizacao.txt'
	entrada='Entradaotimizacao.txt'
	entrad='Entradainformacoes.txt'
	open(ient,file=entrada,status='old',action='read',iostat=ios)
	open(isaida,file=saida,status='replace',action='write',iostat=ios)
	open(ent,file=entrad,status='old',action='read',iostat=ios)


	!definicao de diversas variaveis
	nprop =7		 !numero de propriedades variadas
	write(*,*)'Obtendo numero de modos a serem otimizados...'
	read(ent,*) mv
	write(*,*)'Obtendo numero de elementos do modelo numerico...'
	read(ent,*) ne
	write(*,*)'Obtendo variacoes...'
	read(ent,*) nvar
	read(ent,*) isolve
	read(ent,*) nce
	read(ent,*) ncc
	read(ent,*) noti
	nlin=mv*ne*nvar	   !n de linhas
	ncol=4+nprop+nce+nce*6+ncc*6+1	   !n de colunas


	!alocacao de vetores e matrizes
	allocate (fe(mv))
	allocate (fn(nlin,ncol))
	allocate (res(nvar))
	allocate (mim(mv,noti))
	allocate (auxa(mv,noti))
	allocate (auxb(noti))
	allocate (minimo(noti))
	allocate (var(noti))
	allocate (aux(noti))
	allocate (auxv(noti))

	!contador pra zerar vetores e matrizes
      do i=1,mv
		fe(i)=0
	enddo
	
	do i=1,noti
		minimo(i)=0
		var(i)=0
	enddo

	do i=1,mv
		do j=1,noti
			mim(i,j)=0
		enddo
	enddo

	do i=1,nlin
		do j=1,ncol
			fn(i,j) = 0
		enddo
	enddo

	do i=1,nvar
			res(i)=0
	enddo


	!obtenção de frequencias experimentais
	write(*,*)'Recebendo frequencias naturais experimentais...'
	do i=1,mv
		write(*,*) 'Frequencia...',i
		read(ent,*) fe(i)
	enddo
	close(ent)


	!leitura da primeira linha de caracteres do ient
	read(ient,*)


	!obtenção da matriz de informações do modelo numerico
	do i=1,nlin
			read(ient,75) (fn(i,j),j=1,ncol)
	enddo
	close(ient)


	!selecao das frequencias numericas mais proximas das experimentais
	select case(isolve)

	case(1) !Selecao de residuos minimos individuais de cada modo
	j=1
	do i=1,mv

		if(i.ne.1) j=j+ne

		k=0

			do ii=j,nlin,(mv*ne)

				k=k+1
				res(k)=fe(i)-fn(ii,3)

					if(res(k).lt.0) then

						res(k)=-res(k)

					endif

					do iii=1,noti

						if(k.eq.1) then

	    					mim(i,iii)=999999
							var(iii)=k

						endif

					enddo
					
					do iii=1,noti

						if(res(k).lt.mim(i,iii)) then

							do jj=1,noti

								auxa(i,jj)=mim(i,jj)
								auxb(jj)=var(jj)
											
							enddo

							if(iii.ne.noti) then

								do jjj=iii,noti

									if(jjj.ne.noti) then

										mim(i,jjj+1)=auxa(i,jjj)
										var(jjj+1)=auxb(jjj)

									endif

								enddo

							endif

							mim(i,iii)=res(k)
							var(iii)=k
							goto 11

						endif

					enddo
11			continue
			enddo

	!Escrita das saidas no executavel
			write(*,*)
			write(*,*)
 			write(*,*)

	        do iii=1,noti

			write(*,*)
			write(*,*)		
			write(*,71)'O',iii,'menor residuo do modo',i,
	1'tem o valor',mim(i,iii),'na variacao',var(iii)

	write(*,*)'  As propriedades dos elementos e nos neste caso sao:'
			write(*,73) 'ELE','E','A','I2','I3','It',
	1'G','ro'

			do zz=var(iii),var(iii)+ne-1
				write(*,72) (fn(zz,kk), kk=4,4+nprop )
			enddo

			if(nce.ne.0) then
				 write(*,101) 'NO','Kx','Ky','Kz',
     2			 'Krx','Kry','Krz'
				write(*,102)(fn(zz,kk),kk=4+nprop+1,5+nprop+(nce*6+1))
			endif

			if(ncc.ne.0) then
				write(*,101) 'NO','Fx','Fy','Fz',
     2			'Mx.','My.','Mz.'
			   write(*,102)(fn((zz*ne*mv-i),j),j=5+nprop+nce*6+2,ncol)
			endif
			enddo
 	!Escrita das saidas no txt
			write(isaida,*)
			write(isaida,*)
 			write(isaida,*)

	        do iii=1,noti

			write(isaida,*)
			write(isaida,*)		
			write(isaida,71)'O',iii,'menor residuo do modo',i,
	1'tem o valor',mim(i,iii),'na variacao',var(iii)

	write(isaida,*)'  As prop. dos elementos e nos neste caso sao:'
			write(isaida,73) 'ELE','E','A','I2','I3','It',
	1'G','ro'

			do zz=var(iii),var(iii)+ne-1
				write(isaida,72) (fn(zz,kk), kk=4,4+nprop )
			enddo

			if(nce.ne.0) then
				 write(isaida,101) 'NO','Kx','Ky','Kz',
     2			 'Krx','Kry','Krz'
		   write(isaida,102)(fn(zz,kk),kk=4+nprop+1,5+nprop+(nce*6+1))
			endif

			if(ncc.ne.0) then
				write(isaida,101) 'NO','Fx','Fy','Fz',
     2			'Mx.','My.','Mz.'
		  write(isaida,102)(fn((zz*ne*mv-i),j),j=5+nprop+nce*6+2,ncol)
			endif
			enddo 
	enddo
	

	case(2) !Selecao de residuos minimosde todos os modos
	j=0
	k=0
	z=0
	vv=0

	do i=1,noti
		minimo(i)=999999
	enddo

	do k=1,nvar

		ii=0
		j=1+(k-1)*ne*mv

		do i=j,(mv*ne*k),ne

			ii=ii+1   !ii é a linha das frequencias experi.	 						 
			x=fe(ii)-fn(i,3)	  !x é o residuo individual
			if(x.lt.0)	x=-x	 !pegando o modulo de x
			res(k)=res(k)+x		 !res(k) é residuo acumulado dos modos

		enddo

		do iii=1,noti

		if(res(k).lt.minimo(iii)) then

			do jj=1,noti

				aux(jj)=minimo(jj)
				auxv(jj)=var(jj)
							
			enddo

		   if(iii.ne.noti) then

				do jjj=iii,noti

					if(jjj.ne.noti) then
					
					minimo(jjj+1)=aux(jjj)
					var(jjj+1)=auxv(jjj)

					endif
				enddo

			endif
							
			minimo(iii)=res(k)				 !recebe residuo minimo
			var(iii)=k							     !recebe o nvar

			goto 13

		endif

	    enddo

13    continue
	enddo

	do iii=1,noti
	!impressao de resultado no arquivo executavel (case(2))
	write(*,*)'Segue abaixo as propriedades do',iii,'menor residuo'
	write(*,74) minimo(iii),var(iii)
	write(*,*) 'Segue abaixo as propriedades da variacao'
	write(*,73) 'ELE','E','A','I2','I3','It',
	1'G','ro'
	z=var(iii)
	do i=0,ne-1
		write(*,72) (fn((z*ne*mv-i),j), j=4,4+nprop )
  	enddo

	if(nce.ne.0) then
	write(*,*)'Segue as configuracoes de pontos nodais com molas'
	write(*,101) 'NO','Kx','Ky','Kz','Krx','Kry','Krz'
	write(*,102) (fn((z*ne*mv-i),j), j=4+nprop+1,5+nprop+(nce*6+1))
	endif

	if(ncc.ne.0) then
	write(*,*)'Segue as configuracoes de pontos nodais com forcas'
	write(*,101) 'NO','Fx','Fy','Fz','Mx.','My.','Mz.'
	write(*,102) (fn((z*ne*mv-i),j), j=5+nprop+nce*6+2,ncol)
	endif
	write(*,*)
	write(*,*)
	write(*,*)

	!impressao de resultados no arquivo de saida  (case(2))
	write(isaida,*)'Segue abaixo as prop. do',iii,'menor residuo'
	write(isaida,74) minimo(iii),var(iii)
	write(isaida,*) 'Segue abaixo as propriedades da variacao'
	write(isaida,73) 'ELE','E','A','I2','I3','It',
	1'G','ro'
	z=var(iii)
	do i=0,ne-1
		write(isaida,72) (fn((z*ne*mv-i),j), j=4,4+nprop )
  	enddo

	if(nce.ne.0) then
	write(isaida,*)'Segue as configuracoes de pontos nodais com molas'
	write(isaida,101) 'NO','Kx','Ky','Kz','Krx','Kry','Krz'
	write(isaida,102)(fn((z*ne*mv-i),j),j=4+nprop+1,5+nprop+(nce*6+1))
	endif

	if(ncc.ne.0) then
	write(isaida,*)'Segue as config. de pontos nodais com forcas'
	write(isaida,101) 'NO','Fx','Fy','Fz','Mx.','My.','Mz.'
	write(isaida,102) (fn((z*ne*mv-i),j), j=5+nprop+nce*6+2,ncol)
	endif
	write(isaida,*)
	write(isaida,*)
	write(isaida,*)

	enddo


	!termino do selec case
	end select
	close(isaida)
	deallocate(fe,fn,res,mim)


	!formatos de escrita e leitura
71    format(3x,a1,1x,i2.0,1x,a21,1x,i2.0,1x,a11,1x,f7.3,1x,a11,1x,i4.0)

72	format(2x,f3.0,4x,E8.3e2,4x,
	1 E8.3e2,4x,E8.3e2,4x,
	2 E8.3e2,4x,E8.3e2,4x,E8.3e2,4x,E8.3e2)

73	format(2x,a3,8x,a1,11x,a1,11x,a2,9x,a2,10x,a2
     1 ,10x,a1,11x,a2)

74    format(1x,'O residuo tem o valor de',1x,f7.3,1x,
	1'na variacao',1x,i6.2)

75	format(1x,f6.0,3x,f2.0,3x,f8.3,4x,f2.0,4x,E8.3e2,4x,
	1E8.3e2,4x,E8.3e2,4x,
	2E8.3e2,4x,E8.3e2,4x,E8.3e2,4x,E8.3e2,2x,
     3f1.0,2x,E10.3e2,
     42x,E10.3e2,2x,E10.3e2,2x,E10.3e2,2x,E10.3e2,2x,E10.3e2,2x,
     5f1.0,2x,E10.3e2,
     62x,E10.3e2,2x,E10.3e2,2x,E10.3e2,2x,E10.3e2,2x,E10.3e2,2x,
     7f1.0,2x,E10.3e2,
     82x,E10.3e2,2x,E10.3e2,2x,E10.3e2,2x,E10.3e2,2x,E10.3e2)

76	format(1x,f4.0,3x,f3.0,4x,f7.3,4x,f2.0,4x,E8.3e2,4x,
	1 E8.3e2,4x,E8.3e2,4x,
	2 E8.3e2,4x,E8.3e2,4x,E8.3e2,4x,E8.3e2)

77	format(1x,a3,3x,a4,4x,a4,5x,a3,8x,a1,11x,a1,11x,a2,9x,a2,10x,a2
     1 ,10x,a1,11x,a2)

99    format(3x,f2.0,2x,6E7.2e2,/,3x,f2.0,2x,6E7.2e2)

100   format(3x,a2,4x,a2,5x,a2,5x,a2,5x,a3,4x,a3,4x,a3)

101	format(2x,a3,8x,a2,11x,a2,11x,a2,9x,a3,9x,a3
     1 ,9x,a3)

102   format(3x,f2.0,4x,E8.3e2,4x,
	1 E8.3e2,4x,E8.3e2,4x,
	2 E8.3e2,4x,E8.3e2,4x,E8.3e2)

199   format(1x,31a,1x,i2.0,1x,12a)

      end program





