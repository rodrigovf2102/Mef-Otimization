*####################################################################
*
*	Programacao: Claudio Jose Martins DEC/CEFET-MG, Rodrigo Vieira Fonseca, PPGEC/CEFET-MG
*	e-mail: smc.engenharia@gmail.com, rodrigovf2102@gmail.com
*
*------------------------------------------------------------------------------- 
*       
*	Programa para determinação dos parâmetros modais de porticos 3D lineares
*
*		- entrada do pórtico por coordenadas
*		- e possível subdividir superelementos
*     
*------------------------------------------------------------------------------- 
*------------------------------------------------------------------------------- 
	program Portico3D
cc	use MSIMSL
	implicit none
*
	integer nse,nnp !numero de superelementos e nos principais
	integer nseva,nsevi2,nsevi3,nseve,nsevro
	integer ne,nn,ng,nprop !numero de elementos, nos e graus de liberdade
	integer nde !numero de divisoes dos superelementos
	integer ncn,nce !numero de nos com forcas e	com vinculacoes prescritas
	integer nne,ngn,nge
	integer nd !numero de dimensões (nd=3, X,Y,Z) 
	integer nvar !numero de variação nas propriedades 
	integer nvarA,nvarI2,nvarI3,nvarIt,nvarE,nvarG,nvarro
	integer nvarrigx,nvarrigy,nvarrigz,nvarrigrx,nvarrigry,nvarrigrz
	integer nvarfx,nvarfy,nvarfz,nvarmx,nvarmy,nvarmz
	integer nvarfex
	integer nvarexe,nvarkmax !numero de variacoes a serem executadas
	integer isolve,isolvee !1-linear estático 2-modal,  
	integer nmode !numero de modos a serem plotados 
	integer i,j,ie,ise,ios,ivar, iii, jjj,ll,iee,jj,ii,ivm,j1,j2,j11 
	integer k,kk,ncol,nlin,zz,z,vv,i5,i4,i6 !var. aux. 
	integer ient,iout,isaida,ent !unidades de arquivos 
	integer casovar,zxc,zxs,aux3,aux2 !1elimina variacoes iguais 2nao elimina
	integer casoexe !1 varicoes por combinacao !2variacoes aleatorias
	integer noti,casooti
	integer msim !molas simetricas? 1- Sim, 2- Nao
	integer asim,i2sim,i3sim,esim,rosim	!A,I2,I3,E,ro simetricos?
	real*8 i7, resant
	real*8 ivarreal,ireal,iereal,x,y,fnat
	real*8 tol,res !tolerancia
	real*8 alfa,gama !amortecimento de Rayleigh
	real*8 pi2,shift,c1 !var. aux
	real*8 varmolas,vara,vari2,vari3,vare,varro
	character*30 arqent,arqout, arqoutt,saida,entrad !nome arq.entrada e saida
	integer casoprint
*
	integer, allocatable :: econ(:,:)!conectividade dos elementos
	integer, allocatable :: norig(:),nofex(:) !nos com forças e molas
	integer, allocatable :: scon(:,:)!conectividade dos superelementos
	real*8, allocatable :: coord(:,:)!cordenadas espaciais dos nos
	real*8, allocatable :: pcoord(:,:)!cordenadas espaciais dos nos principais
	real*8, allocatable :: mrot(:,:)!matriz de rotacao 
	real*8, allocatable :: rig(:,:)!rigidezes das vinculacoes externas 
	real*8, allocatable :: fex(:,:)!vetor de forcas nodais 
	real*8, allocatable :: kel(:,:),mel(:,:) !matrizes dos elementos 
	real*8, allocatable :: kgl(:,:),mgl(:,:),fgl(:) !matrizes da estrut 
	real*8, allocatable :: L(:),A(:),I2(:),I3(:),It(:) !prop. geometri
	real*8, allocatable :: E(:),G(:),ro(:),angle(:) !prop. mecanicas
	real*8, allocatable :: b(:,:) !matrizes aux. para subrot. de multm, inv.
      real*8, allocatable :: bsvd(:),csvd(:),tsvd(:) !vetores aux svd
      real*8, allocatable :: s(:) !aulovalores SVD
	real*8, allocatable :: fe(:)
	real*8, allocatable :: fn(:,:), oti(:,:)	 !matriz de informaçoes modais
	double complex, allocatable :: u(:,:),v(:,:) !autovetores SVD
	double complex, allocatable :: modal(:,:) !matriz modal SVD
	real*8, allocatable :: mim(:,:)	!residuo caso1
	real*8, allocatable :: auxa(:,:)	!auxiliar
	integer, allocatable :: auxb(:)	!auxiliar
	real*8, allocatable :: minimo(:) !residuo caso2
	real*8, allocatable :: aux(:) !residuo caso2
	integer, allocatable :: var(:) !auxiliar
	real*8, allocatable :: auxv(:,:) !auxiliar
	real*8, allocatable :: sa(:),sI2(:),sI3(:),sit(:),se(:),sg(:)
	real*8, allocatable :: sro(:),sp(:),agl(:)



*------------------------------------------------------------------------------- 
	write(*,'(80A1)') ('-',i=1,80)
	write(*,*) ' Determinacao dos parametros modais de porticos 3D'
*
	arqent='Entrada.txt' !arquivo de entrada
	ient=14	 
	arqout='Saida.dat' !arquivo de saida
	iout=15
	arqoutt='Entradaotimizacao.txt' !segundo arquivo saida
	isaida=12
	ent=13
      saida='Saidaotimizacao.txt';
	entrad='Entradainformacoes.txt'



     	nd=3
	nne=2
	ngn=6
	nge=nne*ngn
	nprop=7
*
*	Abertura arquivos de entrada e saida  
	write(*,*) 
	open(ient,file=arqent,status='old',action='read',iostat=ios)
	if (ios.ne.0) stop 'ERRO ARQUIVO DE ENTRADA'
	open(iout,file=arqout,status='replace',action='write',iostat=ios)
	if (ios.ne.0) stop 'ERRO ARQUIVO DE SAIDA'
*
	write(iout,'(80A1)') ('-',i=1,80)
	write(iout,*) ' DADOS DE ENTRADA'
	write(iout,*) 
	read(ient,*) nnp
 	write(iout,*) '  Numero de nos principais:',nnp 
	read(ient,*) nse
	write(iout,*) '  Numero de superelementos:',nse 
	read(ient,*) nde
	write(iout,*) '  Numero de elementos por superelemento:',nde 
	read(ient,*) nce
 	write(iout,*) '  Numero de nos com vinculos externos:',nce 
	read(ient,*) msim
 	write(iout,*) '  Molas simetricas? 1-Sim,0-Nao:',msim 
	read(ient,*) asim,i2sim,i3sim,esim,rosim
 	write(iout,*) '  A,I2,I3,E,ro simetricos? 1-Sim, 0-Nao :',asim,
	1i2sim,i3sim,esim,rosim
	read(ient,*) ncn
 	write(iout,*) '  Numero de nos com forcas externas:',ncn
	read(ient,*) nseva,nsevi2,nsevi3,nseve,nsevro 
	read(ient,*) nvarA, nvarI2, nvarI3,nvarE,nvarro
	read(ient,*) nvarrigx,nvarrigy,nvarrigz,nvarrigrx,nvarrigry,
     1nvarrigrz
	read(ient,*) nvarfx,nvarfy,nvarfz,nvarmx,nvarmy,nvarmz
	read(ient,*) nvarexe
 	write(iout,*) '  Numero de variacoes da propriedade A: ', nvarA 
	write(iout,*) '  Numero de variacoes da propriedade I2: ', nvarI2
	write(iout,*) '  Numero de variacoes da propriedade I3: ', nvarI3
	write(iout,*) '  Numero de variacoes da propriedade It: ', nvarIt
      write(iout,*) '  Numero de variacoes da propriedade E: ', nvarE
	write(iout,*) '  Numero de variacoes da propriedade G: ' ,nvarG
	write(iout,*) '  Numero de variacoes da propriedade ro: ' ,nvarro
	write(iout,*) '  Numero de variacoes das molas: ' ,
     1(nvarrigx*nvarrigy*nvarrigz*nvarrigrx*nvarrigry*nvarrigrz)**nce
	write(iout,*) '  Numero de variacoes das forças :', 
     1(nvarfx*nvarfy*nvarfz*nvarmx*nvarmy*nvarmz)**ncn

	if(msim.eq.0) then
		varmolas=((nvarrigx*nvarrigy*nvarrigz)**nce)
		varmolas=varmolas*((nvarrigrx*nvarrigry*nvarrigrz)**nce)
	endif
	if(msim.eq.1) then
		varmolas=(nvarrigx*nvarrigy*nvarrigz)
		varmolas=varmolas*(nvarrigrx*nvarrigry*nvarrigrz)
	endif
	if(asim.eq.0) then
		vara=nvarA**nseva
	endif
	if(asim.eq.1) then
		vara=nvarA
	endif
	if(i2sim.eq.0) then
		vari2=nvari2**nsevi2
	endif
	if(i2sim.eq.1) then
		vari2=nvari2
	endif
	if(i3sim.eq.0) then
		vari3=nvari3**nsevi3
	endif
	if(i3sim.eq.1) then
		vari3=nvari3
	endif
	if(esim.eq.0) then
		vare=nvare**nseve
	endif
	if(esim.eq.1) then
		vare=nvare
	endif		
	if(rosim.eq.0) then
		varro=nvarro**nsevro
	endif
	if(rosim.eq.1) then
		varro=nvarro
	endif		

	nvar=varmolas*vara*vari2*vari3*vare*varro

	write(iout,*) '  Numero de combinaçoes :',nvar
	read(ient,*) isolve
 	write(iout,*) '  Tipo de analise:',isolve 
	read(ient,*) nmode
 	write(iout,*) '  Numero de modos:',nmode 
	read(ient,*) tol
 	write(iout,*) '  Tolerancia:',tol 
	read(ient,*) alfa,gama
	read(ient,*) casovar
	read(ient,*) casoexe
	read(ient,*) casooti
	read(ient,*) casoprint
	if(casoexe.eq.1) nvarexe=nvar
 	write(iout,*) '  Coef. amort. Rayleigh:',sngl(alfa),sngl(gama) 
	ne=nse*nde
	write(iout,*) '  Numero total de elementos:',ne 
	nn=nse*(nde-1)+nnp
	write(iout,*) '  Numero total de nos:',nn 
	ng=ngn*nn
	write(iout,*) '  Numero total de graus de liberdade:',ng 
	ll=0
	nlin=nmode*nse	   !n de linhas	de fn
	ncol=4+nprop+nce+nce*6+ncn*6+ncn
		   !n de colunas de fn

	if(casooti.eq.1) then
		goto 17
	endif

	open(ent,file=entrad,status='old',action='read',iostat=ios)
	
	read(ent,*) ne
	read(ent,*) isolvee
	read(ent,*) noti

	allocate (fe(nmode))
	allocate (mim(nmode,noti))
	allocate (auxa(nmode,noti))
	allocate (auxb(noti))
	allocate (minimo(noti))
	allocate (var(noti))
	allocate (aux(noti))
	allocate (auxv(nlin*noti,ncol))
	allocate (oti(nlin*noti,ncol))

	do i=1,nmode
		read(ent,*) fe(i)
	enddo

	close(ent)

	allocate (fn(nlin,ncol))

 
	do i=1,noti
		minimo(i)=0
		var(i)=0
	enddo

	do i=1,nmode
		do j=1,noti
			mim(i,j)=0
		enddo
	enddo

	do i=1,nlin
		do j=1,ncol
			fn(i,j) = 0
		enddo
	enddo

17	continue

		allocate (scon(nse,nne)) 
		allocate (pcoord(nnp,nd))
		allocate (econ(ne,nne))  
		allocate (coord(nn,nd))  
		allocate (angle(nse))   
		allocate (L(ne),A(ne))	 
		allocate (I2(ne),I3(ne),It(ne))  
		allocate (E(ne),G(ne),ro(ne))   
		allocate (sa(ne),sI2(ne),sI3(ne),sit(ne),se(ne),sg(ne))
		allocate (sro(ne),sp(ne),agl(ne))

	if(nce.ne.0)then
		allocate (fex(nnp,ngn))
		allocate (nofex(ncn))
	endif
	if(nce.ne.0)then
		allocate (rig(nnp,ngn))
		allocate (norig(nce))
	endif

*
		call Entrada(nse,nnp,nn,ne,nd,nne,nde,scon,pcoord,econ
	1,coord,ient,iout,L,sa,sI2,sI3,sit,se,sg,sro,sp,c1,agl,
     2ngn,nprop,angle)

	nvarkmax=nvarrigx
	if(nvarkmax.lt.nvarrigy) nvarkmax=nvarrigy
	if(nvarkmax.lt.nvarrigz) nvarkmax=nvarrigz	
	if(nvarkmax.lt.nvarrigrx) nvarkmax=nvarrigrx
	if(nvarkmax.lt.nvarrigry) nvarkmax=nvarrigry
	if(nvarkmax.lt.nvarrigry) nvarkmax=nvarrigrz

	zxc=0
	ivarreal=0
	kk=0

*
*	Looping sobre a variacao das propriedades fisicas
	do aux2=1,nvarexe
		ivar=ivar+1
		if(ivar.ge.(nvarexe+1)) goto 1777
		ivarreal=ivar

 

		call Parametros(casovar,nse,nnp,ne,msim,
     1A,ro,I2,I3,It,E,G,angle,nce,ncn,ngn,rig,fex,ient,iout
     2,nvarA,nvarI2,nvarI3,nvarIt,nvarE,nvarG,nvarro,aux3,nvarkmax,
     3nvarrigx,nvarrigy,nvarrigz,nvarrigrx,nvarrigry,nvarrigrz,norig,
     4nofex,nvarfx,nvarfy,nvarfz,nvarmx,nvarmy,nvarmz,ivar,nvarexe,
     5sa,sI2,sI3,sit,se,sg,sro,sp,c1,agl,nprop,nlin,ncol,casoexe,
     6nseva,nsevi2,nsevi3,nseve,nsevro,asim,i2sim,i3sim,esim,rosim)


		allocate (kgl(ng,ng),mgl(ng,ng),fgl(ng))
	 
 
*		Calcula matrizes de rigidez e massa da estrutura

		!write(*,*) '  ...Calcula matrizes rigidez e massa da estrutura'

		allocate (kel(nge,nge),mel(nge,nge),mrot(nge,nge))
	 
		do ie=1,nse
*
		!	write(*,*) ' ....matriz de rotacao para superelemento',ise

			call mrotacao(mrot,ie,nge,nse,nnp,nne,nd
	1		,scon,pcoord,angle(ie),tol)
*
c				write(*,*) '	  ......matriz elemento',ie
			 				
				call Lmatriz(nd,nge,A(ie),L(ie),E(ie)
	1			,G(ie),I2(ie),I3(ie),It(ie)
     2			,ro(ie),kel,mel,mrot,tol)


		   call Gmatriz(ngn,nge,nne,ng,ne,ie,econ,kel,mel,kgl,mgl,tol)
				
		enddo
		deallocate(kel,mel,mrot)
*
*		Calcula matrizes de rigidez e massa da estrutura
	
		!write(*,*) '  ...Introducao condicoes de contorno'

		call Contorno(ng,nnp,ngn,rig,fex,kgl,fgl,tol,nvarexe,ivar,mgl)


*
		open(iout,file=arqout,status='old',access='append',iostat=ios)
		if (ios.ne.0) stop 'ERRO ARQUIVO DE SAIDA'


	!do i=1,ngn*nnp
	!	if(i.eq.1) write(iout,1588) (j,j=1,ngn*nnp)
	!	write(iout,1587) i,(kgl(i,j),j=1,ngn*nnp)
	!	write(*,*)
	!enddo
1587  format(1x,i2,1x,e8.2e2,1x,e8.2e2,1x,e8.2e2,1x,e8.2e2,1x,e8.2e2,1x,
	1e8.2e2,1x,e8.2e2,1x,e8.2e2,1x,e8.2e2,1x,e8.2e2,1x,
     2e8.2e2,1x,e8.2e2,1x,e8.2e2,1x,e8.2e2,1x,e8.2e2,1x,
     3e8.2e2,1x,e8.2e2,1x,e8.2e2,1x,e8.2e2,1x,e8.2e2,1x,
     4e8.2e2,1x,e8.2e2,1x,e8.2e2,1x,e8.2e2,1x,e8.2e2,1x,
     5e8.2e2,1x,e8.2e2,1x,e8.2e2,1x,e8.2e2,1x,e8.2e2,1x)
1588	format(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,
	1i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,
     2i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,
     3i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,
     4i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,
     5i8,1x,i8,1x,i8,1x,i8,1x,i8,1x)
     	


		if(casoprint.eq.1) then

			write(iout,*)
			write(iout,*)
			write(iout,*) ' DADOS DE SAIDA, variacao',ivar
			write(iout,*)

		endif
 
*			Analise modal
				allocate(modal(ng,ng),b(ng,ng),bsvd(ng),csvd(ng)
     1				,tsvd(ng))

				shift=0.0d0


				call mshift(ng,kgl,mgl,modal,shift)
				call inverse(modal,b,ng)


				call cmultmm(ng,b,mgl,modal)



				!write(*,*) '  ...Resolvendo autovalor e autovetor'

				allocate(u(ng,ng),v(ng,ng),s(ng))

 
				call dcsvd (modal,ng,ng,0,ng,ng,s,u,v
	1			,bsvd,csvd,tsvd)

				if (nmode.gt.ng .or. nmode.le.0) nmode=ng 
*
				pi2=4.0d0*dasin(1.0d0)
				kk=0
				do i=1,nmode
				kk=kk+1
					ireal=i
					s(i)=dsqrt(1.0d0/(s(i)+shift))/pi2 !converte para Hz

					if(casoprint.eq.1) then
						do iii=1, ng/6
							jjj=iii*6+1
							write(iout,01) i,s(i)
	1						,(dreal(v(j,i)),j=1*(jjj-6),6*iii)
   						enddo
					endif

				if(casooti.eq.1) then	
					write(*,*) i,s(i),ivar
				endif
	!gerando matriz de informacoes modais
				if(casooti.eq.2) then				
					do j=nse,1,-1
				    	k=0
						iereal=j
						if(j.ne.nse) kk=kk+1
						k=k+1
						fn(kk,k)=ivarreal
						k=k+1
						fn(kk,k)=ireal
						k=k+1
						fn(kk,k)=s(i)
						k=k+1
						fn(kk,k)=iereal
						k=k+1
						fn(kk,k)=E(j)
						k=k+1
						fn(kk,k)=A(j)
						k=k+1
						fn(kk,k)=I2(j)
						k=k+1
						fn(kk,k)=I3(j)
						k=k+1 
						fn(kk,k)=It(j)
						k=k+1
						fn(kk,k)=G(j)
						k=k+1
						fn(kk,k)=ro(j)
						do jj=1,nce
							k=k+1
							fn(kk,k)=norig(jj)
							k=k+1	  
							fn(kk,k)=rig(norig(jj),1)
							k=k+1
							fn(kk,k)=rig(norig(jj),2)	
							k=k+1
							fn(kk,k)=rig(norig(jj),3)
							k=k+1  
							fn(kk,k)=rig(norig(jj),4)
							k=k+1
							fn(kk,k)=rig(norig(jj),5)	
							k=k+1
							fn(kk,k)=rig(norig(jj),6)
						enddo
						if(ncn.ne.0)then
						do ii=1,ncn
							k=k+1
							fn(kk,k)=nofex(ii)
							k=k+1
							fn(kk,k)=fex(nofex(ii),1)
							k=k+1
							fn(kk,k)=fex(nofex(ii),2)	
							k=k+1
							fn(kk,k)=fex(nofex(ii),3)
							k=k+1  
							fn(kk,k)=fex(nofex(ii),4)
							k=k+1
							fn(kk,k)=fex(nofex(ii),5)	
							k=k+1
							fn(kk,k)=fex(nofex(ii),6)
						enddo
						endif
					enddo
 				endif
			enddo
			deallocate(modal,b,bsvd,csvd,tsvd,v,s,u) 
			deallocate(kgl,mgl,fgl)
	
     	if(ivar.eq.0) write(*,*)				
	if(ivar.eq.1)then					
 		if(casooti.eq.2) then
	open(isaida,file=saida,status='replace',action='write',iostat=ios)
		endif
	endif

	if(ivar.eq.1) then
		do i=1,20
			write(iout,*)'------------------------------------------'
		enddo
	endif
	j=0
	k=0
	z=0
	vv=0
	!rotina para selecionar residuos otimos
	if(ivar.eq.1)then
		do i=1,noti
			minimo(i)=999999
		enddo
	endif
	if(ivar.eq.1)then
		do j1=1,noti*nlin
			do j2=1,ncol
				auxv(j1,j2)=0
				oti(j1,j2)=0
			enddo
		enddo
	endif
		res=0
		ii=0
		x=0
		y=0
		do i=1,(nmode*ne),ne
			ii=ii+1 !ii é a linha das frequencias experi.	 						 
			x=(fe(ii)-fn(i,3))/fe(ii)
 			if(x.lt.0) x=-x
			res = x + res
			if(ii.eq.nmode) res=res*100/nmode
			if(ivar.eq.1) fnat=100
		enddo

	!impressao no meio da caliracao, para acompanhamento
	zxc=zxc+1
	if((zxc.eq.(2000/nse**2).or.res.le.0.01*fnat).and.ivar.gt.1) then 
		write(*,*)'Variacao',ivar,'   de',nvarexe, '    Residuo',res
		zxc=0
	endif

	!rotinas de otimização (pulo de analises modais de alto residuo)

		i7=res-resant
		resant=res
		if(i7.lt.0) i7=-i7

 
		if(i7.le.0.001*fnat.and.res.ge.0.15*fnat)then
				aux3=1
				do zz=1,(nvarkmax/5)
				ivar=ivar+1
				i5=i5+1
			call Parametros(casovar,nse,nnp,ne,msim,
     1A,ro,I2,I3,It,E,G,angle,nce,ncn,ngn,rig,fex,ient,iout
     2,nvarA,nvarI2,nvarI3,nvarIt,nvarE,nvarG,nvarro,aux3,nvarkmax,
     3nvarrigx,nvarrigy,nvarrigz,nvarrigrx,nvarrigry,nvarrigrz,norig,
     4nofex,nvarfx,nvarfy,nvarfz,nvarmx,nvarmy,nvarmz,ivar,nvarexe,
     5sa,sI2,sI3,sit,se,sg,sro,sp,c1,agl,nprop,nlin,ncol,casoexe,
     6nseva,nsevi2,nsevi3,nseve,nsevro,asim,i2sim,i3sim,esim,rosim)
				enddo
		endif
	
	!rotinas de otimização (pulo de analises modais de alto residuo)
		if(res.gt.(0.20*fnat)) then
			aux3=1
			do zz=1,(nvarkmax/3)
			ivar=ivar+1
			i4=i4+1
			call Parametros(casovar,nse,nnp,ne,msim,
     1A,ro,I2,I3,It,E,G,angle,nce,ncn,ngn,rig,fex,ient,iout
     2,nvarA,nvarI2,nvarI3,nvarIt,nvarE,nvarG,nvarro,aux3,nvarkmax,
     3nvarrigx,nvarrigy,nvarrigz,nvarrigrx,nvarrigry,nvarrigrz,norig,
     4nofex,nvarfx,nvarfy,nvarfz,nvarmx,nvarmy,nvarmz,ivar,nvarexe,
     5sa,sI2,sI3,sit,se,sg,sro,sp,c1,agl,nprop,nlin,ncol,casoexe,
     6nseva,nsevi2,nsevi3,nseve,nsevro,asim,i2sim,i3sim,esim,rosim)
			enddo
		endif
		aux3=0

	!formacao da matriz de parametros, associadas a frequencias otimas
		do iii=1,noti

		if(res.lt.minimo(iii).and.ivar.ne.1) then

			do jj=1,noti
				aux(jj)=minimo(jj)
				auxb(jj)=var(jj)			
			enddo
			do j1=1,nlin*noti,1		
				do j2=1,ncol,1
					auxv(j1,j2)=oti(j1,j2)
				enddo
			enddo

		    if(iii.ne.noti) then
				do j1=1+(iii-1)*nlin,(noti-1)*nlin,1
					do j2=1,ncol,1
						oti(j1+nlin,j2)=auxv(j1,j2)
					enddo												 
				enddo
				do jjj=iii,noti
					if(jjj.ne.noti) then
						minimo(jjj+1)=aux(jjj)
						var(jjj+1)=auxb(jjj)	
					endif
				enddo
			endif

			j11=0
   			do j1=1+(iii-1)*nlin,(iii)*nlin,1
				j11=j11+1
				do j2=1,ncol,1
					if(j11.eq.nlin+1) j11=1
					oti(j1,j2)=fn(j11,j2)
				enddo
			enddo
			minimo(iii)=res				 !recebe residuo minimo
			var(iii)=ivar				 !recebe o nvar
			goto 133
		endif
	    enddo
133    continue
	
	 zxs=zxs+1
	 if(zxs.eq.25000/(nse**2)) then
		zxs=0
		goto 137
	 endif

139	continue

	enddo

137	continue
1777	continue

	do iii=1,noti
	!impressao de resultado no arquivo executavel (case(2))
	write(*,*)'Rotinas de redução do Custo Comp. foram executadas'
	1,i5,i4
	write(*,*)'Segue abaixo as propriedades do',iii,'menor residuo'
	write(*,74) minimo(iii),var(iii)
	write(*,*) 'Segue abaixo as propriedades da variacao'			  
	write(*,774)'Frequencias Naturais'
	1,(oti((iii-1)*nlin+j,3),j=1,nse*nmode,nse)
	write(*,73) 'ELE','E','A','I2','I3','It',
	1'G','ro'
      z=var(iii)
	do i=0,ne-1
		write(*,72) (oti(((iii)*nlin-i),j), j=4,4+nprop )
  	enddo

	if(nce.ne.0) then
	write(*,*)'Segue as configuracoes de pontos nodais com molas'
	write(*,101) 'NO','Kx','Ky','Kz','Krx','Kry','Krz'
	write(*,102)(oti(((iii)*nlin-i),j)
	1,j=4+nprop+1,5+nprop+(nce*6+1))
	endif

	if(ncn.ne.0) then
	write(*,*)'Segue as configuracoes de pontos nodais com forcas'
	write(*,101) 'NO','Fx','Fy','Fz','Mx.','My.','Mz.'
	write(*,102) (oti((iii*nlin-i),j), j=5+nprop+nce*6+2,ncol)
	endif
		
	write(*,*)
	write(*,*)

		if(ivar.lt.nvarexe) goto 139

	enddo


	!impressao de resultados no arquivo de saida  (case(2))
	do iii=1,noti
	write(isaida,*)'Segue abaixo as prop. do',iii,'menor residuo'
	write(isaida,74) minimo(iii),var(iii)
	write(isaida,*) 'Segue abaixo as propriedades da variacao'
	write(isaida,774)'Freq:',(oti((iii-1)*nlin+j,3),j=1,nse*nmode,nse)
	write(isaida,73) 'ELE','E','A','I2','I3','It',
	1'G','ro'
      z=var(iii)
	do i=0,ne-1
		write(isaida,72) (oti((noti*ne*nmode-i),j), j=4,4+nprop )
  	enddo

	if(nce.ne.0) then
	write(isaida,*)'Segue as configuracoes de pontos nodais com molas'
	write(isaida,101) 'NO','Kx','Ky','Kz','Krx','Kry','Krz'
	write(isaida,102)(oti(((iii)*nlin-i),j),
	1j=4+nprop+1,5+nprop+(nce*6+1))
	endif

	if(ncn.ne.0) then
	write(isaida,*)'Segue as config. de pontos nodais com forcas'
	write(isaida,101) 'NO','Fx','Fy','Fz','Mx.','My.','Mz.'
	write(isaida,102)(oti((noti*ne*nmode-i),j),j=5+nprop+nce*6+2,ncol)
	endif
	write(isaida,*)
	write(isaida,*)
	write(isaida,*)

	
	enddo


	read(*,*)

	close(isaida)
	deallocate(fe,fn,mim)



  	close(iout)
	close(ient)



	!formatos de escrita e leitura
71    format(3x,a1,1x,i2.0,1x,a21,1x,i2.0,1x,a11,1x,f7.3,1x,a11,1x,i4.0)

72	format(1x,f4.0,4x,E8.3e2,4x,
	1 E8.3e2,4x,E8.3e2,4x,
	2 E8.3e2,4x,E8.3e2,4x,E8.3e2,4x,E8.3e2)

73	format(2x,a3,8x,a1,11x,a1,11x,a2,9x,a2,10x,a2
     1 ,10x,a1,11x,a2)

74    format(1x,'O residuo tem o valor de',1x,f7.3,1x,
	1'na variacao',1x,i9.2)

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

774   format(1x,a20,3x,f8.4,3x,f8.4,3x,f8.4,3x,f8.4,3x,f8.4,3x,f8.4,3x,
	1f8.4,3x,f8.4,3x,f8.4,3x,f8.4,3x,f8.4,f8.4,3x,f8.4,3x,f8.4,3x,f8.4)




01	format(i3,E13.6E2,<ng>E9.2E2)
	deallocate(coord,econ,pcoord,scon) 
	deallocate(L,A,ro,I2,I3,It,E,G,angle)   
	deallocate(rig,fex) 
	write(*,*) ' Fim das analises'
	write(*,'(80A1)') ('-',i=1,80)
	read(*,*) 
	end
*-------------------------------------------------------------------------------        
*	subrotina que executa entrada de dados 
*
*	obs.: nos e conectividades no sentido horario
	subroutine Entrada(nse,nnp,nn,ne,nd,nne,nde,scon,pcoord,
     1econ,coord,ient,iout,L,sa,sI2,sI3,sit,se,sg,sro,sp,c1,agl
     2,ngn,nprop,angle)

	implicit none


	integer nse,nnp,ne,nn 
	integer nne,nd,nde,ngn,nprop   
	integer i,j,ie,ise,in,ni,nf,k,je	
	integer ient,iout
	real*8 c1
	real*8, allocatable :: dl(:)

	real*8 L(ne)
      integer scon(nse,nne),econ(ne,nne)
	real*8 pcoord(nnp,nd),coord(nn,nd)
	real*8 sa(nse),sI2(nse),sI3(nse),sit(nse),se(nse)
	real*8 sro(nse),sp(nse),agl(nse),angle(nse),sg(nse)

 

*	leitura parâmetros globais
	write(iout,*) 
	write(iout,*) '  Coordenadas (x,y,z) dos nos principais'
	do i=1,nnp
		read(ient,*) (pcoord(i,j),j=1,nd)
		write(iout,*) '     No principal',i,(sngl(pcoord(i,j)),j=1,nd)
	enddo
	
	write(iout,*) 
	write(iout,*) '		Conectividade (ni,nf) dos superelementos'
	do i=1,nse
		read(ient,*) (scon(i,j),j=1,nne)
		write(iout,*) '    Superelemento',i,(scon(i,j),j=1,nne)
	enddo


	write(iout,*) 
	write(iout,*) '  Prop fisicas/geometrias superelementos'
	write(iout,*) '  It,I2,I3 inercias em torno eixos 1(torcao),2,3'

      do ise=1,nse
		write(iout,01)
		read(ient,*) je,sa(ise),sI2(ise),sI3(ise),sIt(ise),
     1	se(ise),sp(ise),sro(ise),agl(ise)
		write(iout,'(i9,2x,8E13.6E2)') je,sa(ise),sI2(ise),
	1	sI3(ise),sit(ise),se(ise),sp(ise),sro(ise),agl(ise)
		angle(ise)=agl(ise)			
		sg=se/(2.0d0*(1.0d0+sp))
		if((je.lt.1).or.(je.gt.nse)) stop 
	enddo
	
		allocate (dl(nd))
		in=0 !contador de nos
		do i=1,nnp
		in=in+1
		do j=1,nd
			coord(in,j)=pcoord(i,j)
		enddo
	enddo


	ie=0 !contador de elementos
	do ise=1,nse

			ni=scon(ise,1)
			nf=scon(ise,2)
			c1=1/dfloat(nde)

			do k=1,nd
				dl(k)=c1*(pcoord(nf,k)-pcoord(ni,k))
			enddo

		ie=ie+1
		econ(ie,1)=ni
		econ(ie,2)=nf


		L(ie)=0.0d0
		do k=1,nd
			L(ie)=L(ie)+dl(k)**2
		enddo
		L(ie)=dsqrt(L(ie))
		
	enddo

	deallocate(dl)

01	format(5x,'Superel.',7x,'Area',4x,'Inerc Min',4x,'Inerc Max',4x, 
	1       'Inerc Tor',3x,'Mod. Elast',6x,'Poisson',3x,'Mass Espec',4x,
	2       'Rot Secao')
	return
	end

*-------------------------------------------------------------------------------        
*	subrotina que executa saida de dados analise estatica
*
	subroutine SaidaCase1(nn,ngn,ng,fpl,iout)
	implicit none
*
	integer nn,ng,ngn   
	integer i,j,gi 
	integer iout  
*
	real*8 fpl(ng) 
	write(iout,*) ' Deslocamentos nodais'
	write(iout,06) 
	do i=1,nn
		gi=ngn*(i-1)
		write(iout,01) i,(fpl(gi+j),j=1,ngn)
	enddo
*
01	format(I5,<ngn>E13.6E2)
06	format(3x'No',6x,'Trans X',6x,'Trans Y',6x,'Trans Z',8x
	1       'Rot X',8x,'Rot Y',8x,'Rot Z')
	return
	end
