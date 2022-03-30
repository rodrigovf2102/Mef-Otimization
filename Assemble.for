*-------------------------------------------------------------------------------        
*	subrotina que determina a matriz de rotacao	do superelemento
*	fonte:			
*     obs: efetua duas rotacoes (r3 e r2) para concidir o eixo 1 (X) com eixo da barra   
	subroutine mrotacao(mrot,ise,nge,nse,nnp,nne,nd
	1,scon,pcoord,agl,tol)
	implicit none
*
	integer nne,nge,nse,nnp,nd 
	integer ise 
	integer i,j,k,ni,nf,ib,iflag 
	real*8 L,aux,teta,beta
	real*8 sr1,cr1,sr2,sr3,cr2,cr3,agl,tol  
*
	integer scon(nse,nne) 
	real*8 pcoord(nnp,nd) 
	real*8 mrot(nge,nge) 
	real*8, allocatable :: dl(:),cdir(:),cf(:),ci(:)  
	real*8, allocatable :: rot(:,:),rot1(:,:),rot2(:,:),rot3(:,:) 							  ''
*
	allocate(rot(nd,nd),rot1(nd,nd),rot2(nd,nd),rot3(nd,nd))
	allocate(dl(nd),cdir(nd),cf(nd),ci(nd))
*
      ni=scon(ise,1)
	nf=scon(ise,2)
	L=0.0d0
	do k=1,nd
		ci(k)=pcoord(ni,k)
		cf(k)=pcoord(nf,k)
		dl(k)=(cf(k)-ci(k))
		L=L+dl(k)**2
	enddo
	L=dsqrt(L)

	do k=1,nd
		cdir(k)=dl(k)/L	!cossenos diretores
	enddo
	!write(*,*) ise
*
	if((dabs(cf(1)).ge.tol).or.(dabs(ci(1)).ge.tol)) then	
*
*	rotacao em torno do eixo 3 (r3=alfa deEuler), ou seja, de x para y
	iflag=0
	cr3=cdir(1) 
	sr3=cdir(2) 
	if(dabs(dl(2)).lt.tol) then  !// XZ
		iflag=1
		cr3=+1.0d0
		sr3=+0.0d0
	endif
*
*	rotacao em torno do eixo 2 (r2=gama de Euler), ou seja, de z para x
	cr2=cdir(1) 
	sr2=-cdir(3) 
	if((dabs(dl(3)).lt.tol).and.(iflag.eq.0)) then	!//XY
		cr2=1.0d0
		sr2=0.0d0
		iflag=2
	endif
*	
	else  !caso estrutura // ZY
		cr3=0.0d0
		sr3=-1.0d0
		cr2=cdir(2)
		sr2=cdir(3)
		iflag=3
	endif
*
	rot3(1,1)=cr3
	rot3(1,2)=sr3
	rot3(1,3)=0.0d0
	rot3(2,1)=-rot3(1,2)
	rot3(2,2)=rot3(1,1)
	rot3(2,3)=0.0d0
	rot3(3,1)=0.0d0
	rot3(3,2)=0.0d0
	rot3(3,3)=1.0d0
*
	rot2(1,1)=cr2
	rot2(1,2)=0.0d0
	rot2(1,3)=-sr2
	rot2(2,1)=0.0d0
	rot2(2,2)=1.0d0
	rot2(2,3)=0.0d0
	rot2(3,1)=-rot2(1,3)
	rot2(3,2)=0.0d0
	rot2(3,3)=rot2(1,1)
*
*	rotacao em torno do eixo 1 (r1=beta de Euler), ou seja, giro da seção transversal
	cr1=dcosd(agl)
	sr1=dsind(agl)
	rot1(1,1)=1.0d0
	rot1(1,2)=0.0d0
	rot1(1,3)=0.0d0
	rot1(2,1)=0.0d0
	rot1(2,2)=cr1
	rot1(2,3)=sr1
	rot1(3,1)=0.0d0
	rot1(3,2)=-rot1(2,3)
	rot1(3,3)=rot1(2,2)
*
*	executa rot=rot1*rot2*rot3
	call multmm(nd,rot2,rot3,rot)
	call multmm(nd,rot1,rot,rot2)
*
*
	do i=1,nge
		mrot(i,i)=0.0d0
		do j=i+1,nge
			mrot(i,j)=0.0d0
			mrot(j,i)=0.0d0
		enddo
	enddo
	do i=1,4  !a matriz do elemento tem 4x4 blocos
		ib=(i-1)*nd
		do j=1,nd
			do k=1,nd
				mrot(j+ib,k+ib)=rot2(j,k)
			enddo
		enddo
	enddo
*
	deallocate (dl,cdir,rot,rot1,rot2,rot3,cf,ci)
	return
	end
*-------------------------------------------------------------------------------        
*	subrotina que calcula matrizes de rigidez e massa local do portico 3D
*	u1=(u1,u2,u3,r1,r2,r3) onde: u1 e r1 axial e torcao
*	https://previa.uclm.es/profesorado/evieira/asignatura/meccomp/book/estructuras/vigas/Porticos_espaciales.htm
*	http://what-when-how.com/the-finite-element-method/fem-for-frames-finite-element-method-part-1/
*	
*
	subroutine Lmatriz(nd,nge,A,L,E,G,I2,I3,It,ro,kel,mel,mrot,tol)
	implicit none
*
	integer nd,nge,i,j,imas 
	real*8 tol,tolk,tolm	
	real*8 c1,c2,c3,c4,c5,sm	
	real*8 L,A,I2,I3,I1,It,E,G,ro  
*
	real*8 kel(nge,nge),mel(nge,nge),mrot(nge,nge) 
	real*8, allocatable :: mtrot(:,:),mkr(:,:) 
*
*	determinacao da matriz de rigidez 
c	write(*,*) '  ......matriz de rigidez' 
	do i=1,nge
		kel(i,i)=0.0d0
		do j=i+1,nge
			kel(i,j)=0.0d0
			kel(j,i)=0.0d0
		enddo
	enddo

	c1=E*A/L
	c2=12*E/(L**3)
	c3=6*E/(L**2)
	c4=G/L
	c5=4*E/L

*
*	linha 1	
	kel(1,1)=c1
	kel(1,7)=-c1
*
*	linha 2
	kel(2,2)=c2*I3
	kel(2,6)=c3*I3
	kel(2,8)=-c2*I3
	kel(2,12)=c3*I3
*
*	linha 3
	kel(3,3)=c2*I2
	kel(3,5)=-c3*I2
	kel(3,9)=-c2*I2
	kel(3,11)=-c3*I2
*
*	linha 4
	kel(4,4)=c4*It
	kel(4,10)=-c4*It
*
*	linha 5
	kel(5,5)=c5*I2
	kel(5,9)=c3*I2
	kel(5,11)=(c5*I2)/2.0d0
*
*	linha 6
	kel(6,6)=c5*I3
	kel(6,8)=-c3*I3
	kel(6,12)=kel(6,6)/2.0d0
*
*	linha 7
	kel(7,7)=kel(1,1)
*
*	linha 8
	kel(8,8)=c2*I3
	kel(8,12)=-c3*I3
*
*	linha 9
	kel(9,9)=c2*I2
	kel(9,11)=c3*I2
*
*	linha 10
	kel(10,10)=c4*It
*
*	linha 11
	kel(11,11)=c5*I2
*
*	linha 12
	kel(12,12)=c5*I3
*	
*	simetria i=linha j=coluna
	
	do j=1,nge
		do i=1,nge
			kel(i,j)=kel(j,i)	
		enddo
	enddo

	
*
*	determinacao da matriz de massa 
c	write(*,*) '  ......matriz de massa' 
	do i=1,nge
		mel(i,i)=0.0d0
		do j=i+1,nge
			mel(i,j)=0.0d0
			mel(j,i)=0.0d0
		enddo
	enddo
	I1=It
	c1=ro*A*L/420.0d0
	c2=c1*(L/2.0d0)
	c3=c1*140.0d0*I1/A
	c4=c1*(L/2.0d0)**2

*
*	linha 1		 
	mel(1,1)= 140.0d0*c1
	mel(1,7)=mel(1,1)/2.0d0
*
*	linha 2
	mel(2,2)= 156.0d0*c1
	mel(2,6)= 044.0d0*c2
	mel(2,8)= 054.0d0*c1
	mel(2,12)=-26.0d0*c2
*
*	linha 3
	mel(3,3)=mel(2,2)
	mel(3,5)=-mel(2,6)
	mel(3,9)=mel(2,8)
	mel(3,11)=-mel(2,12)
*
*	linha 4
	mel(4,4)=c3
	mel(4,10)=-mel(4,4)/2.0d0
*
*	linha 5
	mel(5,5)=16.0d0*c4
	mel(5,9)=mel(2,12)
	mel(5,11)=-12.0d0*c4
*
*	linha 6
	mel(6,6)=mel(5,5)
	mel(6,8)=mel(3,11)
	mel(6,12)=mel(5,11)
*
*	linha 7
	mel(7,7)=mel(1,1)
*
*	linha 8
	mel(8,8)=mel(2,2)																
	mel(8,12)=mel(3,5)
*
*	linha 9
	mel(9,9)=mel(2,2)
	mel(9,11)=mel(2,6)
*
*	linha 10
	mel(10,10)=mel(4,4)
*
*	linha 11
	mel(11,11)=mel(5,5)
*
*	linha 12
	mel(12,12)=mel(6,6)
*	
*	simetria i=linha j=coluna
	do j=1,nge
		do i=j+1,nge
			mel(i,j)=mel(j,i)
		enddo
	enddo	
*
*	calcula a matriz de massa diagonal caso imas.ne.0
*	HRZ Lumping Method
*	http://kis.tu.kielce.pl/mo/COLORADO_FEM/colorado/IFEM.Ch31.pdf
*	VERIFICAR se aplica
	imas=0
	if(imas.ne.0) then
		sm=0.0d0
		do i=1,nd
			sm=sm+mel(i,i)
		enddo
		do i=1,nge
			mel(i,i)=mel(i,i)/sm
			do j=i+1,nge
				mel(i,j)=0.0d0
				mel(j,i)=0.0d0
			enddo
		enddo
	endif
*
*	rotacao matriz de rigidez
	allocate(mtrot(nge,nge),mkr(nge,nge))
*	obtento transposta da rotacao
c	write(*,*) '  ......rotacao matriz de rigidez' 
	do i=1,nge
		mtrot(i,i)=mrot(i,i)
		do j=i+1,nge
			mtrot(i,j)=mrot(j,i)
			mtrot(j,i)=mrot(i,j)
		enddo
	enddo
	call multmm(nge,kel,mrot,mkr)
	call multmm(nge,mtrot,mkr,kel)


*	rotacao matriz de massa
c	write(*,*) '  ......rotacao matriz de massa' 
	call multmm(nge,mel,mrot,mkr)
 	call multmm(nge,mtrot,mkr,mel)
*
*	verificacoes diagonal e simetria
	do i=1,nge
		if(kel(i,i).le.tol) then
			write(*,*) i,kel(i,i)
			stop 'ERRO981 diagonal de kel'
		endif
		if(mel(i,i).le.tol) then
			write(*,*) i,mel(i,i)
			stop 'ERRO99 diagonal de mel'
		endif
		do j=i+1,nge
			tolk=dabs(kel(i,j)-kel(j,i))

			if(tolk.gt.tol) then   !zerar residuo?	 'roubo'
			kel(j,i)= 0
			kel(j,i)=kel(i,j)
			endif

			tolk=dabs(kel(i,j)-kel(j,i))

			if(tolk.gt.tol) then
				write(*,*) i,j,tolk,kel(i,j),kel(j,i)
				stop 'ERRO1 simetria de kel'
			endif


			tolm=dabs(mel(i,j)-mel(j,i))
			if(tolm.gt.tol) then
				write(*,*) i,j,tolm
				stop 'ERRO1 simetria de mel'
			endif
		enddo
	enddo
	deallocate(mtrot,mkr)
	return
	end
*-------------------------------------------------------------------------------        
*	subrotina que efetua montagem das matrizes globais
	subroutine Gmatriz(ngn,nge,nne,ng,ne,ie,econ,kel,mel,kgl,mgl,tol)
	implicit none
*								 
	integer ngn,nge,nne,ne,ng,ie,ni,nf   
	integer i,j,ine,jne,ngi,ngj,ii,jj,i1,j1,iil,jjl,i2,j2
	integer econ(ne,nne)
*
	real*8 tol,tolk,tolm	
	real*8 kel(nge,nge),mel(nge,nge)		 
	real*8 kgl(ng,ng),mgl(ng,ng)  
*
	if(ie.eq.1) then						!ie conta o numero do elemento
		do j=1,ng							!zerando todas as matrizes globais
			kgl(j,j)=0.0d0
			mgl(j,j)=0.0d0
			do i=j+1,ng
				kgl(i,j)=0.0d0
				mgl(i,j)=0.0d0
				kgl(j,i)=kgl(i,j)
				mgl(j,i)=mgl(i,j)
			enddo
		enddo
	endif

*						!ni recebe o numero da conectividade do elemento ie do no1
	ni=econ(ie,1)		!nf recebe o no final do elemento "ie"
	nf=econ(ie,2)
c	write(*,*) '        ...espalhamento na matriz global',ie ,ni,nf
	do ine=1,nne		!contador que vai de 1 até numero de nos do elemento
		ngi=econ(ie,ine) !ngi rcebe o numero da conectividade do elemento ie
		ii=ngn*(ngi-1)	!ii recebe valor multiplo de 6, variand de 0 ate ngl
		iil=ngn*(ine-1)	!iil recebe 0 ou 6
		do i=1,ngn		 !roda de 1 ao 6
			i1=ii+i		 !recebe multiplo de6+(1ate6)
			i2=iil+i	 !recebe (0ou6)+(0ou6) valor de1a12

			do jne=1,nne
				ngj=econ(ie,jne)!ngj rcebe o n da conectividade do elemento ie
				jj=ngn*(ngj-1)!jj rcebe valor multilo de 6,variand de0 ate ngl
				jjl=ngn*(jne-1)	 !jjl recebe 0 ou 6
				do j=1,ngn	 !1a6
					j1=jj+j	   !recebe multiplo de6+(1ate6)
					j2=jjl+j   !recebe (0ou6)+(0ou6)valor de1a12
					kgl(i1,j1)=kgl(i1,j1)+kel(i2,j2)
					mgl(i1,j1)=mgl(i1,j1)+mel(i2,j2)
				enddo
			enddo

		enddo
	enddo
*
	if(ie.eq.ne) then
		do i=1,ng
			if(kgl(i,i).le.tol) then
				write(*,*) 'ERRO diagonal de kgl',i
				stop 
			endif
			if(mgl(i,i).le.tol) then
				write(*,*)  'ERRO diagonal de mgl',i
				stop 
			endif
			do j=i+1,ng
				tolk=dabs(kgl(i,j)-kgl(j,i))

				if(tolk.ge.tol) then	!igualar residuo? "roubo"
				!write(*,*) kgl(j,i),kgl(i,j),i,j
				kgl(j,i)=kgl(i,j)
				endif

				tolk=dabs(kgl(i,j)-kgl(j,i))

				if(tolk.ge.tol) then
					write(*,*) 'ERRO simetria de kgl',i,j,tolk,tol
					stop 
				endif
				tolm=dabs(mgl(i,j)-mgl(j,i))
				if(tolm.ge.tol) then
					write(*,*) 'ERRO simetria de mgl',i,j,tolm,tol
					stop 
				endif
			enddo
		enddo
	endif
*
	return
	end
*-------------------------------------------------------------------------------        
*	subrotina que introduz condicoes de contorno
	subroutine Contorno(ng,nnp,ngn,rig,fex,kgl,fgl,tol,nvarexe,
	1ivar,mgl)
	implicit none
*
	integer ng,nnp,ngn   
	integer i,j,i1,j1,nvarexe,ivar  
	real*8 tol,tolk	
	real*8 rig(nnp,ngn),fex(nnp,ngn)
	real*8 kgl(ng,ng),fgl(ng)
	real*8 mgl(ng,ng)  
*
*	introducao da rigides dos apoios
	do i=1,nnp
		i1=ngn*(i-1)
		do j=1,ngn
			j1=i1+j
			kgl(j1,j1)=kgl(j1,j1)+rig(i,j)
	!write(*,*)	 rig(i,j,ivar), i,j,ivar
		enddo
	enddo

	!do i=1,ngn*nnp
	!	if(i.eq.1) write(*,1588) (j,j=1,ngn*nnp)
	!	write(*,1587) i,(kgl(i,j),j=1,ngn*nnp)
	!	write(*,*)
	!enddo

	
1587     format(1x,i2,1x,e8.2e2,1x,e8.2e2,1x,e8.2e2,1x,e8.2e2,1x,e8.2e2,1x,
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
*
*	introducao das forcas externas
	do i=1,ng
		fgl(i)=0.0d0
	enddo
	do i=1,nnp
		i1=ngn*(i-1)
		do j=1,ngn
			j1=i1+j
			fgl(j1)=fgl(j1)+fex(i,j)
			mgl(j1,j1)=mgl(j1,j1)+fex(i,j)/9.80665
		enddo
	enddo
*
	do i=1,ng
		if(kgl(i,i).le.tol) then
			write(*,*) 'ERRO diagonal de kgl',i
			stop 
		endif
		do j=i+1,ng
			tolk=dabs(kgl(i,j)-kgl(j,i))
			if(tolk.ge.tol) then
				write(*,*) 'ERRO simetria de kgl',i,j,tolk
				stop 
			endif
		enddo
	enddo

cc	open(10,file='temp.txt',status='replace')
cc	do i=1,ng
cc		write(10,1) i,((kgl(i,j)),j=1,ng)
cc	enddo
cc	close(10)
cc1	format(I5,<ng>f13.6)
	return
	end
*-------------------------------------------------------------------------------        
*	subrotina que executa multiplicacao matriz a pelo vetor b
*     obs: multiplicacao por colunas para aumentar performance   
	subroutine multmv(n,a,b,c)
	implicit none
*
	integer n,i,j 
	real*8  a(n,n),b(n),c(n),aux 

	aux=b(1)
	do i=1,n
		c(i)=a(i,1)*aux
	enddo
	do j=2,n
		aux=b(j)
		do i=1,n
			c(i)=c(i)+a(i,j)*aux
		enddo
	enddo
*
	return
	end
*-------------------------------------------------------------------------------        
*	subrotina que executa multiplicacao matriz a pela matriz b	  (MATRIZES REAIS)
*     obs: multiplicacao por colunas para aumentar performance   
	subroutine multmm(n,a,b,c)
	implicit none
*
	integer n,i,j,k 
	real*8  a(n,n),b(n,n),c(n,n),aux 
*
	do k=1,n
		aux=b(1,k)
		do i=1,n
			c(i,k)=a(i,1)*aux
		enddo
		do j=2,n
			aux=b(j,k)
			do i=1,n
				c(i,k)=c(i,k)+a(i,j)*aux
			enddo
		enddo
	enddo
*
	return
	end
*-------------------------------------------------------------------------------        
*	subrotina que executa multiplicacao matriz a pela matriz b	(MATRIZES COMPLEXAS)
*     obs: multiplicacao por colunas para aumentar performance   
	subroutine cmultmm(n,a,b,c)
	implicit none
*
	integer n,i,j,k 
	real*8  a(n,n),b(n,n),aux 
	double complex  c(n,n) 
*
	do k=1,n
		aux=b(1,k)
		do i=1,n
			c(i,k)=dcmplx(a(i,1)*aux,0.0d0)
		enddo
		do j=2,n
			aux=b(j,k)
			do i=1,n
				c(i,k)=c(i,k)+dcmplx(a(i,j)*aux,0.0d0)
			enddo
		enddo
	enddo
*
	return
	end
*-------------------------------------------------------------------------------        
*	subrotina que executa a transformacao K-shift*M
*     resultado armazenado em c, preservando K e M   
	subroutine mshift(n,k,m,c,shift)
	implicit none
*
	integer n,i,j 
	real*8  k(n,n),m(n,n),c(n,n),shift 
*
*	Executa operação preservando M e M
	do i=1,n
		do j=1,n
			c(i,j)=k(i,j)-shift*m(i,j)
		enddo
	enddo
*
	return
	end

*-------------------------------------------------------------------------------        
*	subrotina que executa produto escalar entre vetores a e b
*	prod=|a!*!b!*cos(r), cdir=cos(r)
	subroutine pescalar(n,a,b,prod,cdir)
	implicit none
*
	integer n,i
	real*8  prod,cdir,anorm,bnorm,pi,tol  
	real*8  a(n),b(n)  
*
	pi=4.0d0*datan(1.0d0)
	tol=1.0d-99
*
	anorm=a(1)**2
	bnorm=b(1)**2
	prod=a(1)*b(1)
	do i=2,n
		anorm=anorm+a(i)**2
		bnorm=bnorm+b(i)**2
		prod=prod+a(i)*b(i)
	enddo
	if((anorm.le.tol).or.(anorm.le.tol).or.(prod.le.tol)) then
		prod=0.0d0
		cdir=pi
	else
		cdir=prod/(anorm*bnorm)
	endif
*
	return
	end
