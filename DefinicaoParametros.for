*-------------------------------------------------------------------------------        
*	subrotina que executa parametros 
*
*
	subroutine Parametros(casovar,nse,nnp,ne,msim,
     1A,ro,I2,I3,It,E,G,angle,nce,ncn,ngn,rig,fex,ient,iout
     2,nvarA,nvarI2,nvarI3,nvarIt,nvarE,nvarG,nvarro,aux3,nvarkmax,
     3nvarrigx,nvarrigy,nvarrigz,nvarrigrx,nvarrigry,nvarrigrz,norig,
     4nofex,nvarfx,nvarfy,nvarfz,nvarmx,nvarmy,nvarmz,ivar,nvarexe,
     5sa,sI2,sI3,sit,se,sg,sro,sp,c1,agl,nprop,nlin,ncol,casoexe,
     6nseva,nsevi2,nsevi3,nseve,nsevro,asim,i2sim,i3sim,esim,rosim)

	implicit none


	integer nse,nnp,ne
	integer nne,nd,nde,ngn   
	integer ncn,nce,msim,ncc
	integer asim,i2sim,i3sim,esim,rosim
	integer nseva,nsevi2,nsevi3,nseve,nsevro
	integer i,j,ie,k,m,n,ise,in,je,ni,nf,ii,jj,eq,vv
	integer x,y,z,j21,iii,k1,k2,k3
	integer nvarA, nvarI2, nvarI3, nvarIt, nvarG, nvarE, nvarro
	integer nvarexe,nvarkmax
	integer nvarrigx,nvarrigy,nvarrigz,nvarrigrx,nvarrigry,nvarrigrz
	integer nvarfx,nvarfy,nvarfz,nvarmx,nvarmy,nvarmz
	integer ivar,nprop,iv 
	integer ivargeral,aux3
	integer ivarrigx,ivarrigy,ivarrigz,ivarrigrx,ivarrigry,ivarrigrz
	integer ivarfx,ivarfy,ivarfz,ivarmx,ivarmy,ivarmz
	integer ivarA,ivarE,ivarro,ivarI2,ivarI3
	integer varrigrx,varrigry,varrigrz
	integer varrigx,varrigy,varrigz
	integer varfx,varfy,varfz,varmx,varmy,varmz
	integer varA,varE,varI2,varI3,varro
	integer ient,iout	
	integer npropvar,nlin,ncol
	integer casovar, casoexe
	integer fexx,fexy,fexz,fexmx,fexmy,fexmz
	integer j12
	real*8 c1

	real*8, allocatable :: vif(:),vsf(:),dvf(:)!var inf,sup/incr de var molas
	real*8, allocatable :: vifex(:),vsfex(:),dvfex(:)!var inf,sup/incr.var for
	real*8, allocatable :: vetaux(:) !vetor auxiliar 
	real*8, allocatable :: dl(:),vaux(:,:),vauxfex(:)

	integer norig(nce)
	integer nofex(ncn)
	integer iA(99),iyoung(99),iI2(99),iI3(99),iro(99)
	integer rigrx(99),rigry(99),rigrz(99)
	integer rigx(99),rigy(99),rigz(99)
	integer elea(9999),elei2(9999),elei3(9999)
	integer elee(9999),elero(9999)
	integer aux4(99)
	real*8 rig(nnp,ngn)
	real*8 L(ne),A(ne),I2(ne)
	real*8 I3(ne),It(ne),E(ne)
	real*8 G(ne),ro(ne) 	
	real*8 fex(nnp,ngn)
 	real*8 sa(ne),sI2(ne),sI3(ne),sit(ne),se(ne),sg(ne)
	real*8 sro(ne),sp(ne),agl(ne),angle(ne)
	real*8 vsa(9999),via(9999)
	real*8 vsi2(9999),vii2(9999)
	real*8 vsi3(9999),vii3(9999)
	real*8 vse(9999),vie(9999)
	real*8 vsro(9999),viro(9999)
	real*8 dva(9999),dvi2(9999),dvi3(9999)
	real*8 dve(9999),dvro(9999)	


	if(ivar.eq.1)then
		allocate(vetaux(nprop))
	endif

				
	   !leitura das condicoes de contorno (molas)
	if(ivar.eq.1)then
	allocate(vif(ngn),vsf(ngn),dvf(ngn),vaux(ngn,nce))

		do i=1,nnp
			do j=1,ngn
				rig(i,j)=0.0d0
			enddo
		enddo

	do i=1,nce,1
	  norig(i)=0
	enddo
	endif


	if(ivar.eq.1)then

	do i=1,nce,1
		write(iout,06) 
		read(ient,*) norig(i),(vaux(j,i),j=1,ngn)
		write(iout,05) norig(i),(vaux(k,i),k=1,ngn)								   
		read(ient,*) (vsf(n),n=1,ngn)
		read(ient,*) (vif(m),m=1,ngn)
	enddo
	endif

	do ivarrigx =1,nvarrigx
			iv=1
				if(nvarrigx.le.1) then
					dvf(1)=1.0d0
				else
					dvf(1)=(vsf(1)-vif(1))/dfloat(nvarrigx-1)
				endif

				if (dvf(1).lt.0) stop 'ERRO, variacoes negativas x'
	enddo
	do ivarrigy =1,nvarrigy
			iv=2 
				if(nvarrigy.le.1) then
					dvf(2)=1.0d0
				else
					dvf(2)=(vsf(2)-vif(2))/dfloat(nvarrigy-1)
				endif
				if (dvf(2).lt.0) stop 'ERRO, variacoes negativas1'
	enddo
	do ivarrigz =1,nvarrigz
			iv=3 
				if(nvarrigz.le.1) then
					dvf(3)=1.0d0
				else
					dvf(3)=(vsf(3)-vif(3))/dfloat(nvarrigz-1)
				endif
				if (dvf(3).lt.0) stop 'ERRO, variacoes negativas1'
	enddo
	do ivarrigrx =1,nvarrigrx
			iv=4 
				if(nvarrigrx.le.1) then
					dvf(4)=1.0d0
				else
					dvf(4)=(vsf(4)-vif(4))/dfloat(nvarrigrx-1)
				endif
				if (dvf(4).lt.0) stop 'ERRO, variacoes negativas1'
	enddo
	do ivarrigry =1,nvarrigry
			iv=5 
				if(nvarrigry.le.1) then
					dvf(5)=1.0d0
				else
					dvf(5)=(vsf(5)-vif(5))/dfloat(nvarrigry-1)
				endif
				if (dvf(5).lt.0) stop 'ERRO, variacoes negativas1'
	enddo
	do ivarrigrz =1,nvarrigrz
			iv=6 
				if(nvarrigrz.le.1) then
					dvf(6)=1.0d0
				else
					dvf(6)=(vsf(6)-vif(6))/dfloat(nvarrigrz-1)
				endif
				if (dvf(6).lt.0) stop 'ERRO, variacoes negativas1'

	enddo
 


	!leitura das condicoes de contorno (forças)
 	if(ncn.eq.0) then
		goto 19
	endif


	if(ivar.eq.1)then
	allocate(vifex(ngn),vsfex(ngn),dvfex(ngn),vauxfex(ngn))


		do i=1,nnp
			do j=1,ngn
				fex(i,j)=0.0d0
			enddo
		enddo

	do i=1,ncn,1
	nofex(i)=0
	enddo

 	endif


	if(ivar.eq.1)then
	do i=1,ncn,1
		write(iout,10) 
		read(ient,*) nofex(i),(vauxfex(j),j=1,ngn)
		write(iout,05) nofex(i),(vauxfex(k),k=1,ngn)								   
		read(ient,*) (vsfex(n),n=1,ngn)
		read(ient,*) (vifex(m),m=1,ngn)
	enddo
	endif

	do ivarfx =1,nvarfx
			iv=1
				if(nvarfx.le.1) then
					dvfex(1)=1.0d0
				else
					dvfex(1)=(vsfex(1)-vifex(1))/dfloat(nvarfx-1)
				endif
				if (dvfex(1).lt.0) stop 'ERRO, variacoes negativas x'
	enddo
	do ivarfy =1,nvarfy
			iv=2 
				if(nvarfy.le.1) then
					dvfex(2)=1.0d0
				else
					dvfex(2)=(vsfex(2)-vifex(2))/dfloat(nvarfy-1)
				endif
				if (dvfex(2).lt.0) stop 'ERRO, variacoes negativas1'
	enddo
	do ivarfz =1,nvarfz
			iv=3 
				if(nvarfz.le.1) then
					dvfex(3)=1.0d0
				else
					dvfex(3)=(vsfex(3)-vifex(3))/dfloat(nvarfz-1)
				endif
				if (dvfex(3).lt.0) stop 'ERRO, variacoes negativas1'
		enddo
	  do ivarmx =1,nvarmx
			iv=4 
				if(nvarmx.le.1) then
					dvfex(4)=1.0d0
				else
					dvfex(4)=(vsfex(4)-vifex(4))/dfloat(nvarmx-1)
				endif
				if (dvfex(4).lt.0) stop 'ERRO, variacoes negativas1'
	  enddo
	  do ivarmy =1,nvarmy
			iv=5 
				if(nvarmy.le.1) then
					dvfex(5)=1.0d0
				else
					dvfex(5)=(vsfex(5)-vifex(5))/dfloat(nvarmy-1)
				endif
				if (dvfex(5).lt.0) stop 'ERRO, variacoes negativas1'
	  enddo
	  do ivarmz =1,nvarmz
			iv=6 
				if(nvarmz.le.1) then
					dvfex(6)=1.0d0
				else
					dvfex(6)=(vsfex(6)-vifex(6))/dfloat(nvarmz-1)
				endif
				if (dvfex(6).lt.0) stop 'ERRO, variacoes negativas1'
	  enddo
  19  continue


	if(ivar.eq.1) then
		
		if(nseva.ge.1) then
				read(ient,*) (elea(ise),ise=1,nseva)
				read(ient,*) (vsa(elea(ise)),ise=1,nseva)
				read(ient,*) (via(elea(ise)),ise=1,nseva)
		endif

		if(nsevi2.ge.1) then
				read(ient,*) (elei2(ise),ise=1,nsevi2)
				read(ient,*) (vsi2(elei2(ise)),ise=1,nsevi2)
				read(ient,*) (vii2(elei2(ise)),ise=1,nsevi2)
		endif

		if(nsevi3.ge.1) then
				read(ient,*) (elei3(ise),ise=1,nsevi3)
				read(ient,*) (vsi3(elei3(ise)),ise=1,nsevi3)
				read(ient,*) (vii3(elei3(ise)),ise=1,nsevi3)
		endif
		
		if(nseve.ge.1) then
				read(ient,*) (elee(ise),ise=1,nseve)
				read(ient,*) (vse(elee(ise)),ise=1,nseve)
				read(ient,*) (vie(elee(ise)),ise=1,nseve)
		endif

		if(nsevro.ge.1) then
				read(ient,*) (elero(ise),ise=1,nsevro)
				read(ient,*) (vsro(elero(ise)),ise=1,nsevro)
				read(ient,*) (viro(elero(ise)),ise=1,nsevro)
		endif

	endif

	if(ivar.eq.1) then

		do ise=1,nseva
			if(nvarA.le.1) then
				dva(elea(ise))=1.0d0
			else
				dva(elea(ise))=vsa(elea(ise))-via(elea(ise))
				dva(elea(ise))=dva(elea(ise))/dfloat(nvarA-1)
			endif
			if (dva(elea(ise)).lt.0) stop 'ERRO, variacoes negativas'
		enddo

		do ise=1,nsevi2
			if(nvarI2.le.1) then
				dvi2(elei2(ise))=1.0d0
			else
				dvi2(elei2(ise))=vsi2(elei2(ise))-vii2(elei2(ise))
				dvi2(elei2(ise))=dvi2(elei2(ise))/dfloat(nvarI2-1)
			endif
			if (dvi2(elei2(ise)).lt.0) stop 'ERRO, variacoes negativas'
		enddo

		do ise=1,nsevi3
			if(nvarI3.le.1) then
				dvi3(elei3(ise))=1.0d0
			else
				dvi3(elei3(ise))=vsi3(elei3(ise))-vii3(elei3(ise))
				dvi3(elei3(ise))=dvi3(elei3(ise))/dfloat(nvarI3-1)
			endif
			if (dvi3(elei3(ise)).lt.0) stop 'ERRO, variacoes negativas'
		enddo

		do ise=1,nseve
			if(nvarE.le.1) then
				dve(elee(ise))=1.0d0
			else
				dve(elee(ise))=vse(elee(ise))-vie(elee(ise))
				dve(elee(ise))=dve(elee(ise))/dfloat(nvarE-1)
			endif
			if (dve(elee(ise)).lt.0) stop 'ERRO, variacoes negativas'
		enddo

		do ise=1,nsevro
			if(nvarro.le.1) then
				dvro(ise)=1.0d0
			else
				dvro(elero(ise))=vsro(elero(ise))-viro(elero(ise))
				dvro(elero(ise))=dvro(elero(ise))/dfloat(nvarro-1)
			endif
			if (dvro(elero(ise)).lt.0) stop 'ERRO, variacoes negativas'
		enddo

	endif

	!variacao das propriedades
	ie=0
	i=0

	!Rotina que gera todas as combinacoes de variacoes
	if(casoexe.eq.1) then

	ie=0
	j=1

	if(ivar.eq.1) then
		do i=1,nse
			iA(i)=1
			iyoung(i)=1
			iI2(i)=1
			iI3(i)=1
			iro(i)=1
		enddo
		do i=1,nce
			rigx(i)=1
			rigy(i)=1
			rigz(i)=1
			rigrx(i)=1
			rigry(i)=1
			rigrz(i)=1
		enddo
	endif

	do i=1,nse
		A(i)=0
		I2(i)=0
		I3(i)=0
		It(i)=0
		E(i)=0
		G(i)=0
		ro(i)=0
	enddo


	do i=1,nvarexe

		if(msim.eq.1) ncc=1
		if(msim.eq.0) ncc=nce
		!if(asim.eq.1) nseva=1
		!if(i2sim.eq.1) nsevi2=1
		!if(i3sim.eq.1) nsevi3=1
		!if(esim.eq.1) nseve=1
		!if(rosim.eq.1) nsevro=1

		if(ivar.ne.1)jj=jj+nvarrigrz
		if(jj.ge.nvarrigrz.or.nvarrigrz.eq.1)then
			jj=0
			if(nvarrigrz.ne.1) then
				if(ivar.ne.1)rigrz(1)=rigrz(1)+1
				if(rigrz(1).ne.(nvarrigrz+1)) goto 1133
			endif
		endif
		do j12=1,ncc,1
			if(rigrz(j12).eq.(nvarrigrz+1).or.nvarrigrz.eq.1)then
				if(j12.ne.ncc) rigrz(j12)=1
				if(ivar.ne.1.and.j12.ne.ncc) then
					if(ivar.ne.1)rigrz(j12+1)=rigrz(j12+1)+1
					if(rigrz(j12+1).ne.(nvarrigrz+1)) goto 1133
				endif
			endif
		enddo


		if(rigrz(ncc).eq.(nvarrigrz+1).or.nvarrigrz.eq.1) then
			rigrz(ncc)=1
			if(nvarrigry.ne.1) then
				if(ivar.ne.1)rigry(1)=rigry(1)+1
				if(rigry(1).ne.(nvarrigry+1)) goto 1133
			endif
		endif
		do j12=1,ncc,1
			if(rigry(j12).eq.(nvarrigry+1).or.nvarrigry.eq.1)then
				if(j12.ne.ncc)rigry(j12)=1
				if(ivar.ne.1.and.j12.ne.ncc) then
					if(ivar.ne.1)rigry(j12+1)=rigry(j12+1)+1
					if(rigry(j12+1).ne.(nvarrigrz+1)) goto 1133
				endif
			endif
		enddo


		if(rigry(ncc).eq.(nvarrigry+1).or.nvarrigry.eq.1) then
			rigry(ncc)=1
 			if(nvarrigrx.ne.1) then
				if(ivar.ne.1)rigrx(1)=rigrx(1)+1
				if(rigrx(1).ne.(nvarrigrx+1)) goto 1133
			endif
		endif
		do j12=1,ncc,1
			if(rigrx(j12).eq.(nvarrigrx+1).or.nvarrigrx.eq.1)then
				if(j12.ne.ncc)rigrx(j12)=1
				if(ivar.ne.1.and.j12.ne.ncc) then
					if(ivar.ne.1)rigrx(j12+1)=rigrx(j12+1)+1
					if(rigrx(j12+1).ne.(nvarrigrx+1)) goto 1133
				endif
			endif
		enddo


		if(rigrx(ncc).eq.(nvarrigrx+1).or.nvarrigrx.eq.1) then
			rigrx(ncc)=1
 			if(nvarrigz.ne.1) then
				if(ivar.ne.1)rigz(1)=rigz(1)+1
				if(rigz(1).ne.(nvarrigz+1)) goto 1133
			endif
		endif
		do j12=1,ncc,1
			if(rigz(j12).eq.(nvarrigz+1).or.nvarrigz.eq.1)then
				if(j12.ne.ncc)rigz(j12)=1
				if(ivar.ne.1.and.j12.ne.ncc) then
					if(ivar.ne.1)rigz(j12+1)=rigz(j12+1)+1
					if(rigz(j12+1).ne.(nvarrigz+1)) goto 1133
				endif
			endif
		enddo

		if(rigz(ncc).eq.(nvarrigz+1).or.nvarrigz.eq.1) then
			rigz(ncc)=1
 			if(nvarrigy.ne.1) then
				if(ivar.ne.1)rigy(1)=rigy(1)+1
				if(rigy(1).ne.(nvarrigy+1)) goto 1133
			endif
		endif
		do j12=1,ncc,1
			if(rigy(j12).eq.(nvarrigy+1).or.nvarrigy.eq.1)then
				if(j12.ne.ncc)rigy(j12)=1
				if(ivar.ne.1.and.j12.ne.ncc) then
					if(ivar.ne.1)rigy(j12+1)=rigy(j12+1)+1
					if(rigy(j12+1).ne.(nvarrigy+1)) goto 1133
				endif
			endif
		enddo
		
		if(rigy(ncc).eq.(nvarrigy+1).or.nvarrigy.eq.1) then
			rigy(ncc)=1
 			if(nvarrigx.ne.1) then
				if(ivar.ne.1)rigx(1)=rigx(1)+1
				if(rigx(1).ne.(nvarrigx+1)) goto 1133
			endif
		endif
		do j12=1,ncc,1
			if(rigx(j12).eq.(nvarrigx+1).or.nvarrigx.eq.1)then
				if(j12.ne.ncc)rigx(j12)=1
				if(ivar.ne.1.and.j12.ne.ncc) then
					if(ivar.ne.1)rigx(j12+1)=rigx(j12+1)+1
					if(rigx(j12+1).ne.(nvarrigx+1)) goto 1133
				endif
			endif
		enddo
		
		if(rigx(ncc).eq.(nvarrigx+1).or.nvarrigx.eq.1) then
			rigx(ncc)=1
 			if(nvarA.ne.1) then
				if(ivar.ne.1)iA(elea(1))=iA(elea(1))+1
				if(iA(elea(1)).ne.(nvarA+1)) goto 1133
			endif
		endif
		do j12=1,nseva,1
			if(iA(elea(j12)).eq.(nvarA+1).or.nvarA.eq.1)then
				if(j12.ne.nseva)iA(elea(j12))=1
				if(ivar.ne.1.and.j12.ne.nseva) then
					if(ivar.ne.1)iA(elea(j12+1))=iA(elea(j12+1))+1
					if(iA(elea(j12+1)).ne.(nvarA+1)) goto 1133
				endif
			endif
		enddo

		if(iA(elea(nseva)).eq.(nvarA+1).or.nvarA.eq.1) then
			iA(elea(nseva))=1
 			if(nvarI2.ne.1) then
				if(ivar.ne.1)iI2(elei2(1))=iI2(elei2(1))+1
				if(iI2(elei2(1)).ne.(nvarI2+1)) goto 1133
			endif
		endif
		do j12=1,nsevi2,1
			if(iI2(elei2(j12)).eq.(nvarI2+1).or.nvarI2.eq.1)then
				if(j12.ne.nsevi2)iI2(elei2(j12))=1
				if(ivar.ne.1.and.j12.ne.nsevi2) then
					if(ivar.ne.1)iI2(elei2(j12+1))=iI2(elei2(j12+1))+1
					if(iI2(elei2(j12+1)).ne.(nvarI2+1)) goto 1133
				endif
			endif
		enddo

		if(iI2(elei2(nsevi2)).eq.(nvarI2+1).or.nvarI2.eq.1) then
			iI2(elei2(nsevi2))=1
 			if(nvarI3.ne.1) then
				if(ivar.ne.1)iI3(elei3(1))=iI3(elei3(1))+1
				if(iI3(elei3(1)).ne.(nvarI3+1)) goto 1133
			endif
		endif
		do j12=1,nsevi3,1
			if(iI3(elei3(j12)).eq.(nvarI3+1).or.nvarI3.eq.1)then
				if(j12.ne.nsevi3)iI3(elei3(j12))=1
				if(ivar.ne.1.and.j12.ne.nsevi3) then
					if(ivar.ne.1)iI3(elei3(j12+1))=iI3(elei3(j12+1))+1
					if(iI3(elei3(j12+1)).ne.(nvarI3+1)) goto 1133
				endif
			endif
		enddo


		if(iI3(elei3(nsevi3)).eq.(nvarI3+1).or.nvarI3.eq.1) then
			iI3(elei3(nsevi3))=1
 			if(nvarE.ne.1) then
				if(ivar.ne.1)iyoung(elee(1))=iyoung(elee(1))+1
				if(iyoung(elee(1)).ne.(nvarE+1)) goto 1133
			endif
		endif
		do j12=1,nseve,1
			if(iyoung(elee(j12)).eq.(nvarE+1).or.nvarE.eq.1)then
	         if(j12.ne.nseve)iyoung(elee(j12))=1
			  if(ivar.ne.1.and.j12.ne.nseve) then
			   if(ivar.ne.1)iyoung(elee(j12+1))=iyoung(elee(j12+1))+1
			   if(iyoung(elee(j12+1)).ne.(nvarE+1)) goto 1133
			 endif
			endif
		enddo


		if(iyoung(elee(nseve)).eq.(nvarE+1).or.nvarE.eq.1) then
			iyoung(elee(nseve))=1
 			if(nvarro.ne.1) then
				if(ivar.ne.1)iro(elero(1))=iro(elero(1))+1
				if(iro(elero(1)).ne.(nvarro+1)) goto 1133
			endif
		endif
		do j12=1,nsevro,1
			if(iro(elero(j12)).eq.(nvarro+1).or.nvarro.eq.1)then
				if(j12.ne.nsevro)iro(elero(j12))=1
				if(ivar.ne.1.and.j12.ne.nsevro) then
					if(ivar.ne.1)iro(elero(j12+1))=iro(elero(j12+1))+1
					if(iro(elero(j12+1)).ne.(nvarro+1)) goto 1133
				endif
			endif
		enddo

1133				continue

				k2=0
				do k1=1,nce
					k2=k2+1
					if(msim.eq.1) k2=1
 					c1=vif(1)+(rigx(k2)-1)*dvf(1)
					rig(norig(k1),1)=vaux(1,k2)*c1
				enddo
				k2=0

				do k1=1,nce
					k2=k2+1
					if(msim.eq.1) k2=1
 					c1=vif(2)+(rigy(k2)-1)*dvf(2)
					rig(norig(k1),2)=vaux(2,k2)*c1
				enddo
				k2=0

				do k1=1,nce
					k2=k2+1
					if(msim.eq.1)k2=1
 					c1=vif(3)+(rigz(k2)-1)*dvf(3)
					rig(norig(k1),3)=vaux(3,k2)*c1
				enddo


				k2=0
				do k1=1,nce
					k2=k2+1
					if(msim.eq.1) k2=1
 					c1=vif(4)+(rigrx(k2)-1)*dvf(4)
					rig(norig(k1),4)=vaux(4,k2)*c1
				enddo
				k2=0

				do k1=1,nce
					k2=k2+1
					if(msim.eq.1) k2=1
 					c1=vif(5)+(rigry(k2)-1)*dvf(5)
					rig(norig(k1),5)=vaux(5,k2)*c1
				enddo
				k2=0

				do k1=1,nce
					k2=k2+1
					if(msim.eq.1)k2=1
 					c1=vif(6)+(rigrz(k2)-1)*dvf(6)
					rig(norig(k1),6)=vaux(6,k2)*c1
				enddo
				k2=0

				do k1=1,nse
					if(nseva.eq.0) then
						A(k1)=sA(k1)
					endif
				enddo
				if(nseva.ne.0) then
					do k1=1,nseva
					   k2=k2+1
					   if(asim.eq.1)k2=1						
					   c1=via(elea(k2))+(iA(elea(k2))-1)*dva(elea(k2))
					   A(elea(k1))=sA(elea(k2))*c1
					enddo
					do k1=1,nse
						if(A(k1).eq.0) A(k1)=sa(k1)
					enddo
				endif
				k2=0

				do k1=1,nse
					if(nsevi2.eq.0) then
						I2(k1)=si2(k1)
						It(k1)=I2(k1)
					endif
				enddo
				if(nsevi2.ne.0) then
					do k1=1,nsevi2
						k2=k2+1
						if(i2sim.eq.1)k2=1
				 c1=vii2(elei2(k2))+(iI2(elei2(k2))-1)*dvi2(elei2(k2))
					I2(elei2(k1))=si2(elei2(k2))*c1
					It(elei2(k1))=I2(elei2(k1))
					enddo
					do k1=1,nse
						if(I2(k1).eq.0) then
							I2(k1)=sI2(k1)
							It(k1)=sI2(k1)
						endif
					enddo
				endif
				k2=0
							
				do k1=1,nse
					if(nsevi3.eq.0) then
						I3(k1)=si3(k1)
						It(k1)=I3(k1)+It(k1)
					endif
				enddo
				if(nsevi3.ne.0) then
					do k1=1,nsevi3
					k2=k2+1
					if(i3sim.eq.1)k2=1
				 c1=vii3(elei3(k2))+(iI3(elei3(k2))-1)*dvi3(elei3(k2))
					I3(elei3(k1))=si3(elei3(k2))*c1
					It(elei3(k1))=I3(elei3(k1))+It(elei3(k1))
					enddo
					do k1=1,nse
						if(I3(k1).eq.0) then
							I3(k1)=sI3(k1)
							It(k1)=sIt(k1)+sI3(k1)
						endif
					enddo
				endif
				k2=0

				do k1=1,nse
					if(nseve.eq.0) then
						E(k1)=se(k1)
						G(k1)=se(k1)/(2*(1+sp(k1)))
					endif
				enddo
				if(nseve.ne.0) then
					do k1=1,nseve,1
					k2=k2+1
					if(esim.eq.1) k2=1
				   c1=vie(elee(k2))+(iyoung(elee(k2))-1)*dve(elee(k2))
					E(elee(k1))=se(elee(k2))*c1
					G(elee(k1))=se(elee(k2))*c1/(2*(1+sp(k1)))
					enddo
					do k1=1,nse
						if(E(k1).eq.0) then
							E(k1)=se(k1)
							G(k1)=se(k1)/(2*(1+sp(k1)))
						endif
					enddo
				endif
				k2=0

				do k1=1,nse
					if(nsevro.eq.0) then
						ro(k1)=sro(k1)
					endif
				enddo
				if(nsevro.ne.0) then
					do k1=1,nsevro
					k2=k2+1
					if(rosim.eq.1) k2=1
			     c1=viro(elero(k2))+(iro(elero(k2))-1)*dvro(elero(k2))
					ro(elero(k1))=sro(elero(k2))*c1
					enddo
					do k1=1,nse
						if(ro(k1).eq.0) then
							ro(k1)=sro(k1)
						endif
					enddo
				endif
				k2=0

 				goto 4043
	enddo
	
4043  continue
	!rotina para imprimir a sequencia das variacoes
	k2=1
	do k1=1,ncc
		if(nvarrigrz.gt.1) then
			aux4(k2)=rigrz(k1)
			k2=k2+1
		endif
		if(nvarrigry.gt.1) then
			aux4(k2)=rigry(k1)
			k2=k2+1
		endif
		if(nvarrigrx.gt.1) then
			aux4(k2)=rigrx(k1)
			k2=k2+1
		endif
		if(nvarrigz.gt.1) then
			aux4(k2)=rigz(k1)
			k2=k2+1
		endif
		if(nvarrigy.gt.1) then
			aux4(k2)=rigy(k1)
			k2=k2+1
		endif
		if(nvarrigx.gt.1) then
			aux4(k2)=rigx(k1)
			k2=k2+1
		endif
		if(msim.eq.1) goto 1690
	enddo

1690  continue

	if(nseva.gt.0) then
		do k1=1,nseva
			aux4(k2)=iA(elea(k1))
			k2=k2+1	
			if(asim.eq.1) goto 1691		
		enddo
	endif

1691  continue

	if(nsevi2.gt.0) then
		do k1=1,nsevi2
			aux4(k2)=iI2(elei2(k1))
			k2=k2+1
			if(i2sim.eq.1) goto 1692
		enddo
	endif

1692  continue

	if(nsevi3.gt.0) then
		do k1=1,nsevi3
			aux4(k2)=iI3(elei3(k1))
			k2=k2+1
			if(i3sim.eq.1) goto 1693
		enddo
	endif

1693	continue

	if(nseve.gt.0) then
		do k1=1,nseve
			aux4(k2)=iyoung(elee(k1))
			k2=k2+1
			if(esim.eq.1) goto 1694
		enddo
	endif

1694  continue

	if(nsevro.gt.0) then
		do k1=1,nsevro
			aux4(k2)=iro(elero(k1))
			k2=k2+1
			if(rosim.eq.1) goto 1695
		enddo
	endif

1695  continue

	iii=iii+1
	if(iii.eq.(2000/nse**2).or.ivar.gt.(nvarexe-9).or.ivar.lt.100)then
			write(*,987) (aux4(k3),k3=1,(k2-1))			
			iii=0
	endif

	!do k=1,nse
	!	write(*,*) A(k1),I2(k1),I3(k1)
	!enddo
!	if(msim.eq.0) then
!		iii=iii+1
!		if(iii.eq.100.or.ivar.gt.(nvarexe-100).or.ivar.lt.200) then
!			write(*,987)rigrz(1),rigrz(2),rigry(1),rigry(2),rigrx(1),
!     1		rigrx(2),rigz(1),rigz(2),rigy(1),rigy(2),rigx(1),rigx(2)			
!			iii=0
!	endif
!
!	if(msim.eq.1) then
!		iii=iii+1
!		if(iii.eq.100.or.ivar.gt.(nvarexe-100).or.ivar.lt.200) then
!			write(*,987)rigrz(1),rigry(1),rigrx(1),
!     1		rigz(1),rigy(1),rigx(1),iA(1),iA(2),iA(3),iA(4)			
!			iii=0
!		endif
!	endif

987	format(     1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,
	11x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,
	21x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)

	endif

01	format(5x,'Superel.',7x,'Area',4x,'Inerc Min',4x,'Inerc Max',4x, 
	1       'Inerc Tor',3x,'Mod. Elast',6x,'Poisson',3x,'Mass Espec',4x,
	2       'Rot Secao')
02	format(A11,<nprop>E13.6E2)
03	format(3x,'var',3x,'Prop.',3x,'Area',6x,'Inercia I2',3x,
     1'Inercia I3',2x,'Inerc Tor',3x,'Mod.El.L',4x,'Mod.El.T',3x,
     2'Mass Espec')

11   	format(9x,'Molas',3x,'NO',10x,'Kx',10x,'Ky',10x,'Kz',10x,'Krx',9x,
     1'Kry',9x,'Krz',5x)
13    format(i5,11x,i2,5x,6E12.5E2)
12    format(8x,'Forcas',3x,'NO',8x,'Fx',10x,'Fy',10x,'Fz',10x,'Mx.',9x,
     1'My.', 9x,'Mz.',5x)
   
04	format(A11,<ngn>E13.6E2)
05	format(5x,I4,<ngn>E13.6E2)
06	format(7x'No',6x,'Rig. Kx',6x,'Rig. Ky',6x,'Rig. Kz',5x
	1       'Rig. Krx',5x,'Rig. Kry',5x,'Rig. Krz')
07	 format(7x'No',6x,'Forca X',6x,'Forca Y',6x,'Forca Z',5x
	1       'Moment X',5x,'Moment Y',5x,'Moment Z')
08	format(/5x)
09	format(2x,'Var',I4,<ngn>E13.6E2)
10	format(7x'No',6x,'Fex. Fx',6x,'Fex. Fy',6x,'Fex. Fz',5x
	1       'Fex. Mx.',5x,'Fex. My.',5x,'Fex. Mz.')

	return
	end      