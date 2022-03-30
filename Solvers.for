c*********************************************************************72
c
cc CSVD computes the singular value decomposition of an M by N complex matrix.
c
c  Discussion:
c
c    This routine requires that N <= M.
c
c    The singular value decomposition of a complex M by N matrix A
c    has the form
c
c      A = U S V*
c
c    where 
c
c      U is an M by M unitary matrix,
c      S is an M by N diagonal matrix,
c      V is an N by N unitary matrix.
c
c    Moreover, the entries of S are nonnegative and occur on the diagonal
c    in descending order.
c
c    Thanks to Vladimir Rutsky for pointing out, on 10 January 2013,
c    that the line
c      IF ( N .LT. 2 ) GO TO 570
c    had been incorrectly modified to
c      if ( 2 .lt. n ) then
c    when
c      if ( 1 .lt. n ) then
c    was required.
c
c  Modified:
c
c    10 January 2013
c
c  Author:
c
c    Peter Businger, Gene Golub.
c    Modifications by John Burkardt.
c
c  Reference:
c
c    Peter Businger, Gene Golub,
c    Algorithm 358:
c    Singular Value Decomposition of a Complex Matrix,
c    Communications of the ACM,
c    Volume 12, Number 10, October 1969, pages 564-565.
c
c  Parameters:
c
c    Input/output, complex A(MMAX,*), the M by N matrix, which may be
c    augmented by P extra columns to which the transformation U*
c    is to be applied.  On output, A has been overwritten, and
c    if 0 < P, columns N+1 through N+P have been premultiplied by U*.
c
c    Input, integer M, N, the number of rows and columns in A.
c    It must be the case that 1 <= N <= M.  
c
c    Input, integer P, the number of vectors, stored in A(*,N+1:N+P),
c    to which the transformation U* should be applied.
c
c    Input, integer NU, the number of columns of U to compute.
c
c    Input, integer NV, the number of columns of V to compute.
c
c    Output, real S(N), the computed singular values.
c
c    Output, complex U(m,NU), the first NU columns of U.
c
c    Output, complex V(n,NV), the first NV columns of V.
c
c  Local Parameters:
c
c    Local, real ETA, the relative machine precision.
c    The original text uses ETA = 1.5E-8.
c
c    Local, real TOL, the smallest normalized positive number, divided by ETA.
c    The original test uses TOL = 1.E-31.
c
      subroutine dcsvd ( a, m, n, p, nu, nv, s, u, v, b, c, t )
      implicit none

      integer m
      integer n

      double complex a(m,n)
      double complex u(m,n)
      double complex v(n,n)
      real*8 s(n)

      real*8 b(m)
      real*8 c(m)
      real*8 t(m)
      
	real*8 cs
      real*8 eps
      real*8 eta
      real*8 f
      real*8 g
      real*8 h
      integer i
      integer j
      integer k
      integer kk
      integer k1
      integer l
      integer l1
      integer ll
      integer nu
      integer nv
      integer p
      double complex q
      double complex r
      real*8 sn
      real*8 tol
      real*8 w
      real*8 x
      real*8 y
      real*8 z

      save eta
      save tol

      data eta / 1.0e-030 /	!1.1920929e-017
      data tol / 1.0e-199 /		!1.5e-31
c
c  Householder reduction.
c
      c(1) = 0.0e+00
      k = 1

10    continue

      k1 = k + 1
c
c  Elimination of A(I,K), I = K+1, ..., M.
c
      z = 0.0e+00
      do i = k, m
        z = z + ( dreal ( a(i,k) ) )**2 + ( dimag ( a(i,k) ) )**2
      end do


      b(k) = 0.0e+00

      if ( tol .lt. z ) then

        z = dsqrt ( z )
        b(k) = z
        w = cdabs ( a(k,k) )

        if ( w .eq. 0.0e+00 ) then
          q = dcmplx ( 1.0e+00, 0.0e+00 )
        else
          q = a(k,k) / w
        end if

        a(k,k) = q * ( z + w )

        if ( k .ne. n + p ) then

          do j = k1, n + p

            q =d cmplx ( 0.0e+00, 0.0e+00 )
            do i = k, m
              q = q + dconjg ( a(i,k) ) * a(i,j)
            end do
            q = q / ( z * ( z + w ) )

            do i = k, m
              a(i,j) = a(i,j) - q * a(i,k)
            end do

          end do						   
c
c  Phase transformation.
c
          q = -dconjg ( a(k,k) ) / cdabs ( a(k,k) )

          do j = k1, n + p
            a(k,j) = q * a(k,j)
          end do

        end if

      end if
c
c  Elimination of A(K,J), J = K+2, ..., N
c
      if ( k .eq. n ) then
        go to 140
      end if

      z = 0.0e+00
      do j = k1, n
        z = z + ( dreal ( a(k,j) ) )**2 + ( dimag ( a(k,j) ) )**2
      end do

      c(k1) = 0.0e+00

      if ( tol .lt. z ) then

        z = dsqrt ( z )
        c(k1) = z
        w = cdabs ( a(k,k1) )

        if ( w .eq. 0.0e+00 ) then
          q = dcmplx ( 1.0e+00, 0.0e+00 )
        else
          q = a(k,k1) / w
        end if

        a(k,k1) = q * ( z + w )

        do i = k1, m

          q = dcmplx ( 0.0e+00, 0.0e+00 )

          do j = k1, n
            q = q + dconjg ( a(k,j) ) * a(i,j)
          end do

          q = q / ( z * ( z + w ) )

          do j = k1, n
            a(i,j) = a(i,j) - q * a(k,j)
          end do

        end do
c
c  Phase transformation.
c
        q = -dconjg ( a(k,k1) ) / cdabs ( a(k,k1) )
        do i = k1, m
          a(i,k1) = a(i,k1) * q
        end do

      end if

      k = k1
      go to 10
c
c  Tolerance for negligible elements.
c
140   continue

      eps = 0.0e+00
      do k = 1, n
        s(k) = b(k)
        t(k) = c(k)
        eps = dmax1 ( eps, s(k) + t(k) )
      end do

      eps = eps * eta
c
c  Initialization of U and V.
c
      if ( 0 .lt. nu ) then

       do j = 1, nu
          do i = 1, m
            u(i,j) = dcmplx ( 0.0e+00, 0.0e+00 )
          end do
          u(j,j) = dcmplx ( 1.0e+00, 0.0e+00 )
       end do

      end if

      if ( 0 .lt. nv ) then
       do j = 1, nv
          do i = 1, n
            v(i,j) = dcmplx ( 0.0e+00, 0.0e+00 )
          end do
          v(j,j) = dcmplx ( 1.0e+00, 0.0e+00 )
        end do

      end if
c
c  QR diagonalization.
c
      do kk = 1, n

        k = n + 1 - kk
c
c  Test for split.
c
220     continue

        do ll = 1, k

          l = k + 1 - ll

          if ( dabs ( t(l) ) .le. eps ) then
            go to 290
          end if

          if ( dabs ( s(l-1) ) .le. eps ) then
            go to 240
          end if

        end do
c
c  Cancellation of E(L).
c
240     continue

        cs = 0.0e+00
        sn = 1.0e+00
        l1 = l - 1

        do i = l, k

          f = sn * t(i)
          t(i) = cs * t(i)

          if ( dabs ( f ) .le. eps ) then
            go to 290
          end if

          h = s(i)
          w = dsqrt ( f * f + h * h )
          s(i) = w
          cs = h / w
          sn = - f / w

          if ( 0 .lt. nu ) then

            do j = 1, n
              x = dreal ( u(j,l1) )
              y = dreal ( u(j,i) )
              u(j,l1) = dcmplx ( x * cs + y * sn, 0.0e+00 )
              u(j,i)  = dcmplx ( y * cs - x * sn, 0.0e+00 )
            end do

          end if

          if ( p .ne. 0 ) then

            do j = n + 1, n + p
              q = a(l1,j)
              r = a(i,j)
              a(l1,j) = q * cs + r * sn
              a(i,j)  = r * cs - q * sn
            end do

          end if

        end do
c
c  Test for convergence.
c
290     continue

        w = s(k)

        if ( l .eq. k ) then
          go to 360
        end if
c
c  Origin shift.
c
        x = s(l)
        y = s(k-1)
        g = t(k-1)
        h = t(k)
        f = ( ( y - w ) * ( y + w ) + ( g - h ) * ( g + h ) ) 
     &    / ( 2.0e+00 * h * y )
        g = dsqrt ( f * f + 1.0e+00 )
        if ( f .lt. 0.0e+00 ) then
          g = -g
        end if
        f = ( ( x - w ) * ( x + w ) + ( y / ( f + g ) - h ) * h ) / x
c
c  QR step.
c
        cs = 1.0e+00
        sn = 1.0e+00
        l1 = l + 1

        do i = l1, k

          g = t(i)
          y = s(i)
          h = sn * g
          g = cs * g
          w = dsqrt ( h * h + f * f )
          t(i-1) = w
          cs = f / w
          sn = h / w
          f = x * cs + g * sn
          g = g * cs - x * sn
          h = y * sn
          y = y * cs

          if ( 0 .lt. nv ) then

            do j = 1, n
              x = dreal ( v(j,i-1) )
              w = dreal ( v(j,i) )
              v(j,i-1) = dcmplx ( x * cs + w * sn, 0.0e+00 )
              v(j,i)   = dcmplx ( w * cs - x * sn, 0.0e+00 )
            end do

          end if

          w = dsqrt ( h * h + f * f )
          s(i-1) = w
          cs = f / w
          sn = h / w
          f = cs * g + sn * y
          x = cs * y - sn * g

          if ( 0 .lt. nu ) then

            do j = 1, n
              y = dreal ( u(j,i-1) )
              w = dreal ( u(j,i) )
              u(j,i-1) = dcmplx ( y * cs + w * sn, 0.0e+00 )
              u(j,i)   = dcmplx ( w * cs - y * sn, 0.0e+00 )
            end do

          end if

          if ( p .ne. 0 ) then

            do j = n + 1, n + p
              q = a(i-1,j)
              r = a(i,j)
              a(i-1,j) = q * cs + r * sn
              a(i,j)   = r * cs - q * sn
            end do

          end if

        end do

        t(l) = 0.0e+00
        t(k) = f
        s(k) = x
        go to 220
c
c  Convergence.
c
360     continue

        if ( w .lt. 0.0e+00 ) then

          s(k) = - w

          if ( 0 .lt. nv ) then

            do j = 1, n
              v(j,k) = - v(j,k)
            end do

          end if
 
        end if

      end do
c
c  Sort the singular values.
c
      do k = 1, n

        g = - 1.0e+00
        j = k

        do i = k, n
          if ( g .lt. s(i) ) then
            g = s(i)
            j = i
          end if
        end do

        if ( j .ne. k ) then

          s(j) = s(k)
          s(k) = g
c
c  Interchange V(1:N,J) and V(1:N,K).
c
          if ( 0 .lt. nv ) then

            do i = 1, n
              q      = v(i,j)
              v(i,j) = v(i,k)
              v(i,k) = q
            end do

          end if
c
c  Interchange U(1:N,J) and U(1:N,K).
c
          if ( 0 .lt. nu ) then
  
            do i = 1, n
              q      = u(i,j)
              u(i,j) = u(i,k)
              u(i,k) = q
            end do

          end if
c
c  Interchange A(J,N1:NP) and A(K,N1:NP).
c
          if ( p .ne. 0 ) then

            do i = n + 1, n + p
              q      = a(j,i)
              a(j,i) = a(k,i)
              a(k,i) = q
            end do

          end if

        end if

      end do
c
c  Back transformation.
c
      if ( 0 .lt. nu ) then

        do kk = 1, n

          k = n + 1 - kk

          if ( b(k) .ne. 0.0e+00 ) then

            q = -a(k,k) / cdabs ( a(k,k) )

            do j = 1, nu
              u(k,j) = q * u(k,j)
            end do

            do j = 1, nu

              q = dcmplx ( 0.0e+00, 0.0e+00 )

              do i = k, m
                q = q + dconjg ( a(i,k) ) * u(i,j)
              end do

              q = q / ( cdabs ( a(k,k) ) * b(k) )

              do i = k, m
                u(i,j) = u(i,j) - q * a(i,k)
              end do

            end do

          end if

        end do

      end if

      if ( 0 .lt. nv ) then

        if ( 1 .lt. n ) then

          do kk = 2, n

            k = n + 1 - kk
            k1 = k + 1

            if ( c(k1) .ne. 0.0e+00 ) then

              q = -dconjg ( a(k,k1) ) / cabs ( a(k,k1) )

              do j = 1, nv
                v(k1,j) = q * v(k1,j)
              end do

              do j = 1, nv

                q = dcmplx ( 0.0e+00, 0.0e+00 )

                do i = k1, n
                  q = q + a(k,i) * v(i,j)
                end do

                q = q / ( cdabs ( a(k,k1) ) * c(k1) )

                do i = k1, n
                  v(i,j) = v(i,j) - q * dconjg ( a(k,i) )
                end do

              end do

            end if

          end do

        end if

      end if

**	write(*,*)
**	do i=1,n
**		write(*,*) i,cmplx(a(i,j))   
**		do j=1,n
xx			q=u(i,j)
xx			u(i,j)=v(i,j)
xx			v(i,j)=q
**			write(*,*) cmplx(u(i,j))-cmplx((v(i,j)))   
**		enddo
**	enddo

      return
      end	   !https://slideplayer.com/slide/5166113/
*-------------------------------------------------------------------------------        
*-------------------------------------------------------------------------------        
      subroutine Dsolve(a,b,n,tol)
*-------------------------------------------------------------------------------        
c
c  esta subrrotina resolve o sistema de equacoes lineares atraves
c  do metodo da eliminacao de gauss, com pivotamento parcial
c  mediante troca de linhas
c
c  a: matriz dos coeficientes
c
c  b: vetor dos termos independentes. apos a resolucao armazena
c     o resultado
c    
c  n: ordem do sistema 
c  tol: tolerancia
c
*-------------------------------------------------------------------------------        
	implicit none
	integer n,i,j,k,l,klin,i_troca
      real*8 b(n),a(n,n),pivo,aux,tol,det 

      i_troca=0
      do i=1,n-1
zz		write(*,*) 'eliminando submatriz',i,' de',n
        pivo=a(i,i)
        klin=i
c
c escolhe na coluna i o elemento de maior modulo(pivo)
c 
c klin indica o indice da linha pivotante
c
        do j=i+1,n
          if (dabs(pivo).lt.dabs(a(j,i))) then
            i_troca=i_troca + 1
            pivo=a(j,i)
            klin=j
          endif
        enddo

        if (dabs(pivo).le.tol) then
          write(*,*)' submatriz contendo linha',i,' nula'
          stop
        endif
c
c  executa o pivotamento seguindo mudanca de linha. armazena
c  os multiplicadores na matriz triangular inferior para o
c  caso de varios vetores de termos independentes
c
c  executa a mudanca de linhas
c
        do    j=1,n
          aux=a(klin,j)
          a(klin,j)=a(i,j)
          a(i,j)=aux
        enddo

        aux=b(i)
        b(i)=b(klin)
        b(klin)=aux
c
c  armazena os multiplicadores
c
        do j=i+1,n
          a(j,i)=a(j,i)/pivo
        enddo
c
c  executa a triangularizacao seguindo direcao das colunas
c
       do k=i+1,n
         b(k)=b(k)-a(k,i)*b(i)
         do l=i+1,n
           a(l,k)=a(l,k)-a(l,i)*a(i,k)
         enddo
       enddo

      enddo
c
c  calculo do determinante da matriz
c
      det=a(1,1)
cc      do k=2,n
cc        det=det*a(k,k)
cc      enddo
cc	write(*,*) ' ...determinante da matriz',det
c
c  chama a subrotina para execucao da retrosubstituicao
c
      call r_retrosub(a,b,n)

zz      write(*,*)' numero de troca de linhas ->',i_troca 
      return
      end
*-------------------------------------------------------------------------------        
      subroutine r_retrosub(a,b,n)
c
c  esta subrrotina resolve o sistema triangular superior resultante seguindo
c  a direcao das colunas
c
c  a: matriz dos coeficientes
c
c  b: vetor dos termos independentes. apos a resolucao armazena
c     o resultado
c    
c  n: ordem do sistema 
c
	implicit none
	integer i,j,n
      real*8 a(n,n),b(n)

      do i=n,1,-1
        b(i)=b(i)/a(i,i)
        do j=1,i-1
          b(j)=b(j)-a(j,i)*b(i)
        enddo
      enddo

      return
      end
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed,during the calculation
!===========================================================
	subroutine inverse(a,c,n)
	implicit none 
	integer n
	double precision a(n,n), c(n,n)
	real*8, allocatable ::  L(:,:), U(:,:), b(:), d(:), x(:)
	double precision coeff
	integer i, j, k
	allocate(L(n,n), U(n,n), b(n), d(n), x(n))


!	step 0: initialization for L, U and b
	do k=1,n
		b(k) = 0.0d0
	    do i=1,n
			L(i,k) = 0.0d0
			U(i,k) = 0.0d0
	    end do
	end do

!	step 1: forward elimination
	do k=1, n-1
	   do i=k+1,n
		coeff=a(i,k)/a(k,k)
		L(i,k) = coeff
		do j=k+1,n
		   a(i,j) = a(i,j)-coeff*a(k,j)
		end do
	   end do
	end do

!	Step 2: prepare L and U matrices 
!	L matrix is a matrix of the elimination coefficient, diagonal elements are 1.0
	do i=1,n
	  L(i,i) = 1.0
	end do

!	U matrix is the upper triangular part of A
	do j=1,n
	  do i=1,j
	    U(i,j) = a(i,j)
	  end do
	end do

!	Step 3: compute columns of the inverse matrix C
	do k=1,n
	  b(k)=1.0
	  d(1) = b(1)
!	Step 3a: Solve Ld=b using the forward substitution
	  do i=2,n
	    d(i)=b(i)
	    do j=1,i-1
		d(i) = d(i) - L(i,j)*d(j)
	    end do
	  end do
!	Step 3b: Solve Ux=d using the back substitution
	  x(n)=d(n)/U(n,n)
	  do i = n-1,1,-1
	    x(i) = d(i)
	    do j=n,i+1,-1
		x(i)=x(i)-U(i,j)*x(j)
	    end do
	    x(i) = x(i)/u(i,i)
	  end do
!	Step 3c: fill the solutions x(n) into column k of C
	  do i=1,n
	    c(i,k) = x(i)
	  end do
	  b(k)=0.0
	end do
	deallocate(L, U, b, d, x)
	end subroutine inverse



	subroutine naleatorio(nale,nran)
	implicit none
	integer n
	integer nale
	integer i,j
	integer isemente
	real*8 nran(nale*110)

	isemente=1
	n=abs(isemente)
	 
	do i=1,nale
		n= MOD(8127*n+28417,134453)
		nran(i)=100*real(n)/134453
	enddo
	end subroutine naleatorio


*	subrotina que introduz condicoes de contorno naturais
	subroutine Load(ng,nnp,ngn,fex,fgl,nvm,ivm)
	implicit none
*
	integer ng,nnp,ngn   
	integer i,j,i1,j1,ivm,nvm  
	real*8 fex(nnp,ngn,nvm),fgl(ng)
*
*	introducao das forcas externas
	do i=1,ng
		fgl(i)=0.0d0
	enddo
	do i=1,nnp
		i1=ngn*(i-1)
		do j=1,ngn
			j1=i1+j
			fgl(j1)=fgl(j1)+fex(i,j,ivm)
		enddo
	enddo
*
	return
	end
	

