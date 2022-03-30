!*****************************************************************************80
!
!! JACOBI_EIGENVALUE carries out the Jacobi eigenvalue iteration.
!
!  Discussion:
!
!    This function computes the eigenvalues and eigenvectors of a
!    real symmetric matrix, using Rutishauser's modfications of the classical
!    Jacobi rotation method with threshold pivoting. 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2013
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(N,N), the matrix, which must be square, real,
!    and symmetric.
!
!    Input, integer IT_MAX, the maximum number of iterations.
!
!    Output, real V(N,N), the matrix of eigenvectors.
!
!    Output, real D(N), the eigenvalues, in descending order.
!
!    Output, integer IT_NUM, the total number of iterations.
!
!    Output, integer ROT_NUM, the total number of rotations.
!
	subroutine jacobi (n,a,v,d,it_max,it_num,rot_num,bw,w,zw)
	implicit none
*
	integer n,i,it_max,it_num,rot_num
	integer j,k,l,m,p,q
*
	real*8 a(n,n),v(n,n),d(n)
	real*8 bw(n),w(n),zw(n)
	real*8 c,g,gapq,h,s,t,tau,term,termp,termq,theta,thresh

*
	do j = 1, n
		do i = 1, n
			v(i,j) = 0.0D+00
		end do
		v(j,j) = 1.0D+00
	end do
*
	do i = 1, n
		d(i) = a(i,i)
		bw(i) = d(i)
		zw(i) = 0.0D+00
	end do
*
	it_num = 0
	rot_num = 0
	do while ( it_num < it_max )
		it_num = it_num + 1
*
*  The convergence threshold is based on the size of the elements in
*  the strict upper triangle of the matrix.
*
		thresh = 0.0D+00
		do j = 1, n
			do i = 1, j - 1
				thresh = thresh + a(i,j) ** 2
			end do
		end do
		thresh = dsqrt ( thresh ) / dfloat ( 4 * n )
	    if ( thresh .eq. 0.0D+00 ) exit 
*
		do p = 1, n
			do q = p + 1, n
				gapq = 10.0D+00 * dabs ( a(p,q) )
				termp = gapq + dabs ( d(p) )
				termq = gapq + dabs ( d(q) )
*
*				Annihilate tiny offdiagonal elements.
*
zz				if ( 4 .lt. it_num .and. 
zz	1		        termp.eq.dabs(d(p)) .and. 
zz	1		        termq.eq.dabs(d(q)) ) then
zz					a(p,q) = 0.0D+00
*
*				Otherwise, apply a rotation.
*
zz				else if ( thresh .le. dabs ( a(p,q) ) ) then
					h = d(q) - d(p)
					term = dabs ( h ) + gapq
					if ( term .eq. abs ( h ) ) then
						t = a(p,q) / h
					else
						theta = 0.5D+00 * h / a(p,q)
						t = 1.0D+00 / ( dabs ( theta ) 
	1					+ sqrt ( 1.0D+00 + theta * theta ) )
						if ( theta < 0.0D+00 ) then 
							t = - t
						end if
					end if
					c = 1.0D+00 / dsqrt ( 1.0D+00 + t * t )
					s = t * c
					tau = s / ( 1.0D+00 + c )
					h = t * a(p,q)
*
*					Accumulate corrections to diagonal elements.
*
					zw(p) = zw(p) - h                  
					zw(q) = zw(q) + h
					d(p) = d(p) - h
					d(q) = d(q) + h
					a(p,q) = 0.0D+00
*
*					Rotate, using information from the upper triangle of A only.
*
					do j = 1, p - 1
						g = a(j,p)
						h = a(j,q)
						a(j,p) = g - s * ( h + g * tau )
						a(j,q) = h + s * ( g - h * tau )
					end do
					do j = p + 1, q - 1
						g = a(p,j)
						h = a(j,q)
						a(p,j) = g - s * ( h + g * tau )
						a(j,q) = h + s * ( g - h * tau )
					end do
					do j = q + 1, n
						g = a(p,j)
						h = a(q,j)
						a(p,j) = g - s * ( h + g * tau )
						a(q,j) = h + s * ( g - h * tau )
					end do
*
*					Accumulate information in the eigenvector matrix.
*
					do j = 1, n
						g = v(j,p)
						h = v(j,q)
						v(j,p) = g - s * ( h + g * tau )
						v(j,q) = h + s * ( g - h * tau )
					end do
					rot_num = rot_num + 1
zz				end if
			end do
		end do
		do i = 1, n
			bw(i) = bw(i)+zw(i)
			d(i) = bw(i)
			zw(i) = 0.0D+00
		end do
	end do
*
*	Restore upper triangle of input matrix.
*
	do j = 1, n
		do i = 1, j - 1
			a(i,j) = a(j,i)
		end do
	end do
*
*	Ascending sort the eigenvalues and eigenvectors.
*
	do k = 1, n - 1
		m = k
		do l = k + 1, n
			if ( d(l) .lt. d(m) ) m = l
		end do
		if ( m.ne.k ) then
			t    = d(m)
			d(m) = d(k)
			d(k) = t
 			do i = 1, n
				w(i) = v(i,m)
				v(i,m) = v(i,k)
				v(i,k) = w(i)
			end do
		end if
	end do
	return
	end
