      subroutine subsp(a,b,v,t,g,h,d,dp,dtol,p,nf,nv,neq,imas
     &                 ,shift,tol,prt,its)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2017: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Subspace iteration to extract lowest nf eigenpairs
!               Solves:  (A - shift*B)*V = B*V*d for V and d

!      Inputs:
!         a(*)      - Coefficient matrix (tangent stiffness)
!         b(*)      - Coefficient matrix (mass or geometric stiffness)
!         nf        - Number of pairs to converge
!         nv        - Size of subspace problem > or = nf
!         neq       - Size of A,B
!         imas      - Switch: =1 for consistent B; =2 for diagonal B
!         shift     - Value of shift
!         tol       - Tolerance to converge pairs
!         prt       - Flag, output iteration arrays if true
!         its       - Maximum number of subspace iterations

!      Scratch:
!         t(neq)    - Working vector
!         g(*)      - Projection of A - shift*B matrix
!         h(*)      - Projection of B           matrix
!         dp(*)     - Previous iteration values
!         dtol(*)   - Tolerance of eigenvalue iterations
!         p(nv,*)   - Eigenvectors of G*P = H*P*d

!      Outputs:
!         v(neq,*)  - Eigenvectors
!         d(*)      - Eigenvalues
!-----[--.----+----.----+----.-----------------------------------------]

      implicit  none




      logical   conv,prt,soltyp
      integer   i,j,k,n
      integer   nf,nv,neq,imas,its, itlim, it,itt, num,nmas
      real*4    etime, tary(2)
      real*8    shift,tol, dm,dr

      integer   pdf(1)
      real*8    a(*),b(*),v(neq,*),t(neq)
      real*8    g(*),h(*),d(*),dp(*),dtol(*),p(nv,*),dpp(4)
	real*8 gtan(1) !zzz
	
      save

      data      pdf / 1 /

!     Compute initial iteration vectors

      call pzero(v,nv*neq)
zzz      num    = 0
      nmas   = 0
zzz      soltyp = ittyp.eq.-1
      do n = 1,neq

!       Count number of values less than current shift

zzz        if(soltyp .and. a(n).lt.0.0d0) num = num + 1

!       Count number of nonzero masses

        if(b(n).ne.0.0d0) nmas = nmas + 1

      end do

zzz      if(soltyp) then
zzz        write(iow,2002) num,shift
zzz        if(ior.lt.0) write(*,2002) num,shift
zzz      endif
      nmas = nmas/nv
      i = 0
      j = 1
      do n = 1,neq
        dm = b(n)
        if(dm.ne.0.0d0) then
          v(n,j) = dm
          i      = i + 1
          if(mod(i,nmas).eq.0) j = min(j+1,nv)
        endif
      end do

      do i = 1,nv
        dp(i)   = 0.0
        dtol(i) = 1.0d0
        call scalev(v(1,i),pdf,1,1,neq)
      end do

!     Compute new vectors and project 'a' onto 'g'

      conv = .false.
      itlim = its
      if(nv.eq.nf) itlim = 1
      do it = 1,itlim
        itt = it

!       Project 'b' matrix to form 'h' and compute 'z' vectors

        call sprojb(b,v,t,h,neq,nv,imas)

!       Project 'a' matrix to form 'g'

        call sproja(v,t,g,neq,nv)

!       Solve reduced eigenproblem

zzz        if(imtyp.eq.1) then

!         Eigenproblem: 'g*p = h*p*d'; e.g., vibration modes

          call geig(g,h,d,p,t,nv,1,prt)

zzz        endif

!       Compute new iteration vector 'u' from 'z'

        do i = 1,neq

!         Move row of 'v' into temporary vector 't'

          do j = 1,nv
            t(j) = v(i,j)
          end do

!         Compute new iiteration vector entries

          do j = 1,nv
            v(i,j) = 0.0d0
            do k = 1,nv
              v(i,j) = v(i,j) + t(k)*p(k,j)
            end do
          end do
        end do

!       Check for convergence

        do n = 1,nv
          if(d(n).ne.0.0d0) dtol(n) = abs((d(n)-dp(n))/d(n))
          dp(n) = d(n)
        end do

!       Compute maximum error for iteration

        conv = .true.
        dm   =  0.0d0
        do n = 1,nf
          if(dtol(n).gt.tol) conv = .false.
          dm = max(dm,dtol(n))
        end do

zzz        if(prt) then
zzz          if(ior.gt.0) write(iow,2006) it,(d(n),n=1,nv)
zzz          if(ior.lt.0) write(  *,2006) it,(d(n),n=1,nv)
zzz          if(itlim.gt.1) then
zzz            if(ior.gt.0) write(iow,2001) it,(dtol(n),n=1,nv)
zzz            if(ior.lt.0) write(  *,2001) it,(dtol(n),n=1,nv)
zzz          endif
zzz        endif

        if(conv) go to 200
zzz        if(ior.lt.0) write(*,2005) it, etime(tary), dm
      end do

!     Scaled vectors mass orthonormal

200   gtan(1)=1.0d0  !zzz
	dm = 1.d0/gtan(1)
      do n = 1,nv
        d(n) = (1.0/d(n) + shift)*dm
        dp(n)= sqrt(abs(d(n)))
      end do

!     Output solution values

zzz      write(iow,2000) itt,(d(n),n=1,nv)
zzz      if(itt.gt.1) write(iow,2001) itt,(dtol(n),n=1,nv)
zzz      if(imtyp.eq.1) then
zzz        write(iow,2003) (dp(n),n=1,nv)
zzz        dm = 0.5d0/acos(-1.0d0)
zzz        write(iow,2004) (dp(n)*dm,n=1,nv)
zzz        write(iow,2007)
zzz        dr = 1.d0/dm
zzz        do i = 0,nv-1,4
zzz          do j = 1,min(4,nv-i)
zzz            if(abs(dp(i+j)).gt.tol*10.d0) then
zzz              dpp(j) = dr/dp(i+j)
zzz            else
zzz              dpp(j) = 0.0d0
zzz            endif
zzz          end do
zzz          write(iow,2008) (dpp(j),j=1,min(4,nv-i))
zzz        end do
zzz      endif
zzz      if(ior.lt.0) then
zzz        write(*,2000) itt,(d(n),n=1,nv)
zzz        if(itt.gt.1) write(*,2001) itt,(dtol(n),n=1,nv)
zzz        if(imtyp.eq.1) then
zzz          write(*,2003) (dp(n),n=1,nv)
zzz          write(*,2004) (dp(n)*dm,n=1,nv)
zzz          write(*,2007)
zzz          do i = 0,nv-1,4
zzz            do j = 1,min(4,nv-i)
zzz              if(abs(dp(i+j)).gt.tol*10.d0) then
zzz                dpp(j) = dr/dp(i+j)
zzz              else
zzz                dpp(j) = 0.0d0
zzz              endif
zzz            end do
zzz            write(*,2008) (dpp(j),j=1,min(4,nv-i))
zzz          end do
zzz        endif
zzz      endif

!     Formats

zzz2000  format(/'  SUBSPACE: Current eigenvalues, iteration',i4/
zzz     +        (5x,1p,4d17.8))

zzz2001  format( '  SUBSPACE: Current residuals,   iteration',i4/
zzz     +        (5x,1p,4d17.8))

zzz2002  format(5x,'There are',i4,' eigenvalues less than shift',1p,e12.4)

zzz2003  format(/'  SUBSPACE: Square root of eigenvalues (rad/t)'/
zzz     +        (5x,1p,4d17.8))

zzz2004  format(/'  SUBSPACE: Square root of eigenvalues (Hz.)'/
zzz     +        (5x,1p,4d17.8))

zzz2005  format('  Completed subspace iteration',i4,' Time =',f9.2,
zzz     +       '  Max. tol =',1p,e12.5)

zzz2006  format(/'  SUBSPACE: Inverse eigenvalues, iteration',i4/
zzz     +        (5x,1p,4d17.8))

zzz2007  format(/'  SUBSPACE: Period in units of time      (T)')

zzz2008  format(5x,1p,4d17.8)

      end
***********************************************************************
!$Id:$
      subroutine pzero(v,nn)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2017: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Zero real array of data

!      Inputs:
!         nn     - Length of array

!      Outputs:
!         v(*)   - Array with zero values
!-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      integer   n,nn
      real*8    v(nn)

      save

      do n = 1,nn
        v(n) = 0.0d0
      end do

      end
***********************************************************************
!$Id:$
      subroutine scalev(v,pdf,ndm,ndf,numnp)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2017: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Scale vector to have maximum element of +1.0

!      Inputs:
!         v(ndf,*) - Vector of values
!         pdf(*)   - DOF to scale on
!         ndm      - Space dimension of mesh
!         ndf      - DOF's/node (maximum)
!         numnp    - Number of nodes

!      Outputs:
!         v(ndf,*) - Unit vector of values
!-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      integer   i,n,ndm,ndf,numnp
      integer   pdf(*)
      real*8    v(ndf,*),vmax

      save

!     Locate maximum

      vmax = 0.0d0
      do i = 1,ndm
        if(pdf(i).ge.1.and.pdf(i).le.ndf) then
          do n = 1,numnp
            vmax = max(vmax,abs(v(pdf(i),n)))
          end do ! n
        endif
      end do ! i

!     Perform scaling

      if(vmax.gt.0.0d0) then
        vmax = 1.d0/vmax
        do n = 1,numnp
          do i = 1,ndf
            v(i,n) = v(i,n)*vmax
          end do ! i
        end do ! n
      else
        write(*,*) 'Zero length vector in SCALEV'
      endif

      end
***********************************************************************
!$Id:$
      subroutine sprojb(b,v,t,h,neq,nv,imas)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2017: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Compute subspace projection of 'b' to form 'h'

!      Inputs:
!         b(*)     - Symmetric coefficient matrix for eigenproblem
!         v(neq,*) - Set of iteration vectors
!         neq      - Number of equations in B
!         nv       - Size of projected matrix
!         imas     - Mass type: 1 = consistent; 2 = diagonal.

!      Scratch:
!         t(neq)   - Working vector

!      Outputs:
!         h(*)     - Projected matrix V_trans * B * V
!-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   neq,nv,imas, i,j,k
      real*8    b(*),v(neq,*),t(*),h(*), dot

      save

!     Compute 'z' and 'b' projection to form 'h'

      do j = 1,nv

!       Consistent mass

        if(imas.eq.1) then
          call pzero(t,neq)
		call multmv(b,v(1,j),t,neq)

!       Lumped mass

        else
          do i = 1,neq
            t(i) = v(i,j)*b(i)
          end do
        endif

!       Project 'z' and 'v' vectors to form 'h'

        k = j*(j+1)/2
        do i = j,nv
          h(k) = dot(t,v(1,i),neq)
          k = k + i
        end do
        do i = 1,neq
          v(i,j) = t(i)
        end do
      end do

      end
***********************************************************************
!$Id:$
      subroutine sproja(v,t,g,neq,nv)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2017: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Compute subspace projection of 'a' to form 'g'

!      Inputs:
!         v(neq,*) - Set of iteration vectors
!         neq      - Number of equations in A
!         nv       - Size of projected matrix

!      Scratch:
!         t(neq)   - Working vector

!      Outputs:
!         g(*)     - Projected matrix V_trans * A * V
!-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   neq,nv, i,j,k
      real*8    en

      real*8    v(neq,*),t(neq),g(*)

      real*8    dot

      save

!     Forward reduce eigenvector estimates

      k = 0
      do j = 1,nv

!     Copy vector 'v' into 'z' and solve equations

        call pmove (v(1,j),t(1),neq)
        call dasol (hr(np(1)+neq),hr(np(1)+neq),hr(np(1)),
     &              v(1,j),mr(np(21)),neq,neq,en)

!     Compute projection of stiffness

        do i = 1,j
          k = k + 1
          g(k) = dot(v(1,i),t(1),neq)
        end do
      end do

      end
***********************************************************************
!$Id:$
      subroutine pmove(a,b,nn)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2017: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Move real array a into b

!      Inputs:
!         a(*)      - Array to move
!         nn        - Length of array to move

!      Outputs:
!         b(*)      - Moved array
!-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      integer   n,nn
      real*8    a(nn),b(nn)

      save

!     Move a-array into b-array

      do n = 1,nn
        b(n) = a(n)
      end do

      end
***********************************************************************
!$Id:$
      subroutine dasol(al,au,ad,b,jp,neqs, neqt, energy)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2017: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Solution of algebraic equations stored in profile form

!         Equations '   1  ' to 'neqs' are symmetric.
!         Equations 'neqs+1' to 'neqt' are unsymmetric.

!         Use:
!          a.) All equations are unsymmetric       : neqs = 1 (or 0)
!              N.B.  The top 1 x 1 submatrix is always symmetric.
!                    Both 'al' and 'au' must be provided.

!          b.) All equations are symmetric         : neqs = neqt
!              N.B.  In this case the array 'al' is not used.

!          c.) First 'neqs' equations are symmetric: 1 < neqs < neqt
!              N.B.  Storage of 'al' for unsymmetric equations only.

!      Coefficient matrix must be decomposed into its triangular
!      factors using 'DATRI' before using 'DASOL'.

!      Inputs:
!         al(*)  - Lower triangular factors of A
!         au(*)  - Upper triangular factors of A
!         ad(*)  - Diagonal factors of A
!         b(*)   - Right hand side vector
!         jp(*)  - Pointer array for row/columns of 'al', 'au'.
!         neqs   - Number of symmetric equations
!         neqt   - Number of equations to solve.

!      Outputs:
!         b(*)   - Solution vector, x.
!         energy - Energy of solution: x*A*x
!-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   jp(*)
      integer   is, j, jr, jh, neqs, neqt, neq
      real*8    bd, energy, al(*),au(*),ad(*),b(*), dot

      save

!     Initialize energy

      energy = 0.0d0

!     Find first non-zero entry in right hand side

      do is = 1,neqt
        if(b(is).ne.0.0d0) go to 100
      end do
      if(ior.gt.0) write(iow,2000)
      if(ior.lt.0) write(*,2000)
      return

!     Reduce right hand side

!     Do symmetric part

100   neq = max(1, neqs)
      do j = is+1,neq
        jr = jp(j-1)
        jh = jp(j) - jr
        if(jh.gt.0) then
          b(j) = b(j) - dot(au(jr+1),b(j-jh),jh)
        endif
      end do

!     Do unsymmetric part

      do j = max(is,neq)+1,neqt
        jr = jp(j-1)
        jh = jp(j) - jr
        if(jh.gt.0) then
          jr   = jr   - jp(neq)
          b(j) = b(j) - dot(al(jr+1),b(j-jh),jh)
        endif
      end do

!     Multiply by inverse of diagonal elements

      do j = is,neqt
        bd = b(j)
        b(j) = b(j)*ad(j)
        energy = energy + bd*b(j)
      end do

!     Symmetric and unsymmetric backsubstitution

      do j = neqt,2,-1
        jr = jp(j-1)
        jh = jp(j) - jr
        if(jh.gt.0) then
          call colred(au(jr+1),b(j),jh, b(j-jh))
        endif
      end do

!     Warning format

2000  format(' *WARNING* Zero right-hand-side vector')

      end
***********************************************************************
!$Id:$
      subroutine colred(au,xj,nn, b)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2017: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Columnwise reduction for back substitution

!      Inputs:
!         au(*)   - Upper column of reduced array A
!         xj      - Solution of reduced column
!         nn      - Length to reduce
!         b(*)    - Unreduced column

!      Outputs:
!         b(*)    - Reduced column
!-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   n,nn
      real*8    xj, au(*),b(*)

      do n = 1,nn
        b(n) = b(n) - au(n)*xj
      end do

      end
***********************************************************************
!$Id:$
      subroutine geig(g,h,d,p,t,nv,nvs,prt)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2017: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Solve general eigenproblem 'g*p = h*p*d'

!      Inputs:
!         g(*)   - Left hand projected array
!         h(*)   - Right hand projected array
!         nv     - Size of problem
!         nvs    - Eigenvectors h orthogonal if positive
!         prt    - Output computations if true

!      Outputs:
!         d(*)   - Eigenvalues of problem
!         p(*,*) - Eigenvectors of problem
!-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      logical   prt
      integer   nv,nvs,ir, i, j

      real*8    g(*),h(*),d(*),p(nv,*),t(*)

      save

!     Output projected matrices

      if(prt) then
        call wprojm(g,nv,'Projected G')
        call wprojm(h,nv,'Projected H')
      endif

!     Compute standard eigenvalue problem matrix 'c'

      call chlfwd(h,g,p,nv)

!     Perform eignfunction decomposition of 'c'

      call eisql(g,d,t,p,nv,ir)

!     Compute vectors of original problem

      call chlbac(h,p,nv)

      if(prt) call mprint(p,nv,nv,nv,'vectors  p')

!     Divide eigenvectors by eigenvalue to prevent overflows

      ir = 0
      do j = 1,nv
        ir = ir + j
        if(nvs.gt.0 .and. d(j).ne.0.0d0) then
          t(1) = 1.0d0/d(j)
        elseif(nvs.le.0 .and. h(ir).ne.0.0d0) then
          if(d(j).ne.0.0d0) then
            t(1) = abs(d(j))/sqrt(abs(h(ir)))
          else
            t(1) = 1.0d0/sqrt(abs(h(ir)))
          endif
        else
          t(1) = 1.0d0
        endif
        if(p(j,j).lt.-0.00001d0) t(1) = -t(1)
        do i = 1,nv
          p(i,j) = p(i,j)*t(1)
        end do
      end do

      end
***********************************************************************
!$Id:$
      subroutine wprojm(a,nn,ah)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2017: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Outputs projected subspace arrays: G and H

!      Inputs:
!         a(*)        - Array to output
!         nn          - Number row/columns in array
!         ah          - Header to write (G or H)

!      Outputs:
!         none
!-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      character ah*(*)
      real*8    a(*)
      integer   nn,i,k,n

      save

      i = 0
      write(iow,2000) ah
      do n = 1,nn
        write(iow,2001) (a(i+k),k=1,n)
        i = i + n
      end do
      if(ior.lt.0) then
        i = 0
        write(  *,2000) ah
        do n = 1,nn
          write(  *,2001) (a(i+k),k=1,n)
          i = i + n
        end do
      endif

!     Formats

2000  format(' ',a,'-Matrix ')

2001  format(1p,8d10.2)

      end
***********************************************************************
!$Id:$
      subroutine chlfwd(u,g,s,nn)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2017: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Use Cholesky factors to project onto a standard eigenproblem

!      Inputs:
!         g(*)  - Symmetric projected matrix
!         u(*)  - Upper factor for projection
!         nn    - Size of arrays

!      Outputs:
!         s(*,*) - Projected array

!-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,j,id,im,jd,nn
      real*8    u(*),g(*),s(nn,nn), dot

!     Choleski factorization of a symmetric, positive definite matrix

      u(1) = 1.d0/sqrt(abs(u(1)))
      jd   = 1
      do j = 2,nn
        id = 0
        do i = 1,j-1
          if(i.gt.1) u(jd+i) = u(jd+i) - dot(u(id+1),u(jd+1),i-1)
          id = id + i
          u(jd+i) = u(jd+i)*u(id)
        end do
        u(jd+j) = 1.d0/sqrt(abs(u(jd+j) - dot(u(jd+1),u(jd+1),j-1)))
        jd = jd + j
      end do

!     Perform forward solutions to get projected matrix

      s(1,1) = g(1)*u(1)
      id = 1
      do i = 2,nn
        s(1,i) = g(id+1)*u(1)
        im = i - 1
        jd = 0
        do j = 1,im
         s(i,j) = (g(id+j) - dot(u(id+1),s(1,j),im))*u(id+i)
         if(j.gt.1) s(j,i) = (g(id+j)-dot(u(jd+1),s(1,i),j-1))*u(jd+j)
         jd = jd + j
        end do
        id = id + i
        s(i,i) = (g(id) - dot(u(id-im),s(1,i),im))*u(id)
      end do

!     Complete projection

      g(1) = s(1,1)*u(1)
      jd = 2
      do j = 2,nn
        g(jd) = s(j,1)*u(1)
        id = 2
        do i = 2,j
          im = i - 1
          g(jd+im) = (s(j,i) - dot(u(id),g(jd),im))*u(id+im)
          id = id + i
        end do
        jd = jd + j
      end do

      end
***********************************************************************
!$Id:$
      subroutine eisql(a,d,e,z,n,ierr)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2017: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Compute eigen-pairs for standard eigenproblem.

!      Inputs:
!         a(*)   - Matrix for wanted eigenvalues
!         n      - size of eigenproblem

!      Outputs:
!         d(n)   - Eigenvalues
!         z(n,n) - Eigenvectors
!         ierr   - Error indicator

!      Scratch:
!         e(*)   - Working vector
!-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,j,k,l,m,n,ierr,jp1,ii,l1,mml,n2
      real*8    b,c,f,g,h,hh,p,r,s,scale

      real*8    a(*),d(*),e(*),z(n,n)

!     Eispac QL algorithm adapted from 'tred2' and 'tql2'

      n2 = 0
      do i = 1,n
        do j = 1,i
          n2 = n2 + 1
          z(i,j) = a(n2)
        end do
      end do
      if(n.eq.1) go to 300
      n2 = n + 2
      do ii = 2,n
      i = n2 - ii
      l = i - 1
      h = 0.0d0
      scale = 0.0d0
      if(l.lt.2) then
        e(i) = z(i,l)
      else
        do k = 1,l
          scale = scale + abs(z(i,k))
        end do
        if(scale.ne.0.0d0) go to 100
        e(i) = z(i,l)
      endif
      go to 200
100   do k = 1,l
        z(i,k) = z(i,k)/scale
        h = h + z(i,k)*z(i,k)
      end do
      f = z(i,l)
      g = -sign(dsqrt(h),f)
      e(i) = scale*g
      h = h - f*g
      z(i,l) = f - g
      f = 0.0d0
      do j = 1,l
      z(j,i) = z(i,j)/h
      g = 0.0d0
      do k = 1,j
        g = g + z(j,k)*z(i,k)
      end do
      jp1 = j + 1
      if(l.ge.jp1) then
        do k = jp1,l
          g = g + z(k,j)*z(i,k)
        end do
      endif
      e(j) = g/h
      f = f + e(j)*z(i,j)
      end do
      hh = f/(h+h)
      do j = 1,l
        f = z(i,j)
        g = e(j) - hh*f
        e(j) = g
        do k = 1,j
          z(j,k) = z(j,k) - f*e(k) - g*z(i,k)
        end do
      end do
200   d(i) = h
      end do

!     Set transformation array for ql

300   d(1) = z(1,1)
      z(1,1) = 1.0d0
      e(1) = 0.0d0
      ierr = 0
      if(n.eq.1) return
      do i = 2,n
        l = i - 1
        if(d(i).ne.0.0d0) then
          do j = 1,l
            g = 0.0d0
            do k = 1,l
              g = g + z(i,k)*z(k,j)
            end do
            do k = 1,l
              z(k,j) = z(k,j) - g*z(k,i)
            end do
          end do
        endif
        d(i) = z(i,i)
        z(i,i) = 1.0d0
        do j = 1,l
        z(i,j) = 0.0d0
        z(j,i) = 0.0d0
        end do
      end do

!     Begin 'QL' algorithm on tridagonal matrix now stored in 'd' and 'e

      do i = 2,n
        e(i-1) = e(i)
      end do
      f = 0.0d0
      b = 0.0d0
      e(n) = 0.0d0
      do l = 1,n
        j = 0
        h = epmac*(abs(d(l)) + abs(e(l)))
        if(b.lt.h) b = h
        do m = l,n
          if(abs(e(m)).le.b) go to 400
        end do
400     if(m.ne.l) then
410       if(j.eq.30) go to 500
          j = j + 1
          l1 = l + 1
          g = d(l)
          p = (d(l1)-g)/(e(l)+e(l))
          r = dsqrt(p*p+1.0d0)
          d(l) = e(l)/(p+sign(r,p))
          h = g - d(l)
          do i = l1,n
            d(i) = d(i) - h
          end do
          f = f + h
          p = d(m)
          c = 1.0d0
          s = 0.0d0
          mml = m - l
          do ii = 1,mml
          i = m - ii
          g = c*e(i)
          h = c*p
          if(abs(p).ge.abs(e(i))) then
            c = e(i)/p
            r = dsqrt(c*c+1.0d0)
            e(i+1) = s*p*r
            s = c/r
            c = 1.0d0/r
          else
            c = p/e(i)
            r = dsqrt(c*c+1.0d0)
            e(i+1) = s*e(i)*r
            s = 1.0d0/r
            c = c*s
          endif
          p = c*d(i) - s*g
          d(i+1) = h + s*(c*g + s*d(i))
          do k = 1,n
            h = z(k,i+1)
            z(k,i+1) = s*z(k,i) + c*h
            z(k,i  ) = c*z(k,i) - s*h
          end do
        end do
        e(l) = s*p
        d(l) = c*p
        if(abs(e(l)).gt.b) go to 410
      endif
      d(l) = d(l) + f
      end do
      do ii = 2,n
        i = ii - 1
        k = i
        p = d(i)
        do j = ii,n
          if(abs(d(j)).gt.abs(p)) then
            k = j
            p = d(j)
          endif
        end do
        if(k.ne.i) then
          d(k) = d(i)
          d(i) = p
          do j = 1,n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
          end do
        end if
      end do

      return

500   ierr = l

      end
***********************************************************************
!$Id:$
      subroutine chlbac(u,s,nn)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2017: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Back substitution for Cholesky factors in eigen
!               solutions.

!      Inputs:
!        u(*)   - Unreduced array
!        s(*,*) - Factored array of matrix
!        nn     - Size of arrays

!      Outputs:
!        u(*)   - Solution after back substitution

!-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,j,jd, nn
      real*8    u(*),s(nn,nn)

      save

!     Compute eigenvalues of general linear problem by backsubstitution

      j  = nn
      jd = nn*(nn+1)/2
      do i = 1,nn
        s(nn,i) = s(nn,i)*u(jd)
      end do

      do j = nn,2,-1
        jd = jd - j
        do i = 1,nn
          call colbac(u(jd+1),s(1,i),u(jd),j-1)
        end do
      end do

      end
***********************************************************************
!$Id:$
      subroutine colbac(u,s,d,jj)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2017: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Backsubstitution macro for eigen solution

!      Inputs:
!         s(*)  - Unreduced column
!         u(*)  - Column of upper array already reduced
!         d     - Solution value for 'u' column
!         jj    - Length to reduce

!      Outputs:
!         s(*)  - Reduced column
!-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   j,jj
      real*8    d,u(*),s(*)

      do j = 1,jj
        s(j) = s(j) - u(j)*s(jj+1)
      end do
      s(jj) = s(jj)*d

      end
***********************************************************************
!$Id:$
      subroutine mprint(a,ii,jj,mm,name)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2017: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Output array of integer values

!      Inputs:
!         a(mm,*)  - Array to output
!         ii       - Number of rows to output
!         jj       - Number of columns to output
!         mm       - Dimension of array
!         name     - Name of array to appear with outputs

!      Outputs:
!         none
!-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      character name*(*), aname*30
      integer   ii,jj,mm,nn, ja,jb,j,i,n,irow
      real*8    a(mm,*)

!     Print an ii x jj array whose dimension is mm for first subscript

      aname = name
      nn    = (jj+5)/6
      jb    = 0

      do n = 1,nn

        ja = jb + 1
        jb = min(jj,ja + 5)

!       Output to O file

        if(n.eq.1 .or. ii.gt.1) then
          write(iow,2000) aname,(j,j=ja,jb)
        endif
        do i = 1,ii
          if(ii.eq.1) then
            irow = n
          else
            irow = i
          endif
          write(iow,2001) irow,(a(i,j),j=ja,jb)
        end do

!       Output to screen if interactive

        if(ior.lt.0) then
          if(n.eq.1 .or. ii.gt.1) then
            write(*,2000) aname,(j,j=ja,jb)
          endif
          do i = 1,ii
            if(ii.eq.1) then
              irow = n
            else
              irow = i
            endif
            write(*,2001) irow,(a(i,j),j=ja,jb)
          end do
        endif

      end do

!     Formats

 2000 format(/4x,'Matrix: ',a30/4x,'row/col',i6,5i12)
 2001 format(i8,1p,6e12.4)

      end
***********************************************************************
***********************************************************************
***********************************************************************
***********************************************************************
***********************************************************************
***********************************************************************
***********************************************************************
