C
C Compute all of the eigenvalues and eigenvectors of a real matrix.
C
C In this example, the eigenvalues and eigenvectors of a real matrix
C are computed and printed.  The performance index is also computed
C and printed.  This serves as a check on the computations.  For more
C details, see IMSL routine EPIRG.
C
C Output:
C
C                        EVAL
C                1                2                3
C  ( 2.000, 4.000)  ( 2.000,-4.000)  ( 1.000, 0.000)
C
C                             EVEC
C                     1                  2                  3
C  1  ( 0.3162, 0.3162)  ( 0.3162,-0.3162)  ( 0.4082, 0.0000)
C  2  ( 0.0000, 0.6325)  ( 0.0000,-0.6325)  ( 0.8165, 0.0000)
C  3  ( 0.6325, 0.0000)  ( 0.6325, 0.0000)  ( 0.4082, 0.0000)
C
C  Performance index =  0.037
C
      USE MSIMSLMS
C                                  Declare variables
      INTEGER    LDA, LDEVEC, N
      PARAMETER  (N=3, LDA=N, LDEVEC=N)
 
      INTEGER    NOUT
      REAL       PI
      COMPLEX    EVAL(N), EVEC(LDEVEC,N)
 
      REAL       A(LDA,N)

C                                  Define values of A:
C
C                                  A = (  8.0   -1.0   -5.0  )
C                                      ( -4.0    4.0   -2.0  )
C                                      ( 18.0   -5.0   -7.0  )
C
      DATA A/8.0, -4.0, 18.0, -1.0, 4.0, -5.0, -5.0, -2.0, -7.0/
C
C                                  Find eigenvalues and vectors of A
      CALL EVCRG (N, A, LDA, EVAL, EVEC, LDEVEC)
C                                  Compute performance index
      PI = EPIRG(N,N,A,LDA,EVAL,EVEC,LDEVEC)
C                                  Print results
      CALL UMACH (2, NOUT)
      CALL WRCRN ('EVAL', 1, N, EVAL, 1, 0)
      CALL WRCRN ('EVEC', N, N, EVEC, LDEVEC, 0)
      WRITE (NOUT,'(/,A,F6.3)') ' Performance index = ', PI
      END
