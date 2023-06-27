PROGRAM zheevr_example
    implicit none
!     .. Parameters ..
      INTEGER          N
      PARAMETER        ( N = 4 )
      INTEGER          LDA, LDZ
      PARAMETER        ( LDA = N, LDZ = N )
      INTEGER          LWMAX
      PARAMETER        ( LWMAX = 1000 )
!
!     .. Local Scalars ..
      INTEGER          INFO, LWORK, LRWORK, LIWORK, IL, IU, M
      DOUBLE PRECISION ABSTOL, VL, VU
!
!     .. Local Arrays ..
      INTEGER          ISUPPZ( N ), IWORK( LWMAX )
      DOUBLE PRECISION W( N ), RWORK( LWMAX )
      COMPLEX*16       A( LDA, N ), Z( LDZ, N ), WORK( LWMAX )
      DATA             A/&
       (-2.16, 0.00),(-0.16, 4.86),(-7.23, 9.38),(-0.04,-6.86),&
        ( 0.00, 0.00),( 7.45, 0.00),( 4.39,-6.29),(-8.11, 4.41),&
         ( 0.00, 0.00),( 0.00, 0.00),(-9.03, 0.00),(-6.89, 7.66),&
          ( 0.00, 0.00),( 0.00, 0.00),( 0.00, 0.00),( 7.76, 0.00)&
                            /
!
!     .. External Subroutines ..
      EXTERNAL         ZHEEVR
      EXTERNAL         PRINT_MATRIX, PRINT_RMATRIX
!
!     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
!
!     .. Executable Statements ..
      WRITE(*,*)'ZHEEVR Example Program Results'
!     Negative ABSTOL means using the default value
      ABSTOL = -1.0
!     Set VL, VU to compute eigenvalues in half-open (VL,VU] interval
      VL = -5.0
      VU = 5.0
!
!     Query the optimal workspace.
!
      LWORK = -1
      LRWORK = -1
      LIWORK = -1
      CALL ZHEEVR( 'Vectors', 'Values', 'Lower', N, A, LDA, VL, VU, IL,&
                   IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, RWORK,&
                                LRWORK, IWORK, LIWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      LRWORK = MIN( LWMAX, INT( RWORK( 1 ) ) )
      LIWORK = MIN( LWMAX, IWORK( 1 ) )
!
!     Solve eigenproblem.
!
      CALL ZHEEVR( 'Vectors', 'Values', 'Lower', N, A, LDA, VL, VU, IL,&
                   IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, RWORK,&
                                LRWORK, IWORK, LIWORK, INFO )
!
!     Check for convergence.
!
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      END IF
!
!     Print the number of eigenvalues found.
!
      WRITE(*,'(/A,I2)')' The total number of eigenvalues found:', M
!
!     Print eigenvalues.
!
      CALL PRINT_RMATRIX( 'Selected eigenvalues', 1, M, W, 1 )
!
!     Print eigenvectors.
!
      CALL PRINT_MATRIX( 'Selected eigenvectors (stored columnwise)',&
                         N, M, Z, LDZ )
      STOP
END PROGRAM




SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      COMPLEX*16       A( LDA, * )
!
      INTEGER          I, J
!
      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, M
         WRITE(*,9998) ( A( I, J ), J = 1, N )
      END DO
!
 9998 FORMAT( 11(:,1X,'(',F6.2,',',F6.2,')') )
      RETURN
END SUBROUTINE
!
!     Auxiliary routine: printing a real matrix.
!
SUBROUTINE PRINT_RMATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      DOUBLE PRECISION A( LDA, * )
!
      INTEGER          I, J
!
      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, M
         WRITE(*,9998) ( A( I, J ), J = 1, N )
      END DO
!
 9998 FORMAT( 11(:,1X,F6.2) )
      RETURN
END SUBROUTINE
