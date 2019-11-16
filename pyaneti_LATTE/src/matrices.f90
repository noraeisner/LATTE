 !https://ww2.odu.edu/~agodunov/computing/programs/book2/Ch06/Inverse.f90
  subroutine inverse(m,cc,n)
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
  ! the original matrix a(n,n) will be destroyed
  ! during the calculation
!===========================================================
  use constants
  implicit none
    integer, intent(in) :: n
    real(kind=mireal), intent(in) :: m(0:n-1,0:n-1)
    real(kind=mireal), intent(out) :: cc(0:n-1,0:n-1)
    !
    real(kind=mireal) :: a(n,n), c(n,n), L(n,n), U(n,n), b(n), d(n), x(n)
    real(kind=mireal) :: coeff
    integer :: i, j, k

    a = m

    ! step 0: initialization for matrices L and U and b
    ! Fortran 90/95 aloows such operations on matrices
    L=0.0
    U=0.0
    b=0.0

    !  step 1: forward elimination
    do k=1, n-1
      do i=k+1,n
        coeff=a(i,k)/a(k,k)
        L(i,k) = coeff
        do j=k+1,n
           a(i,j) = a(i,j)-coeff*a(k,j)
        end do
      end do
    end do

    ! Step 2: prepare L and U matrices
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i=1,n
      L(i,i) = 1.0
    end do
    ! U matrix is the upper triangular part of A
    do j=1,n
      do i=1,j
         U(i,j) = a(i,j)
      end do
    end do

    ! Step 3: compute columns of the inverse matrix C
    do k=1,n
      b(k)=1.0
      d(1) = b(1)
    ! Step 3a: Solve Ld=b using the forward substitution
      do i=2,n
        d(i)=b(i)
        do j=1,i-1
          d(i) = d(i) - L(i,j)*d(j)
        end do
      end do
    ! Step 3b: Solve Ux=d using the back substitution
      x(n)=d(n)/U(n,n)
      do i = n-1,1,-1
        x(i) = d(i)
        do j=n,i+1,-1
          x(i)=x(i)-U(i,j)*x(j)
        end do
        x(i) = x(i)/u(i,i)
      end do
    ! Step 3c: fill the solutions x(n) into column k of C
      do i=1,n
        c(i,k) = x(i)
      end do
      b(k)=0.0
    end do

    cc = c

  end subroutine inverse

!Subroutine to find the determinant of a square matrix
!Taken from https://github.com/ashwith/workspace/blob/master/fortran/determinant.f90
!Modified as a subroutine
!----------------------------------------------------------------------------------------
!Function to find the determinant of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Description: The subroutine is based on two key points:
!1] A determinant is unaltered when row operations are performed: Hence, using this principle,
!row operations (column operations would work as well) are used
!to convert the matrix into upper traingular form
!2]The determinant of a triangular matrix is obtained by finding the product of the diagonal elements
!
  subroutine findlogddet(a,det,n)
  use constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  real(kind=mireal), DIMENSION(n,n), intent(in) :: a
  real(kind=mireal), intent (out) :: det
  real(kind=mireal) :: m, temp
  real(kind=mireal), DIMENSION(n,n) :: matrix
  INTEGER :: i, j, k, l
  LOGICAL :: DetExists = .TRUE.
  matrix = a
  l = 1
  !Convert to upper triangular form
  DO k = 1, n-1
    IF (matrix(k,k) == 0) THEN
      DetExists = .FALSE.
      DO i = k+1, n
        IF (matrix(i,k) /= 0) THEN
          DO j = 1, n
            temp = matrix(i,j)
            matrix(i,j)= matrix(k,j)
            matrix(k,j) = temp
          END DO
          DetExists = .TRUE.
          l=-l
          EXIT
        ENDIF
      END DO
      IF (DetExists .EQV. .FALSE.) THEN
        det = 0
        return
      END IF
    ENDIF
    DO j = k+1, n
      m = matrix(j,k)/matrix(k,k)
      DO i = k+1, n
        matrix(j,i) = matrix(j,i) - m*matrix(k,i)
      END DO
    END DO
  END DO

  !Calculate determinant by finding product of diagonal elements
  det = log(abs(real(l)))
  DO i = 1, n
    det = det + log(abs(matrix(i,i)))
  END DO

  end subroutine findlogddet

  subroutine fill_diag(v,M,n)
  use constants
  implicit none
  !
  integer, intent(in) :: n
  real(kind=mireal), intent(in) :: v(0:n-1)
  real(kind=mireal), intent(out) :: M(0:n-1,0:n-1)
  !
  integer :: i

  M = 0.d0
  do i = 0, n - 1
    M(i,i) = v(i)
  end do

  end subroutine fill_diag

  ! Returns the inverse of a matrix calculated by finding the LU
  ! decomposition.  Depends on LAPACK.
  ! Taken and modified from http://fortranwiki.org/fortran/show/Matrix+inversion
  subroutine cholesky_inv_det(A,Ainv,detlog,n)
   use constants
    integer, intent(in) :: n
    real(kind=mireal), dimension(n,n), intent(in) :: A
    real(kind=mireal), dimension(n,n), intent(out) :: Ainv
    real(kind=mireal), intent(out) :: detlog
    !
    integer :: info, i

    ! External procedures defined in LAPACK
    external DPOTRF
    external DPOTRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! https://www.math.utah.edu/software/lapack/lapack-d/dpotrf.html
    call DPOTRF('U',n,Ainv,n,info)

    !Now Ainv = L, such that A = L L*
    !The determinant of det A is then det A = det L x det L* = (det L)^2
    !and det L = mult_1^N L(i,i)

    !Since we are interested on log det A (for the likelihood calculation)
    !log(det(A)) = 2 (log (det (L)))

    detlog = 0.d0
    do i = 1, n
      detlog = detlog + log((Ainv(i,i)))
    end do
    detlog = 2.*detlog

    !DPOTRI calculates the inverse matrix using as input the output of DPOTRF
    !http://www.math.utah.edu/software/lapack/lapack-d/dpotri.html
    call DPOTRI('U',n,Ainv,n,info)

    do i = 2, n
  !    do j = 1, i
  !      Ainv(i,j) = Ainv(j,i)
        Ainv(i,1:i) = Ainv(1:i,i)
  !    end do
    end do

  end subroutine cholesky_inv_det

  subroutine mtops(M,V,n)
  use constants
  implicit none
    integer, intent(in) :: n
    real(kind=mireal), dimension(n,n), intent(in) :: M
    real(kind=mireal), dimension(n*(n+1)/2), intent(out) :: V
    !
    integer :: i, j

    do j = 1, n
      do i = 1, j
        V(i-1 + j*(j-1)/2) = M(i,j)
      end do
    end do

  end subroutine mtops


  subroutine pstom(V,M,n)
  use constants
  implicit none
    integer, intent(in) :: n
    real(kind=mireal), dimension(n*(n+1)/2), intent(in) :: V
    real(kind=mireal), dimension(n,n), intent(out) :: M
    !
    integer :: i, j

    M = 0.d0
    do j = 1, n
      do i = 1, j
        M(i,j) = V(i-1 + j*(j-1)/2)
        if (i == j) M(i,j) = M(i,j)/2.
      end do
    end do

    M = M + transpose(M)

  end subroutine pstom


