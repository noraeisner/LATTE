subroutine get_QP_gammas(x1,x2,le,lp,P,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1,nx2)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2
  real(kind=mireal), intent(in) :: x1(0:nx1-1), x2(0:nx2-1), le, lp, P
  real(kind=mireal), intent(out), dimension(0:nx1-1,0:nx2-1) ::  gamma_g_g, gamma_dg_dg, gamma_dg_g, gamma_g_dg
  !
  real(kind=mireal), dimension(0:nx1-1,0:nx2-1) :: titj, phi, sinphi

  call fcdist(x1,x2,titj,nx1,nx2)

  phi = 2.*pi*titj/P

  sinphi = sin(phi)

  gamma_g_g = - (sin(phi/2.))**2/2./lp**2 - titj**2/2./le**2
  gamma_g_g = exp(gamma_g_g)

  gamma_g_dg = - pi*sinphi/(2.*P*lp**2) - titj/le**2
  gamma_g_dg = gamma_g_g * gamma_g_dg

  gamma_dg_g = - gamma_g_dg

  gamma_dg_dg = - pi**2*sinphi**2/(4.*P**2*lp**4) &
                + pi**2*cos(phi)/P**2/lp**2         &
                - phi*sinphi/(2.*lp**2*le**2)     &
                - titj**2/le**4                     &
                + 1./le**2
  gamma_dg_dg = gamma_g_g * gamma_dg_dg

end subroutine get_QP_gammas

!This subroutine computes the multi-dimensional GP framework as it appears in
!Rajpaul et al., (2015, https://academic.oup.com/mnras/article/452/3/2269/1079217)
!for a Quasi-Periodic Kernel
!Note: There is an error in eq. (21) of Rajpaul et al., (2015), corrected here
subroutine R15MultiKernel(pars,x1,x2,cov,nx1,nx2)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2
  real(kind=mireal), intent(in) :: pars(0:8) !There are only two parameters for this kernel
  real(kind=mireal), intent(in) :: x1(0:nx1-1), x2(0:nx2-1)
  real(kind=mireal), intent(out) :: cov(0:nx1-1,0:nx2-1)
  !
  integer, parameter :: m = 3
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: gamma_g_g, gamma_dg_dg, gamma_dg_g, gamma_g_dg
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k11, k12, k13
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k21, k22, k23
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k31, k32, k33
  real(kind=mireal) :: Vc, Vr, Lc, Lr, Bc, Br
  real(kind=mireal) :: lp, le, P

  Vc = pars(0); Vr = pars(1); Lc = pars(2); Lr = pars(3); Bc = pars(4); Br = pars(5); le = pars(6); lp = pars(7); P  = pars(8)

  call get_QP_gammas(x1(0:nx1/m-1),x2(0:nx2/m-1),le,lp,P,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1/m,nx2/m)

  k11 = Vc**2 * gamma_g_g + Vr**2 * gamma_dg_dg

  k33 = Bc**2 * gamma_g_g + Br**2 * gamma_dg_dg

  k22 = Lc**2 * gamma_g_g !+ Lr**2 * gamma_dg_dg

  k12 = Vc*Lc*gamma_g_g + Vr*Lc*gamma_dg_g

  k13 = Vc*Bc*gamma_g_g + Vr*Br*gamma_dg_dg + (Vc*Br - Vr*Bc)*gamma_g_dg

  k23 = Lc*Bc*gamma_g_g + Lc*Br*gamma_g_dg

  if (nx1 == nx2) then
    k21 = transpose(k12)
    k31 = transpose(k13)
    k32 = transpose(k23)
  else
    k21 = Vc*Lc*gamma_g_g + Vr*Lc*gamma_g_dg
    k31 = Vc*Bc*gamma_g_g + (Vr*Bc-Vc*Br)*gamma_g_dg + Vr*Br*gamma_dg_dg
    k32 = Bc*Lc*gamma_g_g - Lc*Br*gamma_g_dg
  end if


  cov(0*nx1/m:1*nx1/m-1,0*nx2/m:1*nx2/m-1) = k11
  cov(0*nx1/m:1*nx1/m-1,1*nx2/m:2*nx2/m-1) = k12
  cov(0*nx1/m:1*nx1/m-1,2*nx2/m:3*nx2/m-1) = k13
  !
  cov(1*nx1/m:2*nx1/m-1,0*nx2/m:1*nx2/m-1) = k21
  cov(1*nx1/m:2*nx1/m-1,1*nx2/m:2*nx2/m-1) = k22
  cov(1*nx1/m:2*nx1/m-1,2*nx2/m:3*nx2/m-1) = k23
  !
  cov(2*nx1/m:3*nx1/m-1,0*nx2/m:1*nx2/m-1) = k31
  cov(2*nx1/m:3*nx1/m-1,1*nx2/m:2*nx2/m-1) = k32
  cov(2*nx1/m:3*nx1/m-1,2*nx2/m:3*nx2/m-1) = k33


end subroutine R15MultiKernel
!

!This subroutine computes the multi-dimensional GP matrix for one time-serie
!With a QP kernel
subroutine MultiQP1(pars,x1,x2,cov,nx1,nx2)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2
  real(kind=mireal), intent(in) :: pars(0:4) !There are only two parameters for this kernel
  real(kind=mireal), intent(in) :: x1(0:nx1-1)
  real(kind=mireal), intent(in) :: x2(0:nx2-1)
  real(kind=mireal), intent(out) :: cov(0:nx1-1,0:nx2-1)
  !A = pars(0), lambda_p = pars(1), lambda_q = pars(2), P = pars(3)
  !
  real(kind=mireal), dimension(0:nx1-1,0:nx2-1) :: titj, phi, sinphi, gamma_dg_dg, gamma_g_g
  real(kind=mireal) :: Ac, Ar, lp, le, P

  Ac = pars(0)
  Ar = pars(1)
  le = pars(2)
  lp = pars(3)
  P  = pars(4)

  !Get the x_i - x_j
  call fcdist(x1,x2,titj,nx1,nx2)

  phi = 2.*pi*titj/P

  sinphi = sin(phi)

  gamma_g_g = - (sin(phi/2.))**2/2./lp**2 - titj**2/2./le**2
  gamma_g_g = exp(gamma_g_g)

  gamma_dg_dg = - pi**2*sinphi**2/(4.*P**2*lp**4) &
                + pi**2*cos(phi)/P**2/lp**2         &
                - phi*sinphi/(2.*lp**2*le**2)     &
                - titj**2/le**4                     &
                + 1./le**2

  gamma_dg_dg = gamma_g_g * gamma_dg_dg

  cov = Ac*Ac*gamma_g_g + Ar*Ar*gamma_dg_dg

end subroutine MultiQP1

!This subroutine computes the multi-dimensional GP matrix for two time-series
!With a QP kernel
subroutine MultiQP2(pars,x1,x2,cov,nx1,nx2)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2
  real(kind=mireal), intent(in) :: pars(0:8) !There are only two parameters for this kernel
  real(kind=mireal), intent(in) :: x1(0:nx1-1), x2(0:nx2-1)
  real(kind=mireal), intent(out) :: cov(0:nx1-1,0:nx2-1)
  !
  integer, parameter :: m = 2
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: gamma_g_g, gamma_dg_dg, gamma_dg_g, gamma_g_dg
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k11, k12
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k21, k22
  real(kind=mireal) :: Vc, Vr, Lc, Lr
  real(kind=mireal) :: lp, le, P

  Vc = pars(0); Vr = pars(1); Lc = pars(2); Lr = pars(3); le = pars(4); lp = pars(5); P  = pars(6)

  call get_QP_gammas(x1(0:nx1/m-1),x2(0:nx2/m-1),le,lp,P,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1/m,nx2/m)

  k11 = Vc**2 * gamma_g_g + Vr**2 * gamma_dg_dg

  k22 = Lc**2 * gamma_g_g + Lr**2 * gamma_dg_dg

  k12 = Vc*Lc*gamma_g_g + Vr*Lr*gamma_dg_dg + (Vc*Lr - Vr*Lc)*gamma_g_dg

  if (nx1 == nx2) then
    k21 = transpose(k12)
  else
    k21 = Lc*Vc*gamma_g_g + Lr*Vr*gamma_dg_dg + (Lc*Vr - Lr*Vc)*gamma_g_dg
  end if

  !
  cov(0*nx1/m:1*nx1/m-1,0*nx2/m:1*nx2/m-1) = k11
  cov(0*nx1/m:1*nx1/m-1,1*nx2/m:2*nx2/m-1) = k12
  !
  cov(1*nx1/m:2*nx1/m-1,0*nx2/m:1*nx2/m-1) = k21
  cov(1*nx1/m:2*nx1/m-1,1*nx2/m:2*nx2/m-1) = k22
  !

end subroutine MultiQP2

!
subroutine MultiQP3(pars,x1,x2,cov,nx1,nx2)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2
  real(kind=mireal), intent(in) :: pars(0:8) !There are only two parameters for this kernel
  real(kind=mireal), intent(in) :: x1(0:nx1-1), x2(0:nx2-1)
  real(kind=mireal), intent(out) :: cov(0:nx1-1,0:nx2-1)
  !
  integer, parameter :: m = 3
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: gamma_g_g, gamma_dg_dg, gamma_dg_g, gamma_g_dg
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k11, k12, k13
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k21, k22, k23
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k31, k32, k33
  real(kind=mireal) :: a1, b1, a2, b2, a3, b3
  real(kind=mireal) :: lp, le, P

  a1 = pars(0); b1 = pars(1); a2 = pars(2); b2 = pars(3); a3 = pars(4); b3 = pars(5); le = pars(6); lp = pars(7); P  = pars(8)

  call get_QP_gammas(x1(0:nx1/m-1),x2(0:nx2/m-1),le,lp,P,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1/m,nx2/m)

  k11 = a1**2 * gamma_g_g + b1**2 * gamma_dg_dg

  k22 = a2**2 * gamma_g_g + b2**2 * gamma_dg_dg

  k33 = a3**2 * gamma_g_g + b3**2 * gamma_dg_dg

  k12 = a1*a2*gamma_g_g + b1*b2*gamma_dg_dg + (a1*b2 - a2*b1)*gamma_g_dg

  k13 = a1*a3*gamma_g_g + b1*b3*gamma_dg_dg + (a1*b3 - a3*b1)*gamma_g_dg

  k23 = a2*a3*gamma_g_g + b2*b3*gamma_dg_dg + (a2*b3 - a3*b2)*gamma_g_dg

  if (nx1 == nx2) then
    k21 = transpose(k12)
    k31 = transpose(k13)
    k32 = transpose(k23)
  else
  k21 = a2*a1*gamma_g_g + b2*b1*gamma_dg_dg + (a2*b1 - a1*b2)*gamma_g_dg
  k31 = a3*a1*gamma_g_g + b3*b1*gamma_dg_dg + (a3*b1 - a1*b3)*gamma_g_dg
  k32 = a3*a2*gamma_g_g + b3*b2*gamma_dg_dg + (a3*b2 - a2*b3)*gamma_g_dg
  end if


  cov(0*nx1/m:1*nx1/m-1,0*nx2/m:1*nx2/m-1) = k11
  cov(0*nx1/m:1*nx1/m-1,1*nx2/m:2*nx2/m-1) = k12
  cov(0*nx1/m:1*nx1/m-1,2*nx2/m:3*nx2/m-1) = k13
  !
  cov(1*nx1/m:2*nx1/m-1,0*nx2/m:1*nx2/m-1) = k21
  cov(1*nx1/m:2*nx1/m-1,1*nx2/m:2*nx2/m-1) = k22
  cov(1*nx1/m:2*nx1/m-1,2*nx2/m:3*nx2/m-1) = k23
  !
  cov(2*nx1/m:3*nx1/m-1,0*nx2/m:1*nx2/m-1) = k31
  cov(2*nx1/m:3*nx1/m-1,1*nx2/m:2*nx2/m-1) = k32
  cov(2*nx1/m:3*nx1/m-1,2*nx2/m:3*nx2/m-1) = k33


end subroutine MultiQP3
!!
subroutine MultiQP4(pars,x1,x2,cov,nx1,nx2)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2
  real(kind=mireal), intent(in) :: pars(0:10)
  real(kind=mireal), intent(in) :: x1(0:nx1-1), x2(0:nx2-1)
  real(kind=mireal), intent(out) :: cov(0:nx1-1,0:nx2-1)
  !
  integer, parameter :: m = 4
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: gamma_g_g, gamma_dg_dg, gamma_dg_g, gamma_g_dg
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k11, k12, k13, k14
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k21, k22, k23, k24
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k31, k32, k33, k34
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k41, k42, k43, k44
  real(kind=mireal) :: a1, b1, a2, b2, a3, b3, a4, b4
  real(kind=mireal) :: lp, le, P

  a1 = pars(0); b1 = pars(1); a2 = pars(2); b2 = pars(3); a3 = pars(4); b3 = pars(5); a4 = pars(6); b4 = pars(7)
  le = pars(8); lp = pars(9); P  = pars(10)

  call get_QP_gammas(x1(0:nx1/m-1),x2(0:nx2/m-1),le,lp,P,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1/m,nx2/m)

  k11 = a1**2 * gamma_g_g + b1**2 * gamma_dg_dg

  k22 = a2**2 * gamma_g_g + b2**2 * gamma_dg_dg

  k33 = a3**2 * gamma_g_g + b3**2 * gamma_dg_dg

  k44 = a4**2 * gamma_g_g + b4**2 * gamma_dg_dg

  k12 = a1*a2*gamma_g_g + b1*b2*gamma_dg_dg + (a1*b2 - a2*b1)*gamma_g_dg

  k13 = a1*a3*gamma_g_g + b1*b3*gamma_dg_dg + (a1*b3 - a3*b1)*gamma_g_dg

  k14 = a1*a4*gamma_g_g + b1*b4*gamma_dg_dg + (a1*b4 - a4*b1)*gamma_g_dg

  k23 = a2*a3*gamma_g_g + b2*b3*gamma_dg_dg + (a2*b3 - a3*b2)*gamma_g_dg

  k24 = a2*a4*gamma_g_g + b2*b4*gamma_dg_dg + (a2*b4 - a4*b2)*gamma_g_dg

  k34 = a3*a4*gamma_g_g + b3*b4*gamma_dg_dg + (a3*b4 - a4*b3)*gamma_g_dg

  if (nx1 == nx2) then
    k21 = transpose(k12)
    k31 = transpose(k13)
    k32 = transpose(k23)
    k41 = transpose(k14)
    k42 = transpose(k24)
    k43 = transpose(k34)
  else
    k21 = a2*a1*gamma_g_g + b2*b1*gamma_dg_dg + (a2*b1 - a1*b2)*gamma_g_dg
    k31 = a3*a1*gamma_g_g + b3*b1*gamma_dg_dg + (a3*b1 - a1*b3)*gamma_g_dg
    k32 = a3*a2*gamma_g_g + b3*b2*gamma_dg_dg + (a3*b2 - a2*b3)*gamma_g_dg
    k41 = a4*a1*gamma_g_g + b4*b1*gamma_dg_dg + (a4*b1 - a1*b4)*gamma_g_dg
    k42 = a4*a2*gamma_g_g + b4*b2*gamma_dg_dg + (a4*b2 - a2*b4)*gamma_g_dg
    k43 = a4*a3*gamma_g_g + b4*b3*gamma_dg_dg + (a4*b3 - a3*b4)*gamma_g_dg
  end if


  !
  !Time to fill the covariance matrix
  !
  cov(0*nx1/m:1*nx1/m-1,0*nx2/m:1*nx2/m-1) = k11
  cov(0*nx1/m:1*nx1/m-1,1*nx2/m:2*nx2/m-1) = k12
  cov(0*nx1/m:1*nx1/m-1,2*nx2/m:3*nx2/m-1) = k13
  cov(0*nx1/m:1*nx1/m-1,3*nx2/m:4*nx2/m-1) = k14
  !
  cov(1*nx1/m:2*nx1/m-1,0*nx2/m:1*nx2/m-1) = k21
  cov(1*nx1/m:2*nx1/m-1,1*nx2/m:2*nx2/m-1) = k22
  cov(1*nx1/m:2*nx1/m-1,2*nx2/m:3*nx2/m-1) = k23
  cov(1*nx1/m:2*nx1/m-1,3*nx2/m:4*nx2/m-1) = k24
  !
  cov(2*nx1/m:3*nx1/m-1,0*nx2/m:1*nx2/m-1) = k31
  cov(2*nx1/m:3*nx1/m-1,1*nx2/m:2*nx2/m-1) = k32
  cov(2*nx1/m:3*nx1/m-1,2*nx2/m:3*nx2/m-1) = k33
  cov(2*nx1/m:3*nx1/m-1,3*nx2/m:4*nx2/m-1) = k34
  !
  cov(3*nx1/m:4*nx1/m-1,0*nx2/m:1*nx2/m-1) = k41
  cov(3*nx1/m:4*nx1/m-1,1*nx2/m:2*nx2/m-1) = k42
  cov(3*nx1/m:4*nx1/m-1,2*nx2/m:3*nx2/m-1) = k43
  cov(3*nx1/m:4*nx1/m-1,3*nx2/m:4*nx2/m-1) = k44
  !


end subroutine MultiQP4
!
subroutine MultiQP5(pars,x1,x2,cov,nx1,nx2)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2
  real(kind=mireal), intent(in) :: pars(0:12)
  real(kind=mireal), intent(in) :: x1(0:nx1-1), x2(0:nx2-1)
  real(kind=mireal), intent(out) :: cov(0:nx1-1,0:nx2-1)
  !
  integer, parameter :: m = 5
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: gamma_g_g, gamma_dg_dg, gamma_dg_g, gamma_g_dg
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k11, k12, k13, k14, k15
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k21, k22, k23, k24, k25
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k31, k32, k33, k34, k35
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k41, k42, k43, k44, k45
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k51, k52, k53, k54, k55
  real(kind=mireal) :: a1, b1, a2, b2, a3, b3, a4, b4, a5, b5
  real(kind=mireal) :: lp, le, P

  a1 = pars(0); b1 = pars(1); a2 = pars(2); b2 = pars(3); a3 = pars(4); b3 = pars(5); a4 = pars(6); b4 = pars(7)
  a5 = pars(8); b5 = pars(9)
  le = pars(10); lp = pars(11); P  = pars(12)

  call get_QP_gammas(x1(0:nx1/m-1),x2(0:nx2/m-1),le,lp,P,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1/m,nx2/m)

  k11 = a1**2 * gamma_g_g + b1**2 * gamma_dg_dg

  k22 = a2**2 * gamma_g_g + b2**2 * gamma_dg_dg

  k33 = a3**2 * gamma_g_g + b3**2 * gamma_dg_dg

  k44 = a4**2 * gamma_g_g + b4**2 * gamma_dg_dg

  k55 = a5**2 * gamma_g_g + b5**2 * gamma_dg_dg

  k12 = a1*a2*gamma_g_g + b1*b2*gamma_dg_dg + (a1*b2 - a2*b1)*gamma_g_dg

  k13 = a1*a3*gamma_g_g + b1*b3*gamma_dg_dg + (a1*b3 - a3*b1)*gamma_g_dg

  k14 = a1*a4*gamma_g_g + b1*b4*gamma_dg_dg + (a1*b4 - a4*b1)*gamma_g_dg

  k15 = a1*a5*gamma_g_g + b1*b5*gamma_dg_dg + (a1*b5 - a5*b1)*gamma_g_dg

  k23 = a2*a3*gamma_g_g + b2*b3*gamma_dg_dg + (a2*b3 - a3*b2)*gamma_g_dg

  k24 = a2*a4*gamma_g_g + b2*b4*gamma_dg_dg + (a2*b4 - a4*b2)*gamma_g_dg

  k25 = a2*a5*gamma_g_g + b2*b5*gamma_dg_dg + (a2*b5 - a5*b2)*gamma_g_dg

  k34 = a3*a4*gamma_g_g + b3*b4*gamma_dg_dg + (a3*b4 - a4*b3)*gamma_g_dg

  k35 = a3*a5*gamma_g_g + b3*b5*gamma_dg_dg + (a3*b5 - a5*b3)*gamma_g_dg

  k45 = a4*a5*gamma_g_g + b4*b5*gamma_dg_dg + (a4*b5 - a5*b4)*gamma_g_dg

  if (nx1 == nx2) then
    k21 = transpose(k12)
    k31 = transpose(k13)
    k32 = transpose(k23)
    k41 = transpose(k14)
    k42 = transpose(k24)
    k43 = transpose(k34)
    k51 = transpose(k15)
    k52 = transpose(k25)
    k53 = transpose(k35)
    k54 = transpose(k45)
  else
    k21 = a2*a1*gamma_g_g + b2*b1*gamma_dg_dg + (a2*b1 - a1*b2)*gamma_g_dg
    k31 = a3*a1*gamma_g_g + b3*b1*gamma_dg_dg + (a3*b1 - a1*b3)*gamma_g_dg
    k32 = a3*a2*gamma_g_g + b3*b2*gamma_dg_dg + (a3*b2 - a2*b3)*gamma_g_dg
    k41 = a4*a1*gamma_g_g + b4*b1*gamma_dg_dg + (a4*b1 - a1*b4)*gamma_g_dg
    k42 = a4*a2*gamma_g_g + b4*b2*gamma_dg_dg + (a4*b2 - a2*b4)*gamma_g_dg
    k43 = a4*a3*gamma_g_g + b4*b3*gamma_dg_dg + (a4*b3 - a3*b4)*gamma_g_dg
    k51 = a5*a1*gamma_g_g + b5*b1*gamma_dg_dg + (a5*b1 - a1*b5)*gamma_g_dg
    k52 = a5*a2*gamma_g_g + b5*b2*gamma_dg_dg + (a5*b2 - a2*b5)*gamma_g_dg
    k53 = a5*a3*gamma_g_g + b5*b3*gamma_dg_dg + (a5*b3 - a3*b5)*gamma_g_dg
    k54 = a5*a4*gamma_g_g + b5*b4*gamma_dg_dg + (a5*b4 - a4*b5)*gamma_g_dg
  end if
  !
  !Time to fill the covariance matrix
  !
  cov(0*nx1/m:1*nx1/m-1,0*nx2/m:1*nx2/m-1) = k11
  cov(0*nx1/m:1*nx1/m-1,1*nx2/m:2*nx2/m-1) = k12
  cov(0*nx1/m:1*nx1/m-1,2*nx2/m:3*nx2/m-1) = k13
  cov(0*nx1/m:1*nx1/m-1,3*nx2/m:4*nx2/m-1) = k14
  cov(0*nx1/m:1*nx1/m-1,4*nx2/m:5*nx2/m-1) = k15
  !
  cov(1*nx1/m:2*nx1/m-1,0*nx2/m:1*nx2/m-1) = k21
  cov(1*nx1/m:2*nx1/m-1,1*nx2/m:2*nx2/m-1) = k22
  cov(1*nx1/m:2*nx1/m-1,2*nx2/m:3*nx2/m-1) = k23
  cov(1*nx1/m:2*nx1/m-1,3*nx2/m:4*nx2/m-1) = k24
  cov(1*nx1/m:2*nx1/m-1,4*nx2/m:5*nx2/m-1) = k25
  !
  cov(2*nx1/m:3*nx1/m-1,0*nx2/m:1*nx2/m-1) = k31
  cov(2*nx1/m:3*nx1/m-1,1*nx2/m:2*nx2/m-1) = k32
  cov(2*nx1/m:3*nx1/m-1,2*nx2/m:3*nx2/m-1) = k33
  cov(2*nx1/m:3*nx1/m-1,3*nx2/m:4*nx2/m-1) = k34
  cov(2*nx1/m:3*nx1/m-1,4*nx2/m:5*nx2/m-1) = k35
  !
  cov(3*nx1/m:4*nx1/m-1,0*nx2/m:1*nx2/m-1) = k41
  cov(3*nx1/m:4*nx1/m-1,1*nx2/m:2*nx2/m-1) = k42
  cov(3*nx1/m:4*nx1/m-1,2*nx2/m:3*nx2/m-1) = k43
  cov(3*nx1/m:4*nx1/m-1,3*nx2/m:4*nx2/m-1) = k44
  cov(3*nx1/m:4*nx1/m-1,4*nx2/m:5*nx2/m-1) = k45
  !
  cov(4*nx1/m:5*nx1/m-1,0*nx2/m:1*nx2/m-1) = k51
  cov(4*nx1/m:5*nx1/m-1,1*nx2/m:2*nx2/m-1) = k52
  cov(4*nx1/m:5*nx1/m-1,2*nx2/m:3*nx2/m-1) = k53
  cov(4*nx1/m:5*nx1/m-1,3*nx2/m:4*nx2/m-1) = k54
  cov(4*nx1/m:5*nx1/m-1,4*nx2/m:5*nx2/m-1) = k55
  !


end subroutine MultiQP5
!
subroutine QPMultiKernel(pars,x1,x2,cov,nx1,nx2,ndim)
  use constants
  implicit none
  !The ndim variable is defined inside constants.f90
  !
  integer, intent(in) :: nx1, nx2, ndim
  real(kind=mireal), intent(in) :: pars(0:ndim*2+3-1) !m*2 amplitudes, 3 parameters for the QP kernel
  real(kind=mireal), intent(in) :: x1(0:nx1-1), x2(0:nx2-1)
  real(kind=mireal), intent(out) :: cov(0:nx1-1,0:nx2-1)
  !
  real(kind=mireal), dimension(0:nx1/ndim-1,0:nx2/ndim-1) :: gamma_g_g, gamma_dg_dg, gamma_g_dg, gamma_dg_g
  real(kind=mireal), dimension(0:nx1/ndim-1,0:nx2/ndim-1,0:ndim-1,0:ndim-1) :: kas
  real(kind=mireal) :: Amplitudes(0:ndim*2-1)
  real(kind=mireal) :: le, lp, P
  integer :: i, j

  !Read the parameters
  Amplitudes(:) = pars(0:ndim*2-1)
  le = pars(ndim*2)
  lp = pars(ndim*2+1)
   P = pars(ndim*2+2)

  call get_QP_gammas(x1(0:nx1/ndim-1),x2(0:nx2/ndim-1),le,lp,P,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,&
                     nx1/ndim,nx2/ndim)

  if ( nx1 .ne. nx2  ) then !Compute the K's for not squared matrices

  do i = 0, ndim - 1
    do j = 0, ndim - 1
      kas(:,:,i,j) = amplitudes(i*2)*amplitudes(j*2)*gamma_g_g + &
                     amplitudes(i*2+1)*amplitudes(j*2+1)*gamma_dg_dg + &
                     ( amplitudes(i*2)*amplitudes(j*2+1) - amplitudes(i*2+1)*amplitudes(j*2) )* gamma_g_dg
      cov(i*nx1/ndim:(i+1)*nx1/ndim-1,j*nx2/ndim:(j+1)*nx2/ndim-1) = kas(:,:,i,j)
    end do
  end do

  else !Compute the K's for square matrices

  do i = 0, ndim - 1
    do j = i, ndim - 1
      kas(:,:,i,j) = amplitudes(i*2)*amplitudes(j*2)*gamma_g_g + &
                     amplitudes(i*2+1)*amplitudes(j*2+1)*gamma_dg_dg + &
                     ( amplitudes(i*2)*amplitudes(j*2+1) - amplitudes(i*2+1)*amplitudes(j*2) )* gamma_g_dg
      cov(i*nx1/ndim:(i+1)*nx1/ndim-1,j*nx2/ndim:(j+1)*nx2/ndim-1) = kas(:,:,i,j)
    end do
  end do

  do i = 1, ndim - 1
    do j = 0, i - 1
      kas(:,:,i,j) = transpose(kas(:,:,j,i))
      cov(i*nx1/ndim:(i+1)*nx1/ndim-1,j*nx2/ndim:(j+1)*nx2/ndim-1) = kas(:,:,i,j)
    end do
  end do

  end if

end subroutine QPMultiKernel
