!Let us add BLS and EEBLS routine
!Downloaded from https://konkoly.hu/staff/kovacs/bls_code.html
!Modified to f90 to run in pyaneti
!Converted to f90 with https://fortran.uk/
!https://fortran.uk/plusfortonline.php?id=140&id1=3&id2=1&id3=1&id4=-1&id5=1
!*==EEBLS.spg  processed by SPAG 6.72Dc at 10:53 on 12 Jul 2019
subroutine EEBLS(t,x,freqmin,df,qmi,qma,p,bper,bpow,     &
                     & depth,qtran,in1,in2,n,nf,nb)
use constants
!
!------------------------------------------------------------------------
!     >>>>>>>>>>>> This routine computes BLS spectrum <<<<<<<<<<<<<<
!
!         [ see Kovacs, Zucker & Mazeh 2002, A&A, Vol. 391, 369 ]
!
!     This is the slightly modified version of the original BLS routine
!     by considering Edge Effect (EE) as suggested byTOI-128_HARPS_DRS_RVs.dat
!     Peter R. McCullough [ pmcc@stsci.edu ].
!
!     This modification was motivated by considering the cases when
!     the low state (the transit event) happened to be devided between
!     the first and last bins. In these rare cases the original BLS
!     yields lower detection efficiency because of the lower number of
!     data points in the bin(s) covering the low state.
!
!     For further comments/tests see  www.konkoly.hu/staff/kovacs.html
!------------------------------------------------------------------------
!
!     Input parameters:
!     ~~~~~~~~~~~~~~~~~
!
!     n    = number of data points
!     t    = array {t(i)}, containing the time values of the time series
!     x    = array {x(i)}, containing the data values of the time series
!     u    = temporal/work/dummy array, must be dimensioned in the
!            calling program in the same way as  {t(i)}
!     v    = the same as  {u(i)}
!     nf   = number of frequency points in which the spectrum is computed
!     freqmin = minimum frequency (MUST be > 0)
!     df   = frequency step
!     nb   = number of bins in the folded time series at any test period
!     qmi  = minimum fractional transit length to be tested
!     qma  = maximum fractional transit length to be tested
!
!     Output parameters:
!     ~~~~~~~~~~~~~~~~~~
!
!     p    = array {p(i)}, containing the values of the BLS spectrum
!            at the i-th frequency value -- the frequency values are
!            computed as  f = freqmin + (i-1)*df
!     bper = period at the highest peak in the frequency spectrum
!     bpow = value of {p(i)} at the highest peak
!     depth= depth of the transit at   *bper*
!     qtran= fractional transit length  [ T_transit/bper ]
!     in1  = bin index at the start of the transit [ 0 < in1 < nb+1 ]
!     in2  = bin index at the end   of the transit [ 0 < in2 < nb+1 ]
!
!
!     Remarks:
!     ~~~~~~~~
!
!     -- *freqmin* MUST be greater than  *1/total time span*
!     -- *nb*   MUST be lower than  *nbmax*
!     -- Dimensions of arrays {y(i)} and {ibi(i)} MUST be greater than
!        or equal to  *nbmax*.
!     -- The lowest number of points allowed in a single bin is equal
!        to   MAX(minbin,qmi*N),  where   *qmi*  is the minimum transit
!        length/trial period,   *N*  is the total number of data points,
!        *minbin*  is the preset minimum number of the data points per
!        bin.
!
!========================================================================
!
      IMPLICIT NONE
      !In/out variables
      integer, intent(in) :: n, nf, nb
      real(kind=mireal), intent(in), dimension(n) :: t, x
      real(kind=mireal), intent(in) :: freqmin, df, qmi, qma
      !
      integer, intent(out) :: in1, in2
      real(kind=mireal), intent(out), dimension(nf) :: p
      real(kind=mireal), intent(out) :: bper, bpow, depth, qtran
      !Local variables

      real(kind=mireal), dimension(nb) :: y
      real(kind=mireal), dimension(n) ::  u, v
      integer, dimension(nb) :: ibi
      real(kind=mireal) :: f0 , p0 , ph , pow , power ,rn , rn1 , rn3 , s , s3 , t1 , tot
      integer :: i, j, jf , jn1 , jn2 , jnb , k , kk , kkmi , kma , kmi , &
                 minbin, nb1 , nb2 , nbkma, nbmax
!
      minbin = 5
      nbmax = 2000
      IF ( Nb.GT.nbmax ) WRITE (*,*) ' NB > NBMAX !!'
      IF ( Nb.GT.nbmax ) STOP
      tot = T(N) - T(1)
      IF ( freqmin.LT.1.0D0/tot ) WRITE (*,*) ' freqmin < 1/T !!'
      IF ( freqmin.LT.1.0D0/tot ) STOP
!------------------------------------------------------------------------
!
      rn = DFLOAT(N)
      kmi = IDINT(Qmi*DFLOAT(Nb))
      IF ( kmi.LT.1 ) kmi = 1
      kma = IDINT(Qma*DFLOAT(Nb)) + 1
      kkmi = IDINT(rn*Qmi)
      IF ( kkmi.LT.minbin ) kkmi = minbin
      Bpow = 0.0D0
!
!     The following variables are defined for the extension
!     of arrays  ibi()  and  y()  [ see below ]
!
      nb1 = Nb + 1
      nbkma = Nb + kma
!
!=================================
!     Set temporal time series
!=================================
!
      s = 0.0D0
      t1 = T(1)
      DO i = 1 , N
         U(i) = T(i) - t1
         s = s + X(i)
      ENDDO
      s = s/rn
      DO i = 1 , N
         V(i) = X(i) - s
      ENDDO
!
!******************************
!     Start period search     *
!******************************
!

!THIS PART CAN BE PARALLELIZED
      DO jf = 1 , Nf
         f0 = freqmin + Df*DFLOAT(jf-1)
         p0 = 1.0D0/f0
!
!======================================================
!     Compute folded time series with  *p0*  period
!======================================================
!
         DO j = 1 , Nb
            y(j) = 0.d0
            ibi(j) = 0
         ENDDO
!
         DO i = 1 , N
            ph = U(i)*f0
            ph = ph - IDINT(ph)
            j = 1 + IDINT(Nb*ph)
            ibi(j) = ibi(j) + 1
            y(j) = y(j) + V(i)
         ENDDO
!
!-----------------------------------------------
!     Extend the arrays  ibi()  and  y() beyond
!     nb   by  wrapping
!
         DO j = nb1 , nbkma
            jnb = j - Nb
            ibi(j) = ibi(jnb)
            y(j) = y(jnb)
         ENDDO
!-----------------------------------------------
!
!===============================================
!     Compute BLS statistics for this period
!===============================================
!
         power = 0.0D0
!
         DO i = 1 , Nb
            s = 0.0D0
            k = 0
            kk = 0
            nb2 = i + kma
            DO j = i , nb2
               k = k + 1
               kk = kk + ibi(j)
               s = s + y(j)
               IF ( k.GE.kmi ) THEN
                  IF ( kk.GE.kkmi ) THEN
                     rn1 = DFLOAT(kk)
                     pow = s*s/(rn1*(rn-rn1))
                     IF ( pow.GE.power ) THEN
                        power = pow
                        jn1 = i
                        jn2 = j
                        rn3 = rn1
                        s3 = s
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
!
         power = DSQRT(power)
         P(jf) = power
!
         IF ( power.GE.Bpow ) THEN
            Bpow = power
            In1 = jn1
            In2 = jn2
            Qtran = rn3/rn
            Depth = -s3*rn/(rn3*(rn-rn3))
            Bper = p0
         ENDIF
!
      ENDDO
!
!     Edge correction of transit end index
!
      IF ( In2.GT.Nb ) In2 = In2 - Nb
!
end subroutine