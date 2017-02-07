! +-------------------------------------------------------+
! |  program for 1D phononic ctystals with poroelasticity
! |  solved bY PME(plain Expansion Method)
! +-------------------------------------------------------+
!===
!   格子点を作成
!   逆格子空間を作成
!   入射波とパラメタを展開
!   展開した係数から行列を作成
!   LAPACKで固有値，固有ベクトルを出す
!   与えた波数に対応した固有値，固有ベクトルを出力
!
program poroerastic1D_imp
    implicit none
    character :: material*4, nr*3
    integer :: iX, iY, l, k, iter, INFO, i
    integer, parameter :: DD=kind(0d0)
    integer, parameter :: nX=150, nY=150, NE=100  ! 要素数,刻み幅
  !--{点の数} = 2*{正方向の要素数}+{原点}
  !--{点の総数} = {X方向の点の数}*{Y方向の点の数}
    integer, parameter :: ngX=2*nX+1
    integer, parameter :: ngX2=ngX*2, ngX3=ngX*3, ngX4=ngX*4, LWORK=ngX4*2 !dimA=ngX4

    real(DD) :: dkX, dkY, kX, kY, kXmin, kYmin, kXmax, kYmax, X, Y
    real(DD) :: gX(ngX)
    real(DD),dimension(1:2) :: rho11, rho12, rho22, P, Q, R, mub
    real(DD) :: rhoS, rhoF, Ks, Kf, Kb, mu, f, tor, filling, bunbo
    real(DD) :: rhoSp, rhoFp, Ksp, Kfp, Kbp, musp, mubp, fp
    real(DD) :: rhoSh, rhoFh, Ksh, Kfh, Kbh, mush, mubh, fh
    real(DD) :: lX, lY, pi, pi2, RWORK(ngX4*8)
    !--多孔質体(シリカ)の波の伝搬速度
    real(DD), parameter :: csl=5.970d3, cst=3.760d3

    complex :: coeff
    complex(DD),dimension(1:ngX, 1:ngX) :: Pg, Qg, Rg, rho11g, rho12g, rho22g, mubg
    complex(DD),dimension(1:ngX4, 1:ngX4) :: A, B, VL, VR
    complex(DD) :: eigen(ngX4), ALPHA(ngX4), BETA(ngX4), WORK(LWORK)
    complex(DD), parameter :: COM=dcmplx(0.0d0,1.0d0)

    external coeff, DEIGSRT

do i=1,8
write(nr,'(I3.3)') i
!==================================================================
! +-------+
! | INIT. |
! +-------+
!===
    open  (11, file='input1D_'//nr//'.dat', status='unknown')
    !--フォノニック結晶の構造
      read (11,*) lX              ! 格子定数[m]
      read (11,*) filling             ! 充填率　
                                      ! dA=filling*lX       :A層の厚さ
                                      ! dB=(1-filling)*lX   :B層の厚さ
    !--多孔質体(poroelastic)の物質パラメタ
      read (11,*) rhoSp, rhoFp        ! 多孔質体での密度(固体，流体)[kg m-3]
      read (11,*) Ksp, Kfp            ! 体積弾性率（固体，流体）[Pa]
      read (11,*) Kbp, mubp           ! バルクの弾性定数（体積弾性率，せん断弾性率）[Pa]
      read (11,*) fp                  ! 孔隙率
    !--基盤(Host)の物質パラメタ
      read (11,*) rhoSh, rhoFh        ! 多孔質体での密度(固体，流体)[kg m-3]
      read (11,*) Ksh, Kfh            ! 体積弾性率（固体，流体）[Pa]
      read (11,*) Kbh, mubh           ! バルクの弾性定数（体積弾性率，せん断弾性率）[Pa]
      read (11,*) fh                  ! 孔隙率
    close (11)

    pi=4.0d0*atan(1.0d0)
    pi2=2.0d0*pi
!===

!==================================================================
! +--------------------------------+
! | 逆格子ベクトルG1(G),G2(G")の作成 |
! +--------------------------------+
!===
  !--{第一ブリルアンゾーンの境界} = {逆格子ベクトルの辺垂直二等分線}
  !--{刻み幅} = {格子の辺長}/{刻み数}
    kXmin = 0.0d0
    kXmax = pi/lX
    dkX = (kXmax - kXmin) / NE
    lY = lX
    kYmin = 0.0d0
    kYmax = pi/lY
    dkY = (kYmax - kYmin) / NE

  !--{位相空間のベクトルGi(ITER)}= 2*{pi}/{辺長Li} * ITER
    iter=0
    do iX = -nX, nX
     iter = iter+1
     gX(iter) = (pi2/lX)*dble(iX)
    end do
!===

!==================================================================
! +--------------------------------------------+
! | 多孔質体中の物理パラメータを定める  　　　　　 |
! +--------------------------------------------+
!===
  open (25, file='elastic1Dimp_'//nr//'.dat', status='unknown')
  do iter = 1,2
  !--多孔質体中では
    if ( iter .eq. 1 ) then
      f = fp
      rhoS = rhoSp
      rhoF = rhoFp
      Ks = Ksp
      Kf = Kfp
      Kb = Kbp
      mu = mubp
      material = "poro"
      tor=f**(-2.0d0/3.0d0)           ! 迷路度
  !--基盤中では
    else
      f = fh
      rhoS = rhoSh
      rhoF = rhoFh
      Ks = Ksh
      Kf = Kfh
      Kb = Kbh
      mu = mubh
      material = "host"
!********************************************************************
    !--{tor}={f^-2/3}の定義より，ホスト基盤が
    !--f=0(すべて固体)の時に発散しないようにする必要がある
    !--迷路度の定義の妥当性が不明だが元の論文が$30のため保留
      if (f .eq. 0d0) then
        tor=1                                     ! 発散防止
      else
        tor=f**(-2.0d0/3.0d0)                     ! 迷路度
      end if
!*******************************************************************
    end if
    mub(iter) = mu
    rho11(iter) = (1.0d0-f)*rhoS + (tor-1.0d0)*f*rhoF
    rho12(iter) = -(tor-1.0d0)*f*rhoF
    rho22(iter) = tor*f*rhoF
    bunbo=1.0d0-f-Kb/Ks + f*Ks/Kf
    P(iter) = Ks*( (1.0d0-f)*(1.0d0-f-Kb/Ks) + f*Kb/Kf )/bunbo + 4.0d0*mu/3.0d0
    Q(iter) = f*Ks*(1.0d0-f-Kb/Ks)/bunbo
    R(iter) = (f**2 *Ks)/bunbo
    write(25,*) material
    write(25,*) 'rhoS =', rhoS
    write(25,*) 'rhoF =', rhoF
    write(25,*) 'Ks =', Ks
    write(25,*) 'Kf =', Kf
    write(25,*) 'Kb =', Kb
    write(25,*) 'mub =', mu
    write(25,*) 'f =', f
    write(25,*) 'tor =',tor
  end do
    write(25,*) 'P =',P
    write(25,*) 'Q =',Q
    write(25,*) 'R =',R
    write(25,*) 'rho11 =', rho11
    write(25,*) 'rho12 =', rho12
    write(25,*) 'rho22 =', rho22
    write(25,*) 'mub =', mub
  close(25)
!===

!==================================================================
! +---------------------------------------------+
! | 固有方程式のための行列A,Bを作成　　　　　　　　 |
! | {固有値}*{B}.{固有ベクトル}={A}.{固有ベクトル} |
! | ZGGEVを使って固有値求める                     |
! +---------------------------------------------+
!===
  !--係数のマトリクスを作成
  open(10,file='1Dimp_Re_'//nr//'.dat')
  open(20,file='1Dimp_Im_'//nr//'.dat')
    do l=1,ngX
    do k=1,ngX
       Pg(l,k)     = coeff(gX(l)-gX(k),     P(1),     P(2), lX, filling)
       Qg(l,k)     = coeff(gX(l)-gX(k),     Q(1),     Q(2), lX, filling)
       Rg(l,k)     = coeff(gX(l)-gX(k),     R(1),     R(2), lX, filling)
       rho11g(l,k) = coeff(gX(l)-gX(k), rho11(1), rho11(2), lX, filling)
       rho12g(l,k) = coeff(gX(l)-gX(k), rho12(1), rho12(2), lX, filling)
       rho22g(l,k) = coeff(gX(l)-gX(k), rho22(1), rho22(2), lX, filling)
       mubg(l,k)   = coeff(gX(l)-gX(k),   mub(1),   mub(2), lX, filling)
    end do
    end do

  !--y=0とし，x方向の波数を増やす（Γ-X）
    kY = 0.0d0
    do iter = 0, NE-1
     kX = dkX*dble(iter)
     do l=1,ngX
     do k=1,ngX
      A(l    ,k    ) = Pg(l,k)*(kX+gX(k))*(kX+gX(l)) &
                     &+mubg(l,k)*(kY)*(kY)
      A(l    ,k+ngX ) = (Pg(l,k)-2d0*mubg(l,k))*(kY)*(kX+gX(l))&
                     &+mubg(l,k)*(kX+gX(k))*(kY)
      A(l    ,k+ngX2) = Qg(l,k)*(kX+gX(k))*(kX+gX(l))
      A(l    ,k+ngX3) = Qg(l,k)*(kY)*(kX+gX(l))

      A(l+ngX ,k    ) = (Pg(l,k)-2d0*mubg(l,k))*(kX+gX(k))*(kY)&
                     &+mubg(l,k)*(kY)*(kX+gX(l))
      A(l+ngX ,k+ngX ) = Pg(l,k)*(kY)*(kY) &
                   &+mubg(l,k)*((kX+gX(k))*(kX+gX(l)))
      A(l+ngX ,k+ngX2) = Qg(l,k)*(kX+gX(k))*(kY)
      A(l+ngX ,k+ngX3) = Qg(l,k)*(kY)*(kY)

      A(l+ngX2,k    ) = Qg(l,k)*(kX+gX(k))*(kX+gX(l))
      A(l+ngX2,k+ngX ) = Qg(l,k)*(kY)*(kX+gX(l))
      A(l+ngX2,k+ngX2) = Rg(l,k)*(kX+gX(k))*(kX+gX(l))
      A(l+ngX2,k+ngX3) = Rg(l,k)*(kY)*(kX+gX(l))

      A(l+ngX3,k    ) = Qg(l,k)*(kX+gX(k))*(kY)
      A(l+ngX3,k+ngX ) = Qg(l,k)*(kY)*(kY)
      A(l+ngX3,k+ngX2) = Rg(l,k)*(kX+gX(k))*(kY)
      A(l+ngX3,k+ngX3) = Rg(l,k)*(kY)*(kY)

      B(l    ,k    ) = rho11g(l,k)
      B(l    ,k+ngX ) = 0
      B(l    ,k+ngX2) = rho12g(l,k)
      B(l    ,k+ngX3) = 0

      B(l+ngX ,k    ) = 0
      B(l+ngX ,k+ngX ) = rho11g(l,k)
      B(l+ngX ,k+ngX2) = 0
      B(l+ngX ,k+ngX3) = rho12g(l,k)

      B(l+ngX2,k    ) = rho12g(l,k)
      B(l+ngX2,k+ngX ) = 0
      B(l+ngX2,k+ngX2) = rho22g(l,k)
      B(l+ngX2,k+ngX3) = 0

      B(l+ngX3,k    ) = 0
      B(l+ngX3,k+ngX ) = rho12g(l,k)
      B(l+ngX3,k+ngX2) = 0
      B(l+ngX3,k+ngX3) = rho22g(l,k)
     end do
     end do
    call ZGGEV('N', 'V', ngX4, A, ngX4, B, ngX4, ALPHA, BETA, VL, ngX4, VR, ngX4,&
      & WORK, LWORK, RWORK, INFO)
    eigen = SQRT(ALPHA/BETA)*lX/csl
    X=kX*lX/pi
    call DEIGSRT(eigen,VR,ngX4,ngX4)
    write(10,'(1500(e24.10e3,2x),i5)') X, dble(eigen)
    write(20,'(1500(e24.10e3,2x),i5)') X, imag(eigen)
   end do

  ! !--kX=MAXとし，y方向の波数を増やす（X-M）
    kX = kXmax
    do iter = 0, NE-1
     kY = dkY*dble(iter)
      do l=1,ngX
      do k=1,ngX
        A(l    ,k    ) = Pg(l,k)*(kX+gX(k))*(kX+gX(l)) &
                       &+mubg(l,k)*(kY)*(kY)
        A(l    ,k+ngX ) = (Pg(l,k)-2d0*mubg(l,k))*(kY)*(kX+gX(l))&
                       &+mubg(l,k)*(kX+gX(k))*(kY)
        A(l    ,k+ngX2) = Qg(l,k)*(kX+gX(k))*(kX+gX(l))
        A(l    ,k+ngX3) = Qg(l,k)*(kY)*(kX+gX(l))

        A(l+ngX ,k    ) = (Pg(l,k)-2d0*mubg(l,k))*(kX+gX(k))*(kY)&
                       &+mubg(l,k)*(kY)*(kX+gX(l))
        A(l+ngX ,k+ngX ) = Pg(l,k)*(kY)*(kY) &
                     &+mubg(l,k)*((kX+gX(k))*(kX+gX(l)))
        A(l+ngX ,k+ngX2) = Qg(l,k)*(kX+gX(k))*(kY)
        A(l+ngX ,k+ngX3) = Qg(l,k)*(kY)*(kY)

        A(l+ngX2,k    ) = Qg(l,k)*(kX+gX(k))*(kX+gX(l))
        A(l+ngX2,k+ngX ) = Qg(l,k)*(kY)*(kX+gX(l))
        A(l+ngX2,k+ngX2) = Rg(l,k)*(kX+gX(k))*(kX+gX(l))
        A(l+ngX2,k+ngX3) = Rg(l,k)*(kY)*(kX+gX(l))

        A(l+ngX3,k    ) = Qg(l,k)*(kX+gX(k))*(kY)
        A(l+ngX3,k+ngX ) = Qg(l,k)*(kY)*(kY)
        A(l+ngX3,k+ngX2) = Rg(l,k)*(kX+gX(k))*(kY)
        A(l+ngX3,k+ngX3) = Rg(l,k)*(kY)*(kY)

        B(l    ,k    ) = rho11g(l,k)
        B(l    ,k+ngX ) = 0
        B(l    ,k+ngX2) = rho12g(l,k)
        B(l    ,k+ngX3) = 0

        B(l+ngX ,k    ) = 0
        B(l+ngX ,k+ngX ) = rho11g(l,k)
        B(l+ngX ,k+ngX2) = 0
        B(l+ngX ,k+ngX3) = rho12g(l,k)

        B(l+ngX2,k    ) = rho12g(l,k)
        B(l+ngX2,k+ngX ) = 0
        B(l+ngX2,k+ngX2) = rho22g(l,k)
        B(l+ngX2,k+ngX3) = 0

        B(l+ngX3,k    ) = 0
        B(l+ngX3,k+ngX ) = rho12g(l,k)
        B(l+ngX3,k+ngX2) = 0
        B(l+ngX3,k+ngX3) = rho22g(l,k)
      end do
      end do
      call ZGGEV('N', 'V', ngX4, A, ngX4, B, ngX4, ALPHA, BETA, VL, ngX4, VR, ngX4,&
       & WORK, LWORK, RWORK, INFO)
    eigen = SQRT(ALPHA/BETA)*lX/csl
    X=kXmax*lX/pi
    Y=kY*lY/pi
    call DEIGSRT(eigen,VR,ngX4,ngX4)
    write(10,'(1500(e24.10e3,2x),i5)') X+Y, dble(eigen)
    write(20,'(1500(e24.10e3,2x),i5)') X+Y, imag(eigen)
    end do

  ! !--kX=Max,kY=Maxから，原点方向に戻る（M-Γ）
  do iter = 0, NE
   kX = kXmax-dkX*dble(iter)
   kY = kYmax-dkY*dble(iter)
    do l=1,ngX
    do k=1,ngX
      A(l    ,k    ) = Pg(l,k)*(kX+gX(k))*(kX+gX(l)) &
                     &+mubg(l,k)*(kY)*(kY)
      A(l    ,k+ngX ) = (Pg(l,k)-2d0*mubg(l,k))*(kY)*(kX+gX(l))&
                     &+mubg(l,k)*(kX+gX(k))*(kY)
      A(l    ,k+ngX2) = Qg(l,k)*(kX+gX(k))*(kX+gX(l))
      A(l    ,k+ngX3) = Qg(l,k)*(kY)*(kX+gX(l))

      A(l+ngX ,k    ) = (Pg(l,k)-2d0*mubg(l,k))*(kX+gX(k))*(kY)&
                     &+mubg(l,k)*(kY)*(kX+gX(l))
      A(l+ngX ,k+ngX ) = Pg(l,k)*(kY)*(kY) &
                   &+mubg(l,k)*((kX+gX(k))*(kX+gX(l)))
      A(l+ngX ,k+ngX2) = Qg(l,k)*(kX+gX(k))*(kY)
      A(l+ngX ,k+ngX3) = Qg(l,k)*(kY)*(kY)

      A(l+ngX2,k    ) = Qg(l,k)*(kX+gX(k))*(kX+gX(l))
      A(l+ngX2,k+ngX ) = Qg(l,k)*(kY)*(kX+gX(l))
      A(l+ngX2,k+ngX2) = Rg(l,k)*(kX+gX(k))*(kX+gX(l))
      A(l+ngX2,k+ngX3) = Rg(l,k)*(kY)*(kX+gX(l))

      A(l+ngX3,k    ) = Qg(l,k)*(kX+gX(k))*(kY)
      A(l+ngX3,k+ngX ) = Qg(l,k)*(kY)*(kY)
      A(l+ngX3,k+ngX2) = Rg(l,k)*(kX+gX(k))*(kY)
      A(l+ngX3,k+ngX3) = Rg(l,k)*(kY)*(kY)

      B(l    ,k    ) = rho11g(l,k)
      B(l    ,k+ngX ) = 0
      B(l    ,k+ngX2) = rho12g(l,k)
      B(l    ,k+ngX3) = 0

      B(l+ngX ,k    ) = 0
      B(l+ngX ,k+ngX ) = rho11g(l,k)
      B(l+ngX ,k+ngX2) = 0
      B(l+ngX ,k+ngX3) = rho12g(l,k)

      B(l+ngX2,k    ) = rho12g(l,k)
      B(l+ngX2,k+ngX ) = 0
      B(l+ngX2,k+ngX2) = rho22g(l,k)
      B(l+ngX2,k+ngX3) = 0

      B(l+ngX3,k    ) = 0
      B(l+ngX3,k+ngX ) = rho12g(l,k)
      B(l+ngX3,k+ngX2) = 0
      B(l+ngX3,k+ngX3) = rho22g(l,k)
    end do
    end do
    call ZGGEV('N', 'V', ngX4, A, ngX4, B, ngX4, ALPHA, BETA, VL, ngX4, VR, ngX4,&
     & WORK, LWORK, RWORK, INFO)
  eigen = SQRT(ALPHA/BETA)*lX/csl
  X=kXmax*lX/pi
  Y=kY*lY/pi
  call DEIGSRT(eigen,VR,ngX4,ngX4)
  write(10,'(1500(e24.10e3,2x),i5)') X+Y, dble(eigen)
  write(20,'(1500(e24.10e3,2x),i5)') X+Y, imag(eigen)
  end do

  close(10)
  close(20)
!===
end do

end program poroerastic1D_imp

!********************************************************************
!********************************************************************
!==================================================================
! +-------------------------------+
! | フーリエ係数の作成　　　　　　　 |
! +-------------------------------+
!===
function coeff(gX, ap, ag, a, f)
  implicit none
    integer :: k
    integer,parameter::DD=kind(0d0)
    real(DD) :: ap, ag
    real(DD) :: gX !y方向はバルクなのでgYはない
    real(DD) :: a, pi, pi2, f, dA, ph

    complex :: coeff
    complex(DD), parameter :: ai=dcmplx(0.0d0,1.0d0)

    dA=f*a    ! A層の厚さ、lXは周期長、fは充填率
    ph=0.5d0*gX*dA

  !--{a_G} = -i * {ap-ag}/{G*Lx} * (exp[-iG/2Lx] - 1) : G/=0
  !--      = f*{ap}+(1-f)*{ag} : G=0
    if (gX .ne. 0.0d0) then
!      coeff = -(ap-ag)* f * (exp(-ai*gX*lX) - 1.0d0) / (gX * lX)
      coeff = (ap-ag)*f*exp(-ai*ph)*sin(ph)/ph
    else
      coeff = f*ap+(1.0d0-f)*ag
    endif
  return
end function coeff
!===

!***********************************************************************
!.! ROUTINE: EIGSRT
!.! PURPOSE: This routine given the eigenvalues and eigenvectors
!.!          sorts the eigenvalues into ascending order, and rearranges the
!.!          columns of square matrix correspondingly.  The method is straight
!.!          insertion.
!.!
!.!          see pg. 348 w/ explanation pgs.335-376
!.!          Numerical Recipes: The Art of Scientific Programming
!.!          (FORTRAN version), 1st edition
!.!          W.H. Press, B.P. Flannery, S.A. Teukolsky, W.T. Vetterling
!.!          Cambridge Univ. Press., 1986
!.!
!***********************************************************************
SUBROUTINE DEIGSRT(D,V,N,NP)
implicit none
integer*4 :: n, np
complex*16 :: d(np),v(np,np)

integer*4 :: i,j,k
complex*16 :: p

do i=1,n-1
  k = i
  p = d(i)
  do j=i+1,n
    if(dble(d(j)).le.dble(p)) then
      k = j
      p = d(j)
    end if
  enddo
  if( k .ne. i ) then
    d(k) = d(i)
    d(i) = p
    do j = 1,n
      p = v(j,i)
      v(j,i) = v(j,k)
      v(j,k) = p
    enddo
  end if
enddo
return
end
