!C
!C   program for 1D phononic ctystals with poroerasticity
!C   solved bY PME(plain Expansion Method)
!C
!C   格子点を作成
!C   逆格子空間を作成
!C   入射波とパラメタを展開
!C   展開した係数から行列を作成
!C   LAPACKで固有値，固有ベクトルを出す
!C   与えた波数に対応した固有値，固有ベクトルを出力
!C
program Phononic1D
  implicit none
    integer :: iX, ig, l, k, iter, INFO
    integer, parameter :: DD=kind(0d0)
    integer, parameter :: nX=15, NE=100 !C 要素数，刻み数
    integer, parameter :: ngX=2*nX+1, ngx8=ngx*8, LWORK=ngX*4

    real(DD) :: dkX, kX, kXmin, kXmax
    real(DD) :: gX(ngX)
    real(DD) :: rho11p, rho12p, rho22p, Pp, Qp, Rp
    real(DD) :: rhoS, rhoF, Ks, Kf, Kb, mus, mub, f, tor
    real(DD) :: rhog, Kg, mug, Pg     !C 基盤の密度，体積弾性率，せん断弾性率
    real(DD) :: cst, csl, cf, eta, filling, bunbo
    real(DD) :: lX, lY, pi, pi2, RWORK(ngx8*2)

    complex :: coeff
    complex(DD),dimension(1:ngX, 1:ngX) :: P, Q, R, rho11, rho12, rho22
    complex(DD),dimension(1:ngX*2, 1:ngX*2) :: A, B, VL, VR
    complex(DD) :: eigen(ngX*2), ALPHA(ngX*2), BETA(ngX*2), WORK(LWORK)
    complex(DD), parameter :: COM=dcmplx(0.0d0,1.0d0)

    external coeff

!C==================================================================
!C +-------+
!C | INIT. |
!C +-------+
!C===
  ! tor=0.0d0
    open  (11, file='input1D.dat', status='unknown')
      read (11,*) lX, lY              !C 格子定数[m]
      read (11,*) rhoS, rhoF          !C 密度(固体，流体)[kg m-3]
      read (11,*) f, tor              !C 孔隙率，迷路度
      read (11,*) csl, cst, cf        !C 伝搬速度（固体縦波，固体横波，流体）[m s-1]
      read (11,*) Kb, mub             !C バルクの弾性定数（体積弾性率，せん断弾性率）[Pa]
      read (11,*) eta                 !C 流体の粘性係数[Pa s]
      read (11,*) rhog, Kg, mug       !C 基盤の密度[kg m-3]，体積弾性率，せん断弾性率[Pa]
      read (11,*) filling     ! 充填率　
                              ! dA=filling*lX       :A層の厚さ
                              ! dB=(1-filling)*lX   :B層の厚さ
    close (11)

    pi=4.0d0*datan(1.0d0)
    pi2=2.0d0*pi
    ! if ( tor .eq. 0.0d0) then
    ! tor=f**(-2.0d0/3.0d0)                     !C 迷路度
    ! end if


!C===

!C==================================================================
!C +---------------------------------+
!C | 弾性定数を，実験で得た速度から計算 |
!C +---------------------------------+
!C===
  !C--{csl}=[({Ks}+4{mus}/3)/{rhoS}]**(1/2)
  !C--{cst}=({mus}/{rhoS})**(1/2)
  !C--{cf}=({Kf}/{rhoF})**(1/2)
    mus = cst**2 * rhoS
    Ks  = csl**2 * rhoS - (4*mus/3)
    Kf  =  cf**2 * rhoF
  !  write(*,*) 'mus=',mus
  !  write(*,*) 'Ks =',Ks
  !  write(*,*) 'Kf =',Kf
  !  write(*,*) 'tor =',tor
!C===

!C==================================================================
!C +--------------------------------+
!C | 逆格子ベクトルG1(G),G2(G")の作成 |
!C +--------------------------------+
!C===
  !C--{第一ブリルアンゾーンの境界} = {逆格子ベクトルの辺垂直二等分線}
  !C--{刻み幅} = {格子の辺長}/{刻み数}
    kXmin = 0.0d0
    kXmax = pi/lX
    dkX = (kXmax - kXmin) / NE

  !C--{点の数} = 2*{正方向の要素数}+{原点}
  !C--{点の総数} = {X方向の点の数}*{Y方向の点の数}

  !C--{位相空間のベクトルGi(ITER)}= 2*{pi}/{辺長Li} * ITER
    iter=0
    do iX = -nX, nX
     iter = iter+1
     gX(iter) = (pi2/lX)*dble(iX)
    end do

!C==================================================================
!C +--------------------------------------------+
!C | 多孔質体中の物理パラメータを定める  　　　　　 |
!C +--------------------------------------------+
!C===
  !C--多孔質体中では
    rho11p = (1.0d0-f)*rhoS + (tor-1.0d0)*f*rhoF
    rho12p = -(tor-1.0d0)*f*rhoF
    rho22p = tor*f*rhoF
    bunbo=1.0d0-f-Kb/Ks + f*Ks/Kf
    Pp = Ks*( (1.0d0-f)*(1.0d0-f-Kb/Ks) + f*Kb/Kf )/bunbo + 4.0d0*mub/3.0d0
    Qp = f*Ks*(1.0d0-f-Kb/Ks)/bunbo
    Rp = (f**2 *Ks)/bunbo
!    write(*,*) 'Pp =',Pp
!    write(*,*) 'Qp =',Qp
!    write(*,*) 'Rp =',Rp

  !C--基盤中では
  !C--{Pg} = {Kg} + {4/3 mug}
  !C--{Qg} = {Rg} = 0
  !C--{rho11g} = {rhog}
  !C--{rho12g} = {rho22g} = 0
    Pg = (Kg) + (4*mug/3)
!    write(*,*) 'Pg =',Pg

!C==================================================================
!C +---------------------------------------------+
!C | 固有方程式のための行列A,Bを作成　　　　　　　　 |
!C | {固有値}*{B}.{固有ベクトル}={A}.{固有ベクトル} |
!C | ZGGEVを使って固有値求める                     |
!C +---------------------------------------------+
!C===
  !C--y=0とし，x方向の波数を増やす．
  open(10,file='1D_ZGGEV.dat')
  open(20,file='1D_ZGGEV2.dat')
    do l=1,ngX
    do k=1,ngX
       P(l,k)     = coeff(gX(l)-gX(k), Pp,       Pg, lx, filling)
       Q(l,k)     = coeff(gX(l)-gX(k), Qp,        0, lx, filling)
       R(l,k)     = coeff(gX(l)-gX(k), Rp,        0, lx, filling)
       rho11(l,k) = coeff(gX(l)-gX(k), rho11p, rhog, lx, filling)
       rho12(l,k) = coeff(gX(l)-gX(k), rho12p,    0, lx, filling)
       rho22(l,k) = coeff(gX(l)-gX(k), rho22p,    0, lx, filling)
    end do
    end do
  do iter = -NE, NE
     kX = dkX*dble(iter)
     do l=1,ngX
     do k=1,ngX
        A(l    ,k    ) = P(l,k)*(kX+gX(k))*(kX+gX(l))
        A(l    ,k+ngX) = Q(l,k)*(kX+gX(k))*(kX+gX(l))
        A(l+ngX,k    ) = Q(l,k)*(kX+gX(k))*(kX+gX(l))
        A(l+ngX,k+ngX) = R(l,k)*(kX+gX(k))*(kX+gX(l))
        B(l    ,k    ) = rho11(l,k)
        B(l    ,k+ngX) = rho12(l,k)
        B(l+ngX,k    ) = rho12(l,k)
        B(l+ngX,k+ngX) = rho22(l,k)
     end do
     end do

     call ZGGEV('N', 'V', ngX*2, A, ngX*2, B, ngX*2, ALPHA, BETA, VL, ngX*2, VR, ngX*2,&
      & WORK, LWORK, RWORK, INFO)
     eigen = SQRT(ALPHA/BETA)/pi2
     write(10,'(500(e24.10e3,2x),i5)') kX, dble(eigen)
     write(20,'(500(e24.10e3,2x),i5)') kX, imag(eigen)
  end do

end program Phononic1D
!C********************************************************************
!C********************************************************************
!C********************************************************************
!C********************************************************************
!C==================================================================
!C +-------------------------------+
!C | フーリエ係数の作成　　　　　　　 |
!C +-------------------------------+
!C===
function coeff(gX, ap, ag, lx, f)
  implicit none
    integer :: k
    integer,parameter::DD=kind(0d0)

    real(DD) :: ap, ag
    real(DD) :: gX
    real(DD) :: lX, pi, pi2, f, dA, ph

    complex :: coeff
    complex(DD), parameter :: ai=dcmplx(0.0d0,1.0d0)

    dA=f*lX    ! A層の厚さ、lXは周期長、fは充填率
    ph=0.5d0*gX*dA

  !C--{a_G} = -i * {ap-ag}/{G*Lx} * (exp[-iG/2Lx] - 1) : G/=0
  !C--      = {ap-ag}/2 : G=0
    if (gX .ne. 0.0d0) then
!      coeff = -(ap-ag)* f * (exp(-ai*gX*lX) - 1.0d0) / (gX * lX)
      coeff = (ap-ag)*f*exp(-ai*ph)*sin(ph)/ph
    else
      coeff = f*ap+(1.0d0-f)*ag
    endif
  return
end function coeff
