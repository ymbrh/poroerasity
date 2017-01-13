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
    integer :: NE, NE2
    integer, parameter :: DD=kind(0d0)
    integer, parameter :: nX=7
    integer, parameter :: ngX=2*nX+1, ngx8=ngx*8, LWORK=2*2

    real(DD) :: dkX, kX, kXmin, kXmax
    real(DD) :: gX1(ngX), gX2(ngX)
    real(DD) :: rho11p, rho12p, rho22p, Pp, Qp, Rp
    real(DD) :: rhoS, rhoF, Ks, Kf, Kb, mus, mub, f, tor
    real(DD) :: rhog, Kg, mug, Pg     !C 基盤の密度，体積弾性率，せん断弾性率
    real(DD) :: cst, csl, cf, eta
    real(DD) :: lX, lY, pi, pi2, RWORK(16)

    complex :: coeff
    complex(DD),dimension(1:ngX, 1:ngX) :: P, Q, R, rho11, rho12, rho22
    complex(DD),dimension(1:ngX, 1:ngX) :: A, B, VL, VR
    complex(DD) :: eigen(2), ALPHA(2), BETA(2), WORK(500)
    complex(DD), parameter :: COM=dcmplX(0.0d0,1.0d0)

    external coeff

!C==================================================================
!C +-------+
!C | INIT. |
!C +-------+
!C===
    open  (11, file='input.dat', status='unknown')
      read (11,*) NE                  !C 刻み数
      read (11,*) lX, lY              !C 格子定数[m]
      read (11,*) rhoS, rhoF          !C 密度(固体，流体)[kg m-3]
      read (11,*) f, tor              !C 孔隙率，迷路度
      read (11,*) csl, cst, cf        !C 伝搬速度（固体縦波，固体横波，流体）[m s-1]
      read (11,*) Kb, mub             !C バルクの弾性定数（体積弾性率，せん断弾性率）[Pa]
      read (11,*) eta                 !C 流体の粘性係数[Pa s]
      read (11,*) rhog, Kg, mug       !C 基盤の密度[kg m-3]，体積弾性率，せん断弾性率[Pa]
    close (11)

    pi=4.0d0*atan(1.0d0)
    pi2=2.0d0*pi
    NE2=NE*2
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
    Ks = csl**2 * rhoS - (4*mus/3)
    Kf = cf**2 * rhoF
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
    dkX = (kXmax - kXmin) / NE2

  !C--{点の数} = 2*{正方向の要素数}+{原点}
  !C--{点の総数} = {X方向の点の数}*{Y方向の点の数}

  !C--{位相空間のベクトルGi(ITER)}= 2*{pi}/{辺長Li} * ITER
    iter=0
    do iX = -nX, nX
     iter = iter+1
     gX1(iter) = (pi2/lX)*dble(iX)
     gX2(iter) = (pi2/lX)*dble(iX)
    end do

!C==================================================================
!C +--------------------------------------------+
!C | 多孔質体中の物理パラメータを定める  　　　　　 |
!C +--------------------------------------------+
!C===
  !C--多孔質体中では
    rho11p = (1-f)*rhoS + (tor-1)*f*rhoF
    rho12p = -(tor-1)*f*rhoF
    rho22p = tor*f*rhoF
    Pp = Ks*( (1-f)/(1-f-Kb/Ks) + (f*Kb/Kf) )/(1-f-Kb/Ks + f*Ks/Kf) + 4*mub/3
    Qp = f*Ks*(1-f-Kb/Ks)/(1-f-Kb/Ks + f*Ks/Kf)
    Rp = (f**2 *Ks)/(1-f-Kb/Ks + f*Ks/Kf)

  !C--基盤中では
  !C--{Pg} = {Kg} + {4/3 mug}
  !C--{Qg} = {Rg} = 0
  !C--{rho11g} = {rhog}
  !C--{rho12g} = {rho22g} = 0
    Pg = (Kg) + (4*mug/3)

!C==================================================================
!C +---------------------------------------------+
!C | 固有方程式のための行列A,Bを作成　　　　　　　　 |
!C | {固有値}*{B}.{固有ベクトル}={A}.{固有ベクトル} |
!C | ZGGEVを使って固有値求める                     |
!C +---------------------------------------------+
!C===
  !C--y=0とし，x方向の波数を増やす．
  open(10,file='1D_ZGGEV.dat')
    do iter = 0, NE2
     kX = dkX*dble(iter)

     A(1,1)=0
     A(1,2)=0
     A(1,3)=0
     A(1,4)=0
     B(1,1)=0
     B(1,2)=0
     B(1,3)=0
     B(1,4)=0
     do l=1,ngX
      do k=1,ngX
        P(l,k) = coeff(gX2(l)-gX1(k), Pp, Pg, lx)
        Q(l,k) = coeff(gX2(l)-gX1(k), Qp, 0, lx)
        R(l,k) = coeff(gX2(l)-gX1(k), Rp, 0, lx)
        rho11(l,k) = coeff(gX2(l)-gX1(k), rho11p, rhog, lx)
        rho12(l,k) = coeff(gX2(l)-gX1(k), rho12p, 0, lx)
        rho22(l,k) = coeff(gX2(l)-gX1(k), rho22p, 0, lx)

        A(1,1) = A(1,1) + P(l,k)*(kX+gX1(k))*(kX+gX2(l))
        A(1,2) = A(1,2) + Q(l,k)*(kX+gX1(k))*(kX+gX2(l))
        A(2,1) = A(2,1) + Q(l,k)*(kX+gX1(k))*(kX+gX2(l))
        A(2,2) = A(2,2) + R(l,k)*(kX+gX1(k))*(kX+gX2(l))

        B(1,1) = B(1,1) + rho11(l,k)
        B(1,2) = B(1,2) + rho12(l,k)
        B(2,1) = B(2,1) + rho12(l,k)
        B(2,2) = B(2,2) + rho22(l,k)
      end do
     end do

     call ZGGEV('N', 'V', 2, A, 2, B, 2, ALPHA, BETA, VL, 2, VR, 2,&
      & WORK, LWORK, RWORK, INFO)
     eigen = SQRT(ALPHA/BETA)
     write(10,'(25(e24.10e3,2x),i5)') kX, eigen ,ALPHA,BETA
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
function coeff(gX, ap, ag, lx)
  implicit none
    integer :: k
    integer,parameter::DD=kind(0d0)

    real(DD) :: ap, ag
    real(DD) :: kX, gX
    real(DD) :: lX, lX2, pi, pi2

    complex :: coeff
    complex(DD), parameter :: COM=dcmplX(0.0d0,1.0d0)

    lX2=lX*2
  !C--{a_G} = -i * {ap-ag}/{G*Lx} * (exp[-iG/2Lx] - 1) : G/=0
  !C--      = {ap-ag}/2 : G=0
    if (gX .ne. 0) then
      coeff = -(ap-ag) * (exp(-COM*gX / lX2) - 1) / (gX * lX)
    else
      coeff = (ap-ag)/2
    endif
  return
end function coeff
