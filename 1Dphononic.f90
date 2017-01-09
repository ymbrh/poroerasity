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

    integer :: iX, iY, ig, l, k, iter, INFO
    integer :: nX, nY, ngX, ngY, ng, NE
    integer :: LWORK=8
    integer, parameter :: DD=kind(0d0) , NE2=NE*2

    real(DD) :: dkX, dkY, kX, kY, kXmin, kYmin, kXmax, kYmax
    real(DD) :: gX1, gY1, gX2, gY2, gr, X, Y
    real(DD) :: rho11, rho12, rho22, P, Q, R,
    real(DD) :: rhoS, rhoF, Ks, Kf, Kb, mus, mub, f, tor, Kg, mug, rhog
    real(DD) :: rhog, Kg, mug, PP     !C 基盤の密度，体積弾性率，せん断弾性率
    real(DD) :: cst, csl, cf, eta
    real(DD) :: lX, lY, pi, pi2, RWORK(32)

    complex(DD) :: A(4,4), B(4,4), VL(4,4), VR(4,4)
    complex(DD) :: eigen(4), ALPHA(4), BETA(4), WORK(500)
    complex(DD), parameter :: ai=dcmplX(0.0d0,1.0d0)


!C==================================================================
!C +-------+
!C | INIT. |
!C +-------+
!C===
    open  (11, file='input.dat', status='unknown')
      read (11,*) nX, nY              !C 要素数
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
!C +-------------+
!C | 刻み幅を設定 |
!C +-------------+
!C===
  !C--{刻み幅} = {格子の辺長}/{刻み数}
    kXmin = 0.0d0
    kXmax = pi/lX
    dkX = (kXmax - kXmin) / NE2
    kYmin = 0.0d0
    kYmax = pi/lY
    dkY = (kYmax - kYmin) / NE2

!C===

!C==================================================================
!C +----------------------------------------------+
!C | バルクのときの物理パラメータを定める  　　　　　 |
!C | あとでフォノニック結晶にする時は，この部分を元に |
!C | フーリエ級数展開するサブルーチンを追加          |
!C +----------------------------------------------+
!C===
  !C--y=0とし，x方向の波数を増やしたときの行を作成．
    rho11 = (1-f)*rhoS + (tor-1)*f*rhoF
    rho12 = -(tor-1)*f*rhoF
    rho22 = tor*f*rhoF

    P = Ks*( (1-f)/(1-f-Kb/Ks) + (f*Kb/Kf) )/(1-f-Kb/Ks + f*Ks/Kf) + 4*mub/3
    Q = f*Ks*(1-f-Kb/Ks)/(1-f-Kb/Ks + f*Ks/Kf)
    R = (f**2 *Ks)/(1-f-Kb/Ks + f*Ks/Kf)
  !C--弾性定数だけおいとく
    PP = (Kg) + (4*mug/3)

!C==================================================================
!C +---------------------------------------------+
!C | 固有方程式のための行列A,Bを作成　　　　　　　　 |
!C | {固有値}*{B}.{固有ベクトル}={A}.{固有ベクトル} |
!C | ZGGEVを使って固有値求める                     |
!C +---------------------------------------------+
!C===
  !C--y=0とし，x方向の波数を増やす．
  open(10,file='unitXY_ZGGEV.dat')
    do iter = 0, NE2
     kX = dkX*dble(iter)
     kY = 0.0d0

  !C--多孔質体の中ではBiotのモデル
      if ( iter .le. NE ) then
        A(1,1) = P * (kX**2) + mub*(kY**2)
        A(1,2) = mub*(kX*kY)
        A(1,3) = Q*(kX**2)
        A(1,4) = 0
        A(2,1) = mub*(kX*kY)
        A(2,2) = P*(kY**2) + mub*(kX**2)
        A(2,3) = 0
        A(2,4) = Q*(kY**2)
        A(3,1) = Q*(kX**2)
        A(3,2) = 0
        A(3,3) = R*(kX**2)
        A(3,4) = 0
        A(4,1) = 0
        A(4,2) = Q*(kY**2)
        A(4,3) = 0
        A(4,4) = R*(kY**2)

        B(1,1) = rho11
        B(1,2) = 0
        B(1,3) = rho12
        B(1,4) = 0
        B(2,1) = 0
        B(2,2) = rho11
        B(2,3) = 0
        B(2,4) = rho12
        B(3,1) = rho12
        B(3,2) = 0
        B(3,3) = rho22
        B(3,4) = 0
        B(4,1) = 0
        B(4,2) = rho12
        B(4,3) = 0
        B(4,4) = rho22

  !C--基盤中では
  !C--{P} = {Kg} + {4/3 mug}
  !C--{Q} = {R} = 0
  !C--{rho11} = {rhog}
  !C--{rho12} = {rho22} = 0
      else
        A(1,1) = PP * (kX**2) + mub*(kY**2)
        A(1,2) = mub*(kX*kY)
        A(1,3) = Q*(kX**2)
        A(1,4) = 0
        A(2,1) = mub*(kX*kY)
        A(2,2) = P*(kY**2) + mub*(kX**2)
        A(2,3) = 0
        A(2,4) = Q*(kY**2)
        A(3,1) = Q*(kX**2)
        A(3,2) = 0
        A(3,3) = R*(kX**2)
        A(3,4) = 0
        A(4,1) = 0
        A(4,2) = Q*(kY**2)
        A(4,3) = 0
        A(4,4) = R*(kY**2)

        B(1,1) = rho11
        B(1,2) = 0
        B(1,3) = rho12
        B(1,4) = 0
        B(2,1) = 0
        B(2,2) = rho11
        B(2,3) = 0
        B(2,4) = rho12
        B(3,1) = rho12
        B(3,2) = 0
        B(3,3) = rho22
        B(3,4) = 0
        B(4,1) = 0
        B(4,2) = rho12
        B(4,3) = 0
        B(4,4) = rho22
      end if

     call ZGGEV('N', 'V', 4, A, 4, B, 4, ALPHA, BETA, VL, 4, VR, 4,&
      & WORK, LWORK, RWORK, INFO)
     eigen = ALPHA/BETA
     write(10,'(25(e24.10e3,2x),i5)') kX, eigen
    end do

end program Phononic1D
!C********************************************************************
!C********************************************************************
!C********************************************************************
!C********************************************************************
