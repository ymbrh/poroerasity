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
program Phononic1D
  implicit none

    character :: material*4, nr*3

    integer :: iX, l, k, iter, INFO, i
    integer, parameter :: DD=kind(0d0)
    integer, parameter :: nX=150, NE=100 ! 要素数，刻み数
    integer, parameter :: ngX=2*nX+1, ngx8=ngx*8, LWORK=ngX*4

    real(DD) :: dkX, kX, kXmin, kXmax, X
    real(DD) :: gX(ngX)
    real(DD),dimension(1:2) :: rho11, rho12, rho22, P, Q, R
    real(DD) :: rhoS, rhoF, Ks, Kf, Kb, mub, f, tor, filling, bunbo
    real(DD) :: rhoSp, rhoFp, Ksp, Kfp, Kbp, musp, mubp, fp
    real(DD) :: rhoSh, rhoFh, Ksh, Kfh, Kbh, mush, mubh, fh
    real(DD) :: lX, pi, pi2, RWORK(ngx8*2)
    !--多孔質体(シリカ)の波の伝搬速度
    real(DD), parameter :: csl=5.970d3, cst=3.760d3

    complex :: coeff
    complex(DD),dimension(1:ngX, 1:ngX) :: Pg, Qg, Rg, rho11g, rho12g, rho22g
    complex(DD),dimension(1:ngX*2, 1:ngX*2) :: A, B, VL, VR
    complex(DD) :: eigen(ngX*2), ALPHA(ngX*2), BETA(ngX*2), WORK(LWORK)
    complex(DD), parameter :: COM=dcmplx(0.0d0,1.0d0)

    external coeff, DEIGSRT

!--{input_xxxについて連続で計算させるルーチン}
do 100 i=1,8
write(nr,'(I3.3)') i
!==================================================================
! +-------+
! | INIT. |
! +-------+
!===
    open  (11, file='input1D_'//nr//'.dat', status='unknown')
    !--フォノニック結晶の構造
      read (11,*) lX                  ! 格子定数[m]
      read (11,*) filling             ! 充填率　
                                      ! dA=filling*lX       :A層の厚さ
                                      ! dB=(1-filling)*lX   :B層の厚さ
    !--多孔質体(poroelastic)の物質パラメタ
      read (11,*) rhoSp, rhoFp        ! 多孔質体での密度(固体，流体)[kg m-3]
      read (11,*) Ksp, Kfp            ! 体積弾性率（固体，流体）[Pa]
      read (11,*) Kbp, mubp           ! バルクの弾性定数（体積弾性率，せん断弾性率）[Pa]
      read (11,*) fp                  ! 孔隙率
    !--基盤(Host)の物質パラメタ
      read (11,*) rhoSh, rhoFh        ! ホストでの密度(固体，流体)[kg m-3]
      read (11,*) Ksh, Kfh            ! 体積弾性率（固体，流体）[Pa]
      read (11,*) Kbh, mubh           ! バルクの弾性定数（体積弾性率，せん断弾性率）[Pa]
      read (11,*) fh                  ! 孔隙率
    close (11)

     pi=4.0d0*datan(1.0d0)
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

  !--{逆格子点の数} = 2*{正方向の要素数}+{原点}
  !--{位相空間のベクトルGi(r)}= 2*{pi}/{辺長Li} * r
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
  open (25, file='elastic1D_'//nr//'.dat', status='unknown')
  do iter = 1,2
  !--多孔質体中では
    if ( iter .eq. 1 ) then
      f = fp
      rhoS = rhoSp
      rhoF = rhoFp
      Ks = Ksp
      Kf = Kfp
      Kb = Kbp
      mub = mubp
      material = "poro"
      tor=f**(-2.0d0/3.0d0)                     ! 迷路度
  !--基盤中では
    else
      f = fh
      rhoS = rhoSh
      rhoF = rhoFh
      Ks = Ksh
      Kf = Kfh
      Kb = Kbh
      mub = mubh
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
!********************************************************************
    end if
    rho11(iter) = (1.0d0-f)*rhoS + (tor-1.0d0)*f*rhoF
    rho12(iter) = -(tor-1.0d0)*f*rhoF
    rho22(iter) = tor*f*rhoF
    bunbo=1.0d0-f-Kb/Ks + f*Ks/Kf
    P(iter) = Ks*( (1.0d0-f)*(1.0d0-f-Kb/Ks) + f*Kb/Kf )/bunbo + 4.0d0*mub/3.0d0
    Q(iter) = f*Ks*(1.0d0-f-Kb/Ks)/bunbo
    R(iter) = (f**2 *Ks)/bunbo
    write(25,*) material
    write(25,*) 'rhoS =', rhoS
    write(25,*) 'rhoF =', rhoF
    write(25,*) 'Ks =', Ks
    write(25,*) 'Kf =', Kf
    write(25,*) 'Kb =', Kb
    write(25,*) 'mub =', mub
    write(25,*) 'f =', f
    write(25,*) 'tor =',tor
  end do
    write(25,*) 'P =',P
    write(25,*) 'Q =',Q
    write(25,*) 'R =',R
    write(25,*) 'rho11 =', rho11
    write(25,*) 'rho12 =', rho12
    write(25,*) 'rho22 =', rho22
  close(25)


!==================================================================
! +---------------------------------------------+
! | 固有方程式のための行列A,Bを作成　　　　　　　　 |
! | {固有値}*{B}.{固有ベクトル}={A}.{固有ベクトル} |
! | ZGGEVを使って固有値求める                     |
! +---------------------------------------------+
!===
  !--y=0とし，x方向の波数を増やす．
  open(10,file='1D_Re_'//nr//'.dat')
  open(20,file='1D_Im_'//nr//'.dat')
  ! open(26,file='1D_EV_'//nr//'.dat')
    do l=1,ngX
    do k=1,ngX
       Pg(l,k)     = coeff(gX(l)-gX(k),     P(1),     P(2), lx, filling)
       Qg(l,k)     = coeff(gX(l)-gX(k),     Q(1),     Q(2), lx, filling)
       Rg(l,k)     = coeff(gX(l)-gX(k),     R(1),     R(2), lx, filling)
       rho11g(l,k) = coeff(gX(l)-gX(k), rho11(1), rho11(2), lx, filling)
       rho12g(l,k) = coeff(gX(l)-gX(k), rho12(1), rho12(2), lx, filling)
       rho22g(l,k) = coeff(gX(l)-gX(k), rho22(1), rho22(2), lx, filling)
    end do
    end do
  do iter = -NE, NE
     kX = dkX*dble(iter)
     do l=1,ngX
     do k=1,ngX
        A(l    ,k    ) = Pg(l,k)*(kX+gX(k))*(kX+gX(l))
        A(l    ,k+ngX) = Qg(l,k)*(kX+gX(k))*(kX+gX(l))
        A(l+ngX,k    ) = Qg(l,k)*(kX+gX(k))*(kX+gX(l))
        A(l+ngX,k+ngX) = Rg(l,k)*(kX+gX(k))*(kX+gX(l))
        B(l    ,k    ) = rho11g(l,k)
        B(l    ,k+ngX) = rho12g(l,k)
        B(l+ngX,k    ) = rho12g(l,k)
        B(l+ngX,k+ngX) = rho22g(l,k)
     end do
     end do

    call ZGGEV('N', 'V', ngX*2, A, ngX*2, B, ngX*2, ALPHA, BETA, VL, ngX*2, VR, ngX*2,&
      & WORK, LWORK, RWORK, INFO)
    eigen = sqrt(ALPHA/BETA)*lX/csl
    X=kX*lX/pi
    call DEIGSRT(eigen,VR,ngX*2,ngX*2)
    write(10,'(500(e24.10e3,2x),i5)') X, dble(eigen)
    write(20,'(500(e24.10e3,2x),i5)') X, imag(eigen)
    ! write(26,'(500(e24.10e3,2x),i5)') X, dble(VR)
  end do
!===
100 continue

end program Phononic1D
!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************
!==================================================================
! +-------------------------------+
! | フーリエ係数の作成　　　　　　　 |
! +-------------------------------+
!===
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
!.! PURPOSE: This routine given the eigenvalues and eigenvectorps
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
