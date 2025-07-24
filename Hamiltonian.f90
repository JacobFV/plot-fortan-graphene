 module Hamiltonian
 
 use lapack_routines
 
 implicit none
 
 integer, parameter :: dp = 8 ! double precision dp=8, single precision dp=4
 integer(dp), parameter :: numk   = 30 ! number of k-points, use multiples of 6
 integer(dp), parameter :: ntheta = 5 ! commensurate twist angle 
 integer(dp), parameter :: ndim   = 4*(ntheta**2+(ntheta+1)**2+(ntheta+1)*ntheta) ! dimension of Ham
 integer(dp), parameter :: numb = 8 ! total bands to be plot (maximal value is ndim)
 integer(dp), parameter :: ncb = 4 ! conduction bands to be plot (maximal value is ndim/2)
 integer(dp), parameter :: nrelax = 1 ! include lattice relaxation by nrelax=1 (0 otherwise)

! parameters for the Hamiltonian
 real(dp), parameter :: fphase=0.0_dp
 real(dp), parameter :: pressure=1._dp ! tune "twist anlge" by hydrostatic pressure, e.g., 1.096_dp
 real(dp), parameter :: deltaR = 0.1_dp
 real(dp), parameter :: a0 = 1.0_dp/sqrt(3.0_dp)
 real(dp), parameter :: r0 = 0.184_dp
 real(dp), parameter :: d0 = 1.35772_dp/pressure
 real(dp), parameter :: tperp = 1._dp
 real(dp), parameter :: mass1=0._dp ! in eV
 real(dp), parameter :: tz = 1.35772_dp/pressure
 real(dp), parameter :: a1(2) = [ 0.5_dp,sqrt(3.0_dp)*0.5_dp]
 real(dp), parameter :: a2(2) = [-0.5_dp,sqrt(3.0_dp)*0.5_dp]

! general parameters
 real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
 real(dp), parameter :: sq3 = sqrt(3.0_dp)
 real(dp), parameter :: vkd=4.0_dp*pi/3.0_dp

contains


 

!c!c!cC
!c!c!cC
!c!c!cC

subroutine WignerSeitz(coord,t1,t2,t3,cs,sn)

 real(dp), intent(in)    :: t1(2), t2(2), t3(2), cs, sn
 real(dp), intent(inout) :: coord(ndim,3)

 integer(dp) :: na1, nb1, na2, nb2, n1, n2, nind, nrad
 real(dp)    :: rMax, rTemp1, rTemp2, rTemp3, an1, an2, bn1, bn2

nrad=3*ntheta

! initialize running parameters

nind=0
na1=0
nb1=0
na2=0
nb2=0

 rMax = 3.0_dp*ntheta**2+3.0_dp*ntheta+1.0_dp

 do n1=-nrad,nrad
   do n2=-nrad,nrad

! Find coordinate of site A layer 1 is within 1st Wigner-Seitz cell
! n1*a1+n2*a2

     rTemp1=abs(n1*(3.0_dp*ntheta+1)+n2*(3.0_dp*ntheta+2))+ 1e-6
     rTemp2=abs(n1-n2*(3.0_dp*ntheta+1))+ 1e-6
     rTemp3=abs(n1*(3.0_dp*ntheta+2)+n2)+ 1e-6

     if(rTemp1 < rMax .and. rTemp2 < rMax .and. rTemp3 < rMax)then
       na1=na1+1
       nind=nind+1
       coord(nind,1)=n1*a1(1)+n2*a2(1)
       coord(nind,2)=n1*a1(2)+n2*a2(2)
       coord(nind,3)=0.0_dp
     end if

   end do
 end do

 do n1=-nrad,nrad
   do n2=-nrad,nrad

! Find coordinate of site B layer 1 is within 1st Wigner-Seitz cell
! n1*a1+n2*a2

     bn1=n1+1.0_dp/3.0_dp
     bn2=n2+1.0_dp/3.0_dp

     rTemp1=abs(bn1*(3.0_dp*ntheta+1)+bn2*(3.0_dp*ntheta+2))+ 1e-6
     rTemp2=abs(bn1-bn2*(3.0_dp*ntheta+1))+ 1e-6
     rTemp3=abs(bn1*(3.0_dp*ntheta+2)+bn2)+ 1e-6

     if(rTemp1 < rMax .and. rTemp2 < rMax .and. rTemp3 < rMax)then
       nb1=nb1+1
       nind=nind+1
       coord(nind,1)=bn1*a1(1)+bn2*a2(1)
       coord(nind,2)=bn1*a1(2)+bn2*a2(2)
       coord(nind,3)=0.0_dp
     end if

   end do
 end do

!!! Include B atom at zone boundary

nind=nind+1

coord(nind,1)=(t2(1)+t3(1))/3.0_dp
coord(nind,2)=(t2(2)+t3(2))/3.0_dp
coord(nind,3)=0.0_dp

 do n1=-nrad,nrad
   do n2=-nrad,nrad

! Find coordinate of site A layer 2 is within 1st Wigner-Seitz cell
! n1*a1+n2*a2

     an1=n1*(cs-sn/sq3)-2*n2*sn/sq3
     an2=n2*(cs+sn/sq3)+2*n1*sn/sq3

     rTemp1=abs(an1*(3.0_dp*ntheta+1)+an2*(3.0_dp*ntheta+2))+ 1e-6
     rTemp2=abs(an1-an2*(3.0_dp*ntheta+1))+ 1e-6
     rTemp3=abs(an1*(3.0_dp*ntheta+2)+an2)+ 1e-6

     if(rTemp1 < rMax .and. rTemp2 .LT. rMax .and. rTemp3 < rMax)then
       na2=na2+1
       nind=nind+1
       coord(nind,1)=an1*a1(1)+an2*a2(1)
       coord(nind,2)=an1*a1(2)+an2*a2(2)
       coord(nind,3)=tz
     end if

   end do
 end do

! Find coordinate of site B layer 2 is within 1st Wigner-Seitz cell
! n1*a1+n2*a2

 do n1=-nrad,nrad
   do n2=-nrad,nrad

     bn1=(n1+1.0_dp/3.0_dp)*(cs-sn/sq3)-2*(n2+1.0_dp/3.0_dp)*sn/sq3
     bn2=(n2+1.0_dp/3.0_dp)*(cs+sn/sq3)+2*(n1+1.0_dp/3.0_dp)*sn/sq3

     rTemp1=abs(bn1*(3.0_dp*ntheta+1)+bn2*(3.0_dp*ntheta+2))+ 1e-6
     rTemp2=abs(bn1-bn2*(3.0_dp*ntheta+1))+ 1e-6
     rTemp3=abs(bn1*(3.0_dp*ntheta+2)+bn2)+ 1e-6

     if(rTemp1 < rMax .AND. rTemp2 < rMax .AND. rTemp3 < rMax)then
       nb2=nb2+1
       nind=nind+1
       coord(nind,1)=bn1*a1(1)+bn2*a2(1)
       coord(nind,2)=bn1*a1(2)+bn2*a2(2)
       coord(nind,3)=tz
     end if

   end do
 end do

!!! Include B atom at zone boundary

 nind=nind+1

 coord(nind,1)=(t1(1)+t2(1))/3.0_dp
 coord(nind,2)=(t1(2)+t2(2))/3.0_dp
 coord(nind,3)=tz

end subroutine WignerSeitz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine WignerSeitz6(coord,t1,t2,t3,cs,sn)


  real(dp), intent(in)    :: t1(2), t2(2), t3(2), cs, sn

  real(dp), intent(inout) :: coord(ndim,3)

 

  integer(dp) :: na1, nb1, na2, nb2, n1, n2, nind, nrad

  real(dp)    :: rMax, rTemp1, rTemp2, rTemp3, an1, an2, bn1, bn2

 

 nrad=3*ntheta

 

 ! initialize running parameters

 

 nind=0

 na1=0

 nb1=0

 na2=0

 nb2=0

 

  rMax = 3.0_dp*ntheta**2+3.0_dp*ntheta+1.0_dp

 

  do n1=-nrad,nrad

    do n2=-nrad,nrad

 

 ! Find coordinate of site A layer 1 is within 1st Wigner-Seitz cell

 ! n1*a1+n2*a2

      an1=(n1+1.0_dp/3.0_dp)

      an2=(n2-2.0_dp/3.0_dp)

 

      rTemp1=abs(an1*(3.0_dp*ntheta+1)+an2*(3.0_dp*ntheta+2))+ 1e-6

      rTemp2=abs(an1-an2*(3.0_dp*ntheta+1))+ 1e-6

      rTemp3=abs(an1*(3.0_dp*ntheta+2)+an2)+ 1e-6

 

      if(rTemp1 < rMax .and. rTemp2 < rMax .and. rTemp3 < rMax)then

        na1=na1+1

        nind=nind+1

        coord(nind,1)=(n1+1.0_dp/3.0_dp)*a1(1)+(n2-2.0_dp/3.0_dp)*a2(1)

        coord(nind,2)=(n1+1.0_dp/3.0_dp)*a1(2)+(n2-2.0_dp/3.0_dp)*a2(2)

        coord(nind,3)=0.0_dp

      end if

 

    end do

  end do

 

  !! Include A atom at zone boundary

 

 nind=nind+1

 

 coord(nind,1)=(-t2(1)-t1(1))/3.0_dp

 coord(nind,2)=(-t2(2)-t1(2))/3.0_dp

 coord(nind,3)=0.0_dp



  do n1=-nrad,nrad

    do n2=-nrad,nrad

 

 ! Find coordinate of site B layer 1 is within 1st Wigner-Seitz cell

 ! n1*a1+n2*a2

 

      bn1=n1+2.0_dp/3.0_dp

      bn2=n2-1.0_dp/3.0_dp

 

      rTemp1=abs(bn1*(3.0_dp*ntheta+1)+bn2*(3.0_dp*ntheta+2))+ 1e-6

      rTemp2=abs(bn1-bn2*(3.0_dp*ntheta+1))+ 1e-6

      rTemp3=abs(bn1*(3.0_dp*ntheta+2)+bn2)+ 1e-6

 

      if(rTemp1 < rMax .and. rTemp2 < rMax .and. rTemp3 < rMax)then

        nb1=nb1+1

        nind=nind+1

        coord(nind,1)=bn1*a1(1)+bn2*a2(1)

        coord(nind,2)=bn1*a1(2)+bn2*a2(2)

        coord(nind,3)=0.0_dp

      end if

 

    end do

  end do

 

 !! Include B atom at zone boundary

 

 nind=nind+1

 

 coord(nind,1)=(t2(1)+t1(1))/3.0_dp

 coord(nind,2)=(t2(2)+t1(2))/3.0_dp

 coord(nind,3)=0.0_dp

 

  do n1=-nrad,nrad

    do n2=-nrad,nrad

 

 ! Find coordinate of site A layer 2 is within 1st Wigner-Seitz cell

 ! n1*a1+n2*a2

 

      an1=(n1+1.0_dp/3.0_dp)*(cs-sn/sq3)-2*(n2-2.0_dp/3.0_dp)*sn/sq3

      an2=(n2-2.0_dp/3.0_dp)*(cs+sn/sq3)+2*(n1+1.0_dp/3.0_dp)*sn/sq3

 

      rTemp1=abs(an1*(3.0_dp*ntheta+1)+an2*(3.0_dp*ntheta+2))+ 1e-6

      rTemp2=abs(an1-an2*(3.0_dp*ntheta+1))+ 1e-6

      rTemp3=abs(an1*(3.0_dp*ntheta+2)+an2)+ 1e-6

 

      if(rTemp1 < rMax .and. rTemp2 .LT. rMax .and. rTemp3 < rMax)then

        na2=na2+1

        nind=nind+1

        coord(nind,1)=an1*a1(1)+an2*a2(1)

        coord(nind,2)=an1*a1(2)+an2*a2(2)

        coord(nind,3)=tz

      end if

 

    end do

  end do



  !! Include A atom at zone boundary

 

 nind=nind+1

 

 coord(nind,1)=(t1(1)+t2(1))/3.0_dp

 coord(nind,2)=(t1(2)+t2(2))/3.0_dp

 coord(nind,3)=tz

 

 ! Find coordinate of site B layer 2 is within 1st Wigner-Seitz cell

 ! n1*a1+n2*a2

 

  do n1=-nrad,nrad

    do n2=-nrad,nrad

 

      bn1=(n1+2.0_dp/3.0_dp)*(cs-sn/sq3)-2*(n2-1.0_dp/3.0_dp)*sn/sq3

      bn2=(n2-1.0_dp/3.0_dp)*(cs+sn/sq3)+2*(n1+2.0_dp/3.0_dp)*sn/sq3

 

      rTemp1=abs(bn1*(3.0_dp*ntheta+1)+bn2*(3.0_dp*ntheta+2))+ 1e-6

      rTemp2=abs(bn1-bn2*(3.0_dp*ntheta+1))+ 1e-6

      rTemp3=abs(bn1*(3.0_dp*ntheta+2)+bn2)+ 1e-6

 

      if(rTemp1 < rMax .AND. rTemp2 < rMax .AND. rTemp3 < rMax)then

        nb2=nb2+1

        nind=nind+1

        coord(nind,1)=bn1*a1(1)+bn2*a2(1)

        coord(nind,2)=bn1*a1(2)+bn2*a2(2)

        coord(nind,3)=tz

      end if

 

    end do

  end do



 !!! Include B atom at zone boundary

 

  nind=nind+1

 

  coord(nind,1)=(-t1(1)-t2(1))/3.0_dp

  coord(nind,2)=(-t1(2)-t2(2))/3.0_dp

  coord(nind,3)=tz

 

end subroutine WignerSeitz6


 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine ComplexH0(zh,coord,ndim,vk,tn)
 
 integer(dp), intent(in)    :: ndim
 real(dp)   , intent(in)    :: tn(2,3), vk(2), coord(ndim,3)
 complex(dp), intent(out) :: zh(ndim,ndim)
 
 integer(dp) :: i, j, nt
 real(dp)    :: rij, coord_diff
 complex(dp) :: zi_vk_tn,zphase

 
zh = cmplx(0.0_dp,0.0_dp,dp)

 do i = 1 , ndim
   do j = 1 , i

     zphase=exp(cmplx(0.0_dp,-dot_product(vk,coord(i,1:2)-coord(j,1:2))*fphase,dp))
     
     rij = norm2(coord(i,:)-coord(j,:))
 
     coord_diff = (coord(i,3)-coord(j,3))**2
  
     if(coord_diff < deltaR .and. rij > deltaR) then
       zh(i,j) = ft(rij)*zphase
     endif
   
     if(coord_diff > deltaR) then
        zh(i,j) = ftperp(rij)*zphase
     endif
   
!!!!!! Coupling to adjacent Wigner-Seitz cell with tn
     do nt = 1 , 3
       rij = norm2(coord(i,:)-[coord(j,1:2)-tn(:,nt),coord(j,3)])
       zi_vk_tn = cmplx(0.0_dp,dot_product(vk,tn(:,nt)),dp)
       if(coord_diff < deltaR)then
         zh(i,j)=zh(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase
       else
         zh(i,j)=zh(i,j)+ftperp(rij)*exp(-zi_vk_tn)*zphase
       endif
!!!!!! Coupling to adjacent Wigner-Seitz cell with -tn
  rij = norm2(coord(i,:)-[coord(j,1:2)+tn(:,nt),coord(j,3)])
  zi_vk_tn = cmplx(0.0_dp,-dot_product(vk,tn(:,nt)),dp)
  if(coord_diff < deltaR)then
    zh(i,j)=zh(i,j)+ft(rij)*exp(-zi_vk_tn)*zphase
  else
    zh(i,j)=zh(i,j)+ftperp(rij)*exp(-zi_vk_tn)*zphase
  endif
enddo

   end do
 end do

 do i=1,ndim/4
 zh(i,i)=zh(i,i)+mass1
 j=i+ndim/4
 zh(j,j)=zh(j,j)-mass1
 enddo

end subroutine ComplexH0

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure function ftperp(r)

 real(dp), intent(in) :: r
 real(dp) :: ftperp

 real(dp) :: pd

 pd=(d0/r)**2

 ftperp = (-2.7_dp*exp((a0-r)/r0)*(1.0_dp-pd)+0.48_dp*exp((d0*pressure-r)/r0)*pd)*tperp

end function ftperp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure function ft(r)

 real(dp), intent(in) :: r
 real(dp) :: ft
 
! if( r <= 1.01_dp*a0.AND.r>deltaR) then
 if( r>deltaR) then
   ft = -2.7_dp*exp((a0-r)/r0)
 else
   ft = 0.0_dp
 end if

end function ft


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine Relaxation(coord,Gn)

 real(dp), intent(inout) :: coord(ndim,3)
 real(dp), intent(in)    :: Gn(8,2)

 integer(dp) :: i, m, n, ncount, nn
 real(dp) :: ux, uy, rx, ry
 real(dp) :: uqx(5,5), uqy(5,5), delta(2)
 real(dp) :: vk(2)
 real(dp), allocatable :: uq(:,:)

 allocate(uq(ndim,2))
 
uqx=0.0_dp
uqy=0.0_dp 

if(ntheta.EQ.2)then
uqx(2,1)=-0.0007312154977677281
uqy(2,1)=0.0004714203951774816
uqx(3,1)=-1.154897989107238e-06
uqy(3,1)=7.447553144970968e-07
uqx(3,2)=-1.586938089071483e-06
uqy(3,2)=0.0_dp
uqx(4,1)=-2.434676944887857e-09
uqy(4,1)=1.570103901518326e-09
uqx(4,2)=-4.341715993276229e-09
uqy(4,2)=1.000979455471334e-09
uqx(4,3)=-4.341715993276229e-09
uqy(4,3)=-1.000979455471334e-09
uqx(5,1)=-5.701274058766884e-12
uqy(5,1)=3.740789872160164e-12
uqx(5,2)=-1.288484816335295e-11
uqy(5,2)=4.352172282562038e-12
uqx(5,3)=-1.480922507260647e-11
uqy(5,3)=0.0_dp
uqx(5,4)=-1.288484816335295e-11
uqy(5,4)=-4.352172282562038e-12
endif
if(ntheta.EQ.3)then
uqx(2,1)=-0.001435896542049262
uqy(2,1)=0.0008971722346957707
uqx(3,1)=-4.431771204029964e-06
uqy(3,1)=2.770018809111695e-06
uqx(3,2)=-5.980337522539451e-06
uqy(3,2)=0.0_dp
uqx(4,1)=-1.828082429017813e-08
uqy(4,1)=1.142651350640817e-08
uqx(4,2)=-3.192212955644209e-08
uqy(4,2)=7.41532195221355e-09
uqx(4,3)=-3.192212955644209e-08
uqy(4,3)=-7.41532195221355e-09
uqx(5,1)=-8.489859929513374e-11
uqy(5,1)=5.318288726255601e-11
uqx(5,2)=-1.864781523972493e-10
uqy(5,2)=6.321663026418908e-11
uqx(5,3)=-2.126102770742755e-10
uqy(5,3)=0.0_dp
uqx(5,4)=-1.864781523972493e-10
uqy(5,4)=-6.321663026418908e-11
endif
if(ntheta.EQ.4)then
uqx(2,1)=-0.002370539233839467
uqy(2,1)=0.001455154743491223
uqx(3,1)=-1.207488833535449e-05
uqy(3,1)=7.415571455420628e-06
uqx(3,2)=-1.606356439211949e-05
uqy(3,2)=0.0_dp
uqx(4,1)=-8.232092034607089e-08
uqy(4,1)=5.05576290960723e-08
uqx(4,2)=-1.412786099990225e-07
uqy(4,2)=3.296261058695678e-08
uqx(4,3)=-1.412786099990225e-07
uqy(4,3)=-3.296261058695678e-08
uqx(5,1)=-6.338434146746167e-10
uqy(5,1)=3.894401693009466e-10
uqx(5,2)=-1.363506478716756e-09
uqy(5,2)=4.635827193166175e-10
uqx(5,3)=-1.545001808355566e-09
uqy(5,3)=0.0_dp
uqx(5,4)=-1.363506478716756e-09
uqy(5,4)=-4.635827193166175e-10
endif
if(ntheta.EQ.5)then
uqx(2,1)=-0.00352744756912305
uqy(2,1)=0.002140730764691146
uqx(3,1)=-2.68010445377339e-05
uqy(3,1)=1.627423402644026e-05
uqx(3,2)=-3.516363993145865e-05
uqy(3,2)=0.0_dp
uqx(4,1)=-2.729706294022432e-07
uqy(4,1)=1.657616285518053e-07
uqx(4,2)=-4.604302998995516e-07
uqy(4,2)=1.078115792758269e-07
uqx(4,3)=-4.604302998995516e-07
uqy(4,3)=-1.078115792758269e-07
uqx(5,1)=-3.147029295689003e-09
uqy(5,1)=1.911153328301163e-09
uqx(5,2)=-6.629538953061466e-09
uqy(5,2)=2.259053162518054e-09
uqx(5,3)=-7.46551582954267e-09
uqy(5,3)=0.0_dp
uqx(5,4)=-6.629538953061466e-09
uqy(5,4)=-2.259053162518054e-09
endif
if(ntheta.EQ.6)then
uqx(2,1)=-0.004895854840040079
uqy(2,1)=0.002947504859599212
uqx(3,1)=-5.185833778140336e-05
uqy(3,1)=3.124215464415819e-05
uqx(3,2)=-6.704921137114243e-05
uqy(3,2)=0.0_dp
uqx(4,1)=-7.375225497266297e-07
uqy(4,1)=4.443469426422899e-07
uqx(4,2)=-1.22125720193429e-06
uqy(4,2)=2.869558314022102e-07
uqx(4,3)=-1.22125720193429e-06
uqy(4,3)=-2.869558314022102e-07
uqx(5,1)=-1.189937627779894e-08
uqy(5,1)=7.168907495586338e-09
uqx(5,2)=-2.451139179182853e-08
uqy(5,2)=8.36909403937007e-09
uqx(5,3)=-2.741742519309818e-08
uqy(5,3)=0.0_dp
uqx(5,4)=-2.451139179182853e-08
uqy(5,4)=-8.36909403937007e-09
endif
if(ntheta.EQ.7)then
uqx(2,1)=-0.006461449348026135
uqy(2,1)=0.00386704590013867
uqx(3,1)=-9.088953324281387e-05
uqy(3,1)=5.443895214647056e-05
uqx(3,2)=-0.000115650920643144
uqy(3,2)=0.0_dp
uqx(4,1)=-1.719155170137875e-06
uqy(4,1)=1.029772845664776e-06
uqx(4,2)=-2.790313685737744e-06
uqy(4,2)=6.580010547012858e-07
uqx(4,3)=-2.790313685737744e-06
uqy(4,3)=-6.580010547012858e-07
uqx(5,1)=-3.697520081119331e-08
uqy(5,1)=2.214633072942774e-08
uqx(5,2)=-7.434248153831823e-08
uqy(5,2)=2.543295223277813e-08
uqx(5,3)=-8.254619477080724e-08
uqy(5,3)=0.0_dp
uqx(5,4)=-7.434248153831823e-08
uqy(5,4)=-2.543295223277813e-08
endif
if(ntheta.EQ.8)then
uqx(2,1)=-0.008206072919252466
uqy(2,1)=0.004888718057559328
uqx(3,1)=-0.0001477464619472436
uqy(3,1)=8.80996557824374e-05
uqx(3,2)=-0.000184730257310846
uqy(3,2)=0.0_dp
uqx(4,1)=-3.581496930090418e-06
uqy(4,1)=2.135790438804529e-06
uqx(4,2)=-5.688326741335424e-06
uqy(4,2)=1.34655485744528e-06
uqx(4,3)=-5.688326741335424e-06
uqy(4,3)=-1.34655485744528e-06
uqx(5,1)=-9.895113648508195e-08
uqy(5,1)=5.900265937736671e-08
uqx(5,2)=-1.938529162557286e-07
uqy(5,2)=6.645103335141035e-08
uqx(5,3)=-2.135258580530001e-07
uqy(5,3)=0.0_dp
uqx(5,4)=-1.938529162557286e-07
uqy(5,4)=-6.645103335141035e-08
endif
if(ntheta.EQ.9)then
uqx(2,1)=-0.01010772196045918
uqy(2,1)=0.005999689628423322
uqx(3,1)=-0.0002262612975656868
uqy(3,1)=0.0001344419648152459
uqx(3,2)=-0.0002775244610317592
uqy(3,2)=0.0_dp
uqx(4,1)=-6.825276745678818e-06
uqy(4,1)=4.055911695604714e-06
uqx(4,2)=-1.059052135909834e-05
uqy(4,2)=2.517294978569754e-06
uqx(4,3)=-1.059052135909834e-05
uqy(4,3)=-2.517294978569754e-06
uqx(5,1)=-2.352002772937919e-07
uqy(5,1)=1.397520149575992e-07
uqx(5,2)=-4.48274468435581e-07
uqy(5,2)=1.539852461003042e-07
uqx(5,3)=-4.895322080001663e-07
uqy(5,3)=0.0_dp
uqx(5,4)=-4.48274468435581e-07
uqy(5,4)=-1.539852461003042e-07
endif
if(ntheta.EQ.10)then
uqx(2,1)=-0.01214094723296179
uqy(2,1)=0.007185172779869293
uqx(3,1)=-0.0003299929115003222
uqy(3,1)=0.000195519493609438
uqx(3,2)=-0.0003964101199496858
uqy(3,2)=0.0_dp
uqx(4,1)=-1.209151294434386e-05
uqy(4,1)=7.165004163245931e-06
uqx(4,2)=-1.830230180218659e-05
uqy(4,2)=4.369241246891779e-06
uqx(4,3)=-1.830230180218659e-05
uqy(4,3)=-4.369241246891779e-06
uqx(5,1)=-5.072472977662023e-07
uqy(5,1)=3.005414175984017e-07
uqx(5,2)=-9.393528199777028e-07
uqy(5,2)=3.233772873927477e-07
uqx(5,3)=-1.016494567114082e-06
uqy(5,3)=0.0_dp
uqx(5,4)=-9.393528199777028e-07
uqy(5,4)=-3.233772873927477e-07
endif
if(ntheta.EQ.11)then
uqx(2,1)=-0.01427768000921253
uqy(2,1)=0.008428908168456773
uqx(3,1)=-0.0004619778319639407
uqy(3,1)=0.0002730781672413903
uqx(3,2)=-0.0005426360143594519
uqy(3,2)=0.0_dp
uqx(4,1)=-2.014728185682702e-05
uqy(4,1)=1.191074841386241e-05
uqx(4,2)=-2.970961882633886e-05
uqy(4,2)=7.124796626544694e-06
uqx(4,3)=-2.970961882633886e-05
uqy(4,3)=-7.124796626544694e-06
uqx(5,1)=-1.008004944770031e-06
uqy(5,1)=5.958424731241988e-07
uqx(5,2)=-1.811980725955097e-06
uqy(5,2)=6.251978084599946e-07
uqx(5,3)=-1.942211173240885e-06
uqy(5,3)=0.0_dp
uqx(5,4)=-1.811980725955097e-06
uqy(5,4)=-6.251978084599946e-07
endif
if(ntheta.EQ.12)then
uqx(2,1)=-0.01648841976833508
uqy(2,1)=0.00971385525211049
uqx(3,1)=-0.0006245208306885359
uqy(3,1)=0.0003684360208645578
uqx(3,2)=-0.0007161707066312485
uqy(3,2)=0.0_dp
uqx(4,1)=-3.185339267975134e-05
uqy(4,1)=1.879461360555224e-05
uqx(4,2)=-4.570953565691623e-05
uqy(4,2)=1.101356078295402e-05
uqx(4,3)=-4.570953565691623e-05
uqy(4,3)=-1.101356078295402e-05
uqx(5,1)=-1.86730786016839e-06
uqy(5,1)=1.101634412158447e-06
uqx(5,2)=-3.256105072314334e-06
uqy(5,2)=1.12609910984101e-06
uqx(5,3)=-3.456108802116457e-06
uqy(5,3)=0.0_dp
uqx(5,4)=-3.256105072314334e-06
uqy(5,4)=-1.12609910984101e-06
endif
if(ntheta.EQ.13)then
uqx(2,1)=-0.01874363020388613
uqy(2,1)=0.01102299870360727
uqx(3,1)=-0.000819056300796783
uqy(3,1)=0.0004824043177345613
uqx(3,2)=-0.0009156910682583513
uqy(3,2)=0.0_dp
uqx(4,1)=-4.81174507641687e-05
uqy(4,1)=2.834441813303126e-05
uqx(4,2)=-6.713199189077586e-05
uqy(4,2)=1.62532984542057e-05
uqx(4,3)=-6.713199189077586e-05
uqy(4,3)=-1.62532984542057e-05
uqx(5,1)=-3.254182136977055e-06
uqy(5,1)=1.916674289995359e-06
uqx(5,2)=-5.502334381872796e-06
uqy(5,2)=1.907490988568425e-06
uqx(5,3)=-5.782489009202849e-06
uqy(5,3)=0.0_dp
uqx(5,4)=-5.502334381872796e-06
uqy(5,4)=-1.907490988568425e-06
endif
if(ntheta.EQ.14)then
uqx(2,1)=-0.02101513878443223
uqy(2,1)=0.01234015250707122
uqx(3,1)=-0.001046098478254697
uqy(3,1)=0.0006152602359553973
uqx(3,2)=-0.001138709526361856
uqy(3,2)=0.0_dp
uqx(4,1)=-6.983936154117285e-05
uqy(4,1)=4.108260075599009e-05
uqx(4,2)=-9.46663508367406e-05
uqy(4,2)=2.303123144864916e-05
uqx(4,3)=-9.46663508367406e-05
uqy(4,3)=-2.303123144864916e-05
uqx(5,1)=-5.374846134548899e-06
uqy(5,1)=3.161278639182164e-06
uqx(5,2)=-8.811118519862508e-06
uqy(5,2)=3.06194712722794e-06
uqx(5,3)=-9.167731209069543e-06
uqy(5,3)=0.0_dp
uqx(5,4)=-8.811118519862508e-06
uqy(5,4)=-3.06194712722794e-06
endif
if(ntheta.EQ.15)then
uqx(2,1)=-0.0232773402036116
uqy(2,1)=0.01365064721466022
uqx(3,1)=-0.001305280475347721
uqy(3,1)=0.0007667708872687428
uqx(3,2)=-0.001381810113794881
uqy(3,2)=0.0_dp
uqx(4,1)=-9.785794354653753e-05
uqy(4,1)=5.749524297849592e-05
uqx(4,2)=-0.0001288049100968271
uqy(4,2)=3.148869813455972e-05
uqx(4,3)=-0.0001288049100968271
uqy(4,3)=-3.148869813455972e-05
uqx(5,1)=-8.46637610884184e-06
uqy(5,1)=4.973588375291042e-06
uqx(5,2)=-1.345696760939882e-05
uqy(5,2)=4.687804727529729e-06
uqx(5,3)=-1.386349319225168e-05
uqy(5,3)=0.0_dp
uqx(5,4)=-1.345696760939882e-05
uqy(5,4)=-4.687804727529729e-06
endif
if(ntheta.EQ.16)then
uqx(2,1)=-0.02550806137259484
uqy(2,1)=0.01494181959537789
uqx(3,1)=-0.001595465429771206
uqy(3,1)=0.000936258825655531
uqx(3,2)=-0.001640945873196756
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0001329064560728502
uqy(4,1)=7.80063592484412e-05
uqx(4,2)=-0.0001698108360147586
uqy(4,2)=4.171126842036178e-05
uqx(4,3)=-0.0001698108360147586
uqy(4,3)=-4.171126842036178e-05
uqx(5,1)=-1.278694321771364e-05
uqy(5,1)=7.503853321298472e-06
uqx(5,2)=-1.971030956794678e-05
uqy(5,2)=6.882830706202828e-06
uqx(5,3)=-2.010885042278717e-05
uqy(5,3)=0.0_dp
uqx(5,4)=-1.971030956794678e-05
uqy(5,4)=-6.882830706202828e-06
endif
if(ntheta.EQ.17)then
uqx(2,1)=-0.02768902996623525
uqy(2,1)=0.01620327280964095
uqx(3,1)=-0.001914903012504294
uqy(3,1)=0.001122693458825686
uqx(3,2)=-0.001911747148401908
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0001755820633807187
uqy(4,1)=0.0001029603397676965
uqx(4,2)=-0.0002177119608478438
uqy(4,2)=5.372503557832448e-05
uqx(4,3)=-0.0002177119608478438
uqy(4,3)=-5.372503557832448e-05
uqx(5,1)=-1.860418596133933e-05
uqy(5,1)=1.090765125005818e-05
uqx(5,2)=-2.78198841681919e-05
uqy(5,2)=9.737957708444016e-06
uqx(5,3)=-2.811427798376512e-05
uqy(5,3)=0.0_dp
uqx(5,4)=-2.78198841681919e-05
uqy(5,4)=-9.737957708444016e-06
endif
if(ntheta.EQ.18)then
uqx(2,1)=-0.0298059688900068
uqy(2,1)=0.01742692140290434
uqx(3,1)=-0.002261402797319784
uqy(3,1)=0.001324791852622696
uqx(3,2)=-0.00218979975575966
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0002263307798910572
uqy(4,1)=0.0001326134069799172
uqx(4,2)=-0.0002723167306780222
uqy(4,2)=6.749853746454043e-05
uqx(4,3)=-0.0002723167306780222
uqy(4,3)=-6.749853746454043e-05
uqx(5,1)=-2.618344775700557e-05
uqy(5,1)=1.533904799393413e-05
uqx(5,2)=-3.799806596259685e-05
uqy(5,2)=1.333194246897647e-05
uqx(5,3)=-3.804957109168239e-05
uqy(5,3)=0.0_dp
uqx(5,4)=-3.799806596259685e-05
uqy(5,4)=-1.333194246897647e-05
endif
if(ntheta.EQ.19)then
uqx(2,1)=-0.03184839354005421
uqy(2,1)=0.01860686639767005
uqx(3,1)=-0.002632500871417905
uqy(3,1)=0.001541115305242353
uqx(3,2)=-0.00247086793286604
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0002854464623341651
uqy(4,1)=0.0001671332222719948
uqx(4,2)=-0.0003332456229776692
uqy(4,2)=8.294893631494263e-05
uqx(4,3)=-0.0003332456229776692
uqy(4,3)=-8.294893631494263e-05
uqx(5,1)=-3.5777331916109e-05
uqy(5,1)=2.094453998211611e-05
uqx(5,2)=-5.041053364039542e-05
uqy(5,2)=1.772747608189362e-05
uqx(5,3)=-5.003665347329104e-05
uqy(5,3)=0.0_dp
uqx(5,4)=-5.041053364039542e-05
uqy(5,4)=-1.772747608189362e-05
endif
if(ntheta.EQ.20)then
uqx(2,1)=-0.03380921098354144
uqy(2,1)=0.01973915809034029
uqx(3,1)=-0.003025604178510975
uqy(3,1)=0.001770152796829164
uqx(3,2)=-0.002751053010400932
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0003530806004474008
uqy(4,1)=0.0002066047359070089
uqx(4,2)=-0.0003999707148415089
uqy(4,2)=9.995079644234232e-05
uqx(4,3)=-0.0003999707148415089
uqy(4,3)=-9.995079644234232e-05
uqx(5,1)=-4.761749968378625e-05
uqy(5,1)=2.785831220463549e-05
uqx(5,2)=-6.517066822235473e-05
uqy(5,2)=2.296892009302885e-05
uqx(5,3)=-6.414716777612842e-05
uqy(5,3)=0.0_dp
uqx(5,4)=-6.517066822235473e-05
uqy(5,4)=-2.296892009302885e-05
endif
if(ntheta.EQ.21)then
uqx(2,1)=-0.03568421614472853
uqy(2,1)=0.0208215014107083
uqx(3,1)=-0.003438105278930305
uqy(3,1)=0.002010387168015211
uqx(3,2)=-0.003026890903307585
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0004292590682540321
uqy(4,1)=0.0002510400416038367
uqx(4,2)=-0.0004718571203498731
uqy(4,2)=0.0001183459294007554
uqx(4,3)=-0.0004718571203498731
uqy(4,3)=-0.0001183459294007554
uqx(5,1)=-6.190907493630538e-05
uqy(5,1)=3.61990173774744e-05
uqx(5,2)=-8.233825903866238e-05
uqy(5,2)=2.908155778293243e-05
uqx(5,3)=-8.040403158689358e-05
uqy(5,3)=0.0_dp
uqx(5,4)=-8.233825903866238e-05
uqy(5,4)=-2.908155778293243e-05
endif
if(ntheta.EQ.22)then
uqx(2,1)=-0.03747156089826429
uqy(2,1)=0.02185294737277631
uqx(3,1)=-0.003867466589799954
uqy(3,1)=0.002260343560684938
uqx(3,2)=-0.003295398510176307
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0005139023086210218
uqy(4,1)=0.0003003901861872707
uqx(4,2)=-0.0005482018340347233
uqy(4,2)=0.0001379531343067409
uqx(4,3)=-0.0005482018340347233
uqy(4,3)=-0.0001379531343067409
uqx(5,1)=-7.882756813514927e-05
uqy(5,1)=4.606802067573352e-05
uqx(5,2)=-0.0001019216358938264
uqy(5,2)=3.607207962412092e-05
uqx(5,3)=-9.878582619642953e-05
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0001019216358938264
uqy(5,4)=-3.607207962412092e-05
endif
if(ntheta.EQ.23)then
uqx(2,1)=-0.03917124800252712
uqy(2,1)=0.02283360020288869
uqx(3,1)=-0.00431127715331047
uqy(3,1)=0.002518621941051702
uqx(3,2)=-0.003554081631682294
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0006068462026936415
uqy(4,1)=0.0003545573425025866
uqx(4,2)=-0.0006282673934951395
uqy(4,2)=0.00015857707553736
uqx(4,3)=-0.0006282673934951395
uqy(4,3)=-0.00015857707553736
uqx(5,1)=-9.851794930567258e-05
uqy(5,1)=5.754889134915976e-05
uqx(5,2)=-0.0001238821974354894
uqy(5,2)=4.392995908857313e-05
uqx(5,3)=-0.0001192328797561733
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0001238821974354894
uqy(5,4)=-4.392995908857313e-05
endif
if(ntheta.EQ.24)then
uqx(2,1)=-0.04078467983132041
uqy(2,1)=0.02376435704408403
uqx(3,1)=-0.004767286855073334
uqy(3,1)=0.00278391658984503
uqx(3,2)=-0.003800916556729175
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0007078617681184562
uqy(4,1)=0.0004134062745389956
uqx(4,2)=-0.0007113093010572949
uqy(4,2)=0.0001800159036691447
uqx(4,3)=-0.0007113093010572949
uqy(4,3)=-0.0001800159036691447
uqx(5,1)=-0.0001210953749828021
uqy(5,1)=7.070785140253987e-05
uqx(5,2)=-0.0001481403717251601
uqy(5,2)=5.26293883724199e-05
uqx(5,3)=-0.0001416540841691958
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0001481403717251601
uqy(5,4)=-5.26293883724199e-05
endif
if(ntheta.EQ.25)then
uqx(2,1)=-0.04231427495534361
uqy(2,1)=0.02464668743107288
uqx(3,1)=-0.005233423406853966
uqy(3,1)=0.003055025655109403
uqx(3,2)=-0.004034315476632931
uqy(3,2)=0.0_dp
uqx(4,1)=-0.000816672621751051
uqy(4,1)=0.0004767744840141218
uqx(4,2)=-0.0007965971920023176
uqy(4,2)=0.0002020674940364024
uqx(4,3)=-0.0007965971920023176
uqy(4,3)=-0.0002020674940364024
uqx(5,1)=-0.0001466470646981058
uqy(5,1)=8.559488727257606e-05
uqx(5,2)=-0.0001745822216037299
uqy(5,2)=6.213149817457247e-05
uqx(5,3)=-0.0001659337332567695
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0001745822216037299
uqy(5,4)=-6.213149817457247e-05
endif
if(ntheta.EQ.26)then
uqx(2,1)=-0.0437631542688065
uqy(2,1)=0.02548245325111225
uqx(3,1)=-0.005707796948781806
uqy(3,1)=0.003330853585610454
uqx(3,2)=-0.00425308343192316
uqy(3,2)=0.0_dp
uqx(4,1)=-0.0009329697453595618
uqy(4,1)=0.000544480777040355
uqx(4,2)=-0.0008834303440076407
uqy(4,2)=0.0002245343486965423
uqx(4,3)=-0.0008834303440076407
uqy(4,3)=-0.0002245343486965423
uqx(5,1)=-0.0001752348810775276
uqy(5,1)=0.0001022452652849156
uqx(5,2)=-0.0002030661159222708
uqy(5,2)=7.23866540137963e-05
uqx(5,3)=-0.0001919379139721673
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0002030661159222708
uqy(5,4)=-7.23866540137963e-05
endif
if(ntheta.EQ.27)then
uqx(2,1)=-0.04513489169976175
uqy(2,1)=0.02627376613853282
uqx(3,1)=-0.006188696292743292
uqy(3,1)=0.003610408772998517
uqx(3,2)=-0.0044563721863212
uqy(3,2)=0.0_dp
uqx(4,1)=-0.001056423499571107
uqy(4,1)=0.0006163322244658482
uqx(4,2)=-0.0009711483975435052
uqy(4,2)=0.0002472273003311512
uqx(4,3)=-0.0009711483975435052
uqy(4,3)=-0.0002472273003311512
uqx(5,1)=-0.0002068982555902986
uqy(5,1)=0.0001206812433280477
uqx(5,2)=-0.0002334290810502318
uqy(5,2)=8.333668679649319e-05
uqx(5,3)=-0.0002195201853748941
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0002334290810502318
uqy(5,4)=-8.333668679649319e-05
endif
if(ntheta.EQ.28)then
uqx(2,1)=-0.04643332085491404
uqy(2,1)=0.02702287718328207
uqx(3,1)=-0.006674579953450168
uqy(3,1)=0.003892798220305804
uqx(3,2)=-0.004643634524000381
uqy(3,2)=0.0_dp
uqx(4,1)=-0.001186693080114393
uqy(4,1)=0.0006921296326732865
uqx(4,2)=-0.001059138218690306
uqy(4,2)=0.0002699681973246943
uqx(4,3)=-0.001059138218690306
uqy(4,3)=-0.0002699681973246943
uqx(5,1)=-0.0002416571953113832
uqy(5,1)=0.000140913825684749
uqx(5,2)=-0.0002654926052230902
uqy(5,2)=9.491697014352015e-05
uqx(5,3)=-0.000248526431033852
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0002654926052230902
uqy(5,4)=-9.491697014352015e-05
endif
if(ntheta.EQ.29)then
uqx(2,1)=-0.0476623722906518
uqy(2,1)=0.02773208419857225
uqx(3,1)=-0.007164065345384643
uqy(3,1)=0.004177221193478502
uqx(3,2)=-0.004814580805092504
uqy(3,2)=0.0_dp
uqx(4,1)=-0.001323433931627007
uqy(4,1)=0.0007716718280408433
uqx(4,2)=-0.001146837865060967
uqy(4,2)=0.0002925918557552575
uqx(4,3)=-0.001146837865060967
uqy(4,3)=-0.0002925918557552575
uqx(5,1)=-0.0002795152299791238
uqy(5,1)=0.0001629444802374371
uqx(5,2)=-0.0002990678227357772
uqy(5,2)=0.0001070583234817996
uqx(5,3)=-0.0002787988860265953
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0002990678227357772
uqy(5,4)=-0.0001070583234817996
endif
if(ntheta.EQ.30)then
uqx(2,1)=-0.04882563370160641
uqy(2,1)=0.02840347765698087
uqx(3,1)=-0.00765593696089868
uqy(3,1)=0.004462974022425572
uqx(3,2)=-0.004969132784832358
uqy(3,2)=0.0_dp
uqx(4,1)=-0.00146630725510065
uqy(4,1)=0.0008547611790123954
uqx(4,2)=-0.001233739988517768
uqy(4,2)=0.0003149494088037902
uqx(4,3)=-0.001233739988517768
uqy(4,3)=-0.0003149494088037902
uqx(5,1)=-0.0003204630115092076
uqy(5,1)=0.0001867672361436166
uqx(5,2)=-0.0003339607774787402
uqy(5,2)=0.0001196892313188401
uqx(5,3)=-0.0003101796338676208
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0003339607774787402
uqy(5,4)=-0.0001196892313188401
endif
if(ntheta.EQ.31)then
uqx(2,1)=-0.04992193062183704
uqy(2,1)=0.0290363736961492
uqx(3,1)=-0.008149400568775942
uqy(3,1)=0.004749598900816575
uqx(3,2)=-0.005107282626575934
uqy(3,2)=0.0_dp
uqx(4,1)=-0.001615038562818766
uqy(4,1)=0.0009412378657909575
uqx(4,2)=-0.001319412580829482
uqy(4,2)=0.0003369375238227954
uqx(4,3)=-0.001319412580829482
uqy(4,3)=-0.0003369375238227954
uqx(5,1)=-0.0003644932065908745
uqy(5,1)=0.0002123773891528697
uqx(5,2)=-0.0003699867939099062
uqy(5,2)=0.0001327448853399043
uqx(5,3)=-0.0003425164020803855
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0003699867939099062
uqy(5,4)=-0.0001327448853399043
endif
if(ntheta.EQ.32)then
uqx(2,1)=-0.05090329799667861
uqy(2,1)=0.02960289728720748
uqx(3,1)=-0.008646542228826316
uqy(3,1)=0.00503832273936899
uqx(3,2)=-0.005227869327712846
uqy(3,2)=0.0_dp
uqx(4,1)=-0.001769947347072314
uqy(4,1)=0.001031290132427384
uqx(4,2)=-0.001403673062984919
uqy(4,2)=0.0003587752109830483
uqx(4,3)=-0.001403673062984919
uqy(4,3)=-0.0003587752109830483
uqx(5,1)=-0.0004117279745827106
uqy(5,1)=0.0002398460686580247
uqx(5,2)=-0.0004070748751161033
uqy(5,2)=0.0001462439658497053
uqx(5,3)=-0.0003756924781756205
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0004070748751161033
uqy(5,4)=-0.0001462439658497053
endif
if(ntheta.EQ.33)then
uqx(2,1)=-0.0513606096830928
uqy(2,1)=0.02986539511944403
uqx(3,1)=-0.009170005042813718
uqy(3,1)=0.005342398460657084
uqx(3,2)=-0.005318370702859235
uqy(3,2)=0.0_dp
uqx(4,1)=-0.001935879097096427
uqy(4,1)=0.001127757228004547
uqx(4,2)=-0.001487689894477332
uqy(4,2)=0.0003830788831587224
uqx(4,3)=-0.001487689894477332
uqy(4,3)=-0.0003830788831587224
uqx(5,1)=-0.000463394046633049
uqy(5,1)=0.0002698905192682427
uqx(5,2)=-0.0004460132086709931
uqy(5,2)=0.0001608687247518319
uqx(5,3)=-0.0004098000837104587
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0004460132086709931
uqy(5,4)=-0.0001608687247518319
endif
if(ntheta.EQ.34)then
uqx(2,1)=-0.04874144182744936
uqy(2,1)=0.0283417826549473
uqx(3,1)=-0.009858487467627138
uqy(3,1)=0.005742969922777688
uqx(3,2)=-0.005287936990165755
uqy(3,2)=0.0_dp
uqx(4,1)=-0.002144796524776807
uqy(4,1)=0.001249360187141939
uqx(4,2)=-0.001578285629772355
uqy(4,2)=0.0004267732968133611
uqx(4,3)=-0.001578285629772355
uqy(4,3)=-0.0004267732968133611
uqx(5,1)=-0.0005276562550504655
uqy(5,1)=0.000307285033763345
uqx(5,2)=-0.0004925301179564875
uqy(5,2)=0.0001813270867116877
uqx(5,3)=-0.0004456179148133969
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0004925301179564875
uqy(5,4)=-0.0001813270867116877
endif
if(ntheta.EQ.35)then
uqx(2,1)=-0.03297296564955147
uqy(2,1)=0.01918709821552975
uqx(3,1)=-0.01123505415890843
uqy(3,1)=0.006546108167901785
uqx(3,2)=-0.004684706983380724
uqy(3,2)=0.0_dp
uqx(4,1)=-0.002533504540963353
uqy(4,1)=0.001476247589033027
uqx(4,2)=-0.001676066560316715
uqy(4,2)=0.0005604350432240707
uqx(4,3)=-0.001676066560316715
uqy(4,3)=-0.0005604350432240707
uqx(5,1)=-0.0006414180926895605
uqy(5,1)=0.0003736003362540499
uqx(5,2)=-0.0005690638663033982
uqy(5,2)=0.0002272281056786763
uqx(5,3)=-0.0004791697292710052
uqy(5,3)=0.0_dp
uqx(5,4)=-0.0005690638663033982
uqy(5,4)=-0.0002272281056786763
endif

!!!!!!!!!!!!!!!!!!! Layer 1 at z=0
 
 do i=1,ndim/2

   rx=coord(i,1)
   ry=coord(i,2)

   delta(1)=(0.d0,0.d0)
   delta(2)=(0.d0,0.d0)

!write(*,*) n,m,vk(1),vk(2),rx,ry

   ncount=0
     do nn=1,5,2
       do n=2,5
         do m=1,n-1

           vk(1)=(n-1)*Gn(nn,1)+(m-1)*Gn(nn+2,1)
           vk(2)=(n-1)*Gn(nn,2)+(m-1)*Gn(nn+2,2)

           ux=cos((nn-1)*pi/3.0_dp)*uqx(n,m)-sin((nn-1)*pi/3.0_dp)*uqy(n,m)
           uy=sin((nn-1)*pi/3.0_dp)*uqx(n,m)+cos((nn-1)*pi/3.0_dp)*uqy(n,m)

           delta(1)=delta(1)-2*ux*sin(vk(1)*rx+vk(2)*ry)
           delta(2)=delta(2)-2*uy*sin(vk(1)*rx+vk(2)*ry)

         end do
       end do
     end do

!     if(coord(i,3).LT.deltaR) then 
       coord(i,1)=rx+delta(1)/2
       coord(i,2)=ry+delta(2)/2
!     else
 
 
    coord(i+ndim/2,1)=-coord(i,1)!rx-delta(1)/2
    coord(i+ndim/2,2)=coord(i,2) !ry-delta(2)/2
   
    
enddo


end subroutine Relaxation



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine samplePointsBands(npointsBZ,npointsBack,numk,ncount,g1,g2)

integer(dp), intent(in) :: numk
integer(dp), intent(in) :: ncount
integer(dp), intent(inout) :: npointsBZ(:,:),npointsBack(:,:)
real(dp), intent(in) :: g1(2),g2(2)
 
 
 integer(dp) :: i,ivk,ivk1,ivk2,icount
 real(dp) :: vK1x,vK1y,vK2x,vK2y,vM1x,vM1y,vM2x,vM2y
 real(dp) :: a11,a12,a21,a22,b11,b12,b21,b22,vkx,vky,det
 real(dp) :: aGM, aKG, aMK, aT

 !!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!

 vK1x=(g1(1)+g2(1))/3.0_dp
 vK1y=(g1(2)+g2(2))/3.0_dp

 vK2x=2.0_dp*(g1(1)+g2(1))/3.0_dp
 vK2y=2.0_dp*(g1(2)+g2(2))/3.0_dp

 vM1x=g1(1)/2.0_dp
 vM1y=g1(2)/2.0_dp

 vM2x=g2(1)/2.0_dp
 vM2y=g2(2)/2.0_dp


 aKG=sqrt(vK1x**2+vK1y**2)
 aGM=sqrt(vM1x**2+vM1y**2)
 aMK=sqrt((vK1x-vM1x)**2+(vK1y-vM1y)**2)
 aT=aKG+aGM+aMK

 a11=g1(1)/real(numk,dp)
 a12=g2(1)/real(numk,dp)
 a21=g1(2)/real(numk,dp)
 a22=g2(2)/real(numk,dp)

 det=a11*a22-a12*a21
 b11=a22/det
 b12=-a12/det
 b21=-a21/det
 b22=a11/det


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
icount=0
 
  do ivk2=0,numk/3
  !write(*,*) 'Test1',ivk2
   icount=icount+1
 ivk1=ivk2
npointsBZ(icount,1)=ivk1
npointsBZ(icount,2)=ivk2
npointsBZ(icount,3)=0
npointsBZ(icount,4)=1
npointsBack(npointsBZ(icount,1)+1,npointsBZ(icount,2)+1)=icount
enddo
do ivk=numk/6-1,1,-1
vkx=vM1x+(vK1x-vM1x)*ivk/real(numk/6,dp)
vky=vM1y+(vK1y-vM1y)*ivk/real(numk/6,dp)
!write(*,*) 'Test2',vkx,vky
ivk2=nint(b11*vkx+b12*vky)
ivk1=nint(b21*vkx+b22*vky)
!write(*,*) 'Test2',ivk2,ivk1
vkx=vM1x+(vK1x-vM1x)*ivk/real(numk/6,dp)
vky=vM1y+(vK1y-vM1y)*ivk/real(numk/6,dp)
ivk1=nint(b11*vkx+b12*vky)
ivk2=nint(b21*vkx+b22*vky)
 icount=icount+1
npointsBZ(icount,1)=ivk1
npointsBZ(icount,2)=ivk2
npointsBZ(icount,3)=0
npointsBZ(icount,4)=1
npointsBack(npointsBZ(icount,1)+1,npointsBZ(icount,2)+1)=icount
enddo
do ivk1=numk/2,0,-1
!write(*,*) 'Test3',ivk2
ivk2=0
 icount=icount+1
npointsBZ(icount,1)=ivk1
npointsBZ(icount,2)=ivk2
npointsBZ(icount,3)=0
npointsBZ(icount,4)=1
npointsBack(npointsBZ(icount,1)+1,npointsBZ(icount,2)+1)=icount
enddo
 

 
 
end subroutine samplePointsBands

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine plotBands(bandsTB,fmu,npointsBack,ndimEV,numk,ncount,vq1,vq12)
 
 integer(dp), intent(in)    :: ndimEV,numk,ncount
 integer(dp), intent(in)    :: npointsBack(numk,numk)
 real(dp)   , intent(in)    :: fmu,vq1(2), vq12(2), bandsTB(ndimEV,ncount)
 
 integer(dp) :: i,ivk,ivk1,ivk2
 real(dp) :: vK1x,vK1y,vK2x,vK2y,vM1x,vM1y,vM2x,vM2y
 real(dp) :: a11,a12,a21,a22,b11,b12,b21,b22,vkx,vky,det
 real(dp) :: aGM, aKG, aMK, aT

 !!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!

 vK1x=(vq1(1)+vq12(1))/3.0_dp
 vK1y=(vq1(2)+vq12(2))/3.0_dp

 vK2x=2.0_dp*(vq1(1)+vq12(1))/3.0_dp 
 vK2y=2.0_dp*(vq1(2)+vq12(2))/3.0_dp

 vM1x=vq1(1)/2.0_dp
 vM1y=vq1(2)/2.0_dp

 vM2x=vq12(1)/2.0_dp
 vM2y=vq12(2)/2.0_dp


 aKG=sqrt(vK1x**2+vK1y**2)
 aGM=sqrt(vM1x**2+vM1y**2)
 aMK=sqrt((vK1x-vM1x)**2+(vK1y-vM1y)**2)
 aT=aKG+aGM+aMK

 a11=vq1(1)/real(numk,dp)
 a12=vq12(1)/real(numk,dp)
 a21=vq1(2)/real(numk,dp)
 a22=vq12(2)/real(numk,dp)

 det=a11*a22-a12*a21
 b11=a22/det
 b12=-a12/det
 b21=-a21/det
 b22=a11/det

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!

 do ivk2=numk/3,0,-1
     ivk1=ivk2
     !write(*,*) (numk/3-ivk2)*aKG/real(numk/3,dp)/aT,(ivk1*vq1(:)+ivk2*vq12(:))/real(numk)
     write(99,fmt='(F16.8,3X)', advance="no") (numk/3-ivk2)*aKG/real(numk/3,dp)/aT
     do i=1,ndimEV-1 
     write(99,fmt='(F16.8,3X)', advance="no") bandsTB(i,npointsBack(ivk1+1,ivk2+1))-fmu
     enddo
     write(99,fmt='(F16.8,3X)') bandsTB(ndimEV,npointsBack(ivk1+1,ivk2+1))-fmu
 enddo
 do ivk2=1,numk/2
     ivk1=0
     !write(*,*) (aKG+(ivk2)*aGM/real(numk/2,dp))/aT,(ivk1*vq1(:)+ivk2*vq12(:))/real(numk)
     write(99,fmt='(F16.8,3X)', advance="no") (aKG+(ivk2)*aGM/real(numk/2,dp))/aT
     do i=1,ndimEV-1
     write(99,fmt='(F16.8,3X)', advance="no") bandsTB(i,npointsBack(ivk1+1,ivk2+1))-fmu
     enddo
     write(99,fmt='(F16.8,3X)') bandsTB(ndimEV,npointsBack(ivk1+1,ivk2+1))-fmu
 enddo
 do ivk=1,numk/6
     vkx=vM1x+(vK1x-vM1x)*ivk/real(numk/6,dp)
     vky=vM1y+(vK1y-vM1y)*ivk/real(numk/6,dp)
     ivk2=nint(b11*vkx+b12*vky)
     ivk1=nint(b21*vkx+b22*vky)
     !write(*,*) (aKG+aGM+(ivk)*aMK/real(numk/6,dp))/aT,(ivk1*vq1(:)+ivk2*vq12(:))/real(numk)
     write(99,fmt='(F16.8,3X)', advance="no") (aKG+aGM+(ivk)*aMK/real(numk/6,dp))/aT
     do i=1,ndimEV-1
     write(99,fmt='(F16.8,3X)', advance="no") bandsTB(i,npointsBack(ivk1+1,ivk2+1))-fmu
     enddo
     write(99,fmt='(F16.8,3X)') bandsTB(ndimEV,npointsBack(ivk1+1,ivk2+1))-fmu
 enddo

 close(99)
 
end subroutine plotBands

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine plotBandsGamma(bandsTB,fmu,npointsBack,ndimEV,numk,ncount,vq1,vq12)
 
 integer(dp), intent(in)    :: ndimEV,numk,ncount
 integer(dp), intent(in)    :: npointsBack(numk,numk)
 real(dp)   , intent(in)    :: fmu,vq1(2), vq12(2), bandsTB(ndimEV,ncount)
 
 integer(dp) :: i,ivk,ivk1,ivk2
 real(dp) :: vK1x,vK1y,vK2x,vK2y,vM1x,vM1y,vM2x,vM2y
 real(dp) :: a11,a12,a21,a22,b11,b12,b21,b22,vkx,vky,det
 real(dp) :: aGM, aKG, aMK, aT

 !!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!

 vK1x=(vq1(1)+vq12(1))/3.0_dp
 vK1y=(vq1(2)+vq12(2))/3.0_dp

 vK2x=2.0_dp*(vq1(1)+vq12(1))/3.0_dp
 vK2y=2.0_dp*(vq1(2)+vq12(2))/3.0_dp

 vM1x=vq1(1)/2.0_dp
 vM1y=vq1(2)/2.0_dp

 vM2x=vq12(1)/2.0_dp
 vM2y=vq12(2)/2.0_dp


 aKG=sqrt(vK1x**2+vK1y**2)
 aGM=sqrt(vM1x**2+vM1y**2)
 aMK=sqrt((vK1x-vM1x)**2+(vK1y-vM1y)**2)
 aT=aKG+aGM+aMK

 a11=vq1(1)/real(numk,dp)
 a12=vq12(1)/real(numk,dp)
 a21=vq1(2)/real(numk,dp)
 a22=vq12(2)/real(numk,dp)

 det=a11*a22-a12*a21
 b11=a22/det
 b12=-a12/det
 b21=-a21/det
 b22=a11/det

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!

 do ivk2=0,numk/3
ivk1=ivk2
    write(99,fmt='(F16.8,3X)', advance="no") ivk2*aKG/real(numk/3,dp)/aT
    do i=1,ndimEV-1
    write(99,fmt='(F16.8,3X)', advance="no") bandsTB(i,npointsBack(ivk1+1,ivk2+1))-fmu
    enddo
    write(99,fmt='(F16.8,3X)') bandsTB(ndimEV,npointsBack(ivk1+1,ivk2+1))-fmu
enddo 
do ivk=numk/6,0,-1
vkx=vM1x+(vK1x-vM1x)*ivk/real(numk/6,dp)
vky=vM1y+(vK1y-vM1y)*ivk/real(numk/6,dp)
ivk1=nint(b11*vkx+b12*vky)
ivk2=nint(b21*vkx+b22*vky)
    !write(*,*) '2',(aKG+(numk/6-ivk)*aMK/real(numk/6,dp))/aT,(ivk1*vq1(:)+ivk2*vq12(:))/real(numk)
    write(99,fmt='(F16.8,3X)', advance="no") (aKG+(numk/6-ivk)*aMK/real(numk/6,dp))/aT
    do i=1,ndimEV-1
    write(99,fmt='(F16.8,3X)', advance="no") bandsTB(i,npointsBack(ivk1+1,ivk2+1))-fmu
    enddo
    write(99,fmt='(F16.8,3X)') bandsTB(ndimEV,npointsBack(ivk1+1,ivk2+1))-fmu
enddo
do ivk1=numk/2,0,-1
ivk2=0
!write(*,*) ivk1,ivk2
write(99,fmt='(F16.8,3X)', advance="no") 1.-ivk1*aGM/real(numk/2,dp)/aT
    do i=1,ndimEV-1
    write(99,fmt='(F16.8,3X)', advance="no") bandsTB(i,npointsBack(ivk1+1,ivk2+1))-fmu
    enddo
    write(99,fmt='(F16.8,3X)') bandsTB(ndimEV,npointsBack(ivk1+1,ivk2+1))-fmu
enddo 
close(99)
 
end subroutine plotBandsGamma

end module Hamiltonian
