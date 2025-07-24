! gfortran -lm bilayer LapackRoutines.f90 Hamiltonian.f90 plotBands.f90 -llapack -lblas

! ifort -qopenmp -o test LapackRoutines.f90 Hamiltonian.f90 plotBands.f90 -qmkl -lpthread -lm

! see: https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html


program plotBands

use omp_lib
use Hamiltonian
use lapack_routines

 implicit none
 
 character(250) :: filename
 character(220) :: parameters

 integer(dp) :: n,ivk1, ivk2,icount,ncount,ndimEV,il,iu,nNP 

 
 real(dp) :: aMoire,cs,sn,ang
 
 real(dp) :: t1(2), t2(2), t3(2),tn(6,2) 
 real(dp) :: vq1(2), vq12(2),Gn(8,2)
 real(dp) :: ang_mat(2,2)

 real(dp) :: vk(2)
 real(dp) :: fmu


 integer(dp), allocatable :: npointsBZ(:,:),npointsBack(:,:)
 real(dp), allocatable    :: coord(:,:), coord2(:,:)

 real(dp), allocatable :: evals(:)
 complex(dp), allocatable ::  zh(:,:)
 real(dp), allocatable :: bandsTB(:,:)
 complex(dp), allocatable :: zvx(:,:),zvy(:,:)


 !----------------------------------------------------------------------
 ! Allocate fields
 !----------------------------------------------------------------------

 allocate(coord(ndim,3))
 allocate(coord2(ndim,2))
 allocate(evals(ndim))
 allocate(zh(ndim,ndim)) !Hamiltonian
 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Set geometry Brillouin zone
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! area of the supercell

 aMoire = 3.0_dp*ntheta**2 + 3.0_dp*ntheta + 1.0_dp

! reciprocal lattice vectors of the superlattice

 vq1  = vkd/aMoire*(real(3*ntheta+1,dp)*a1+a2)
 vq12 = vkd/aMoire*(real(3*ntheta+2,dp)*a2-a1)


! angle of rotation cs

 cs = 1.0_dp-1.0_dp/(2.0_dp*aMoire)
 ang    = 0.5_dp*acos(cs)
 ang_mat = reshape([cos(ang),-sin(ang),sin(ang),cos(ang)],[2,2])
 sn = sqrt(1.0_dp-cs**2)

 vq1  = matmul(ang_mat,vq1 )
 vq12 = matmul(ang_mat,vq12)

Gn(1,:)=vq1(:)
Gn(2,:)=vq12(:)
Gn(3,:)=vq12(:)-vq1(:)
Gn(4,:)=-Gn(1,:)
Gn(5,:)=-Gn(2,:)
Gn(6,:)=-Gn(3,:)
Gn(7,:)=Gn(1,:)
Gn(8,:)=Gn(2,:)

write(*,*) 'Twist Angle',2*ang*180/pi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Set geometry real space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Moire lattice parameters

t1 = real( ntheta  ,dp)*a1 + real(  ntheta+1,dp)*a2
t2 = real(-ntheta-1,dp)*a1 + real(2*ntheta+1,dp)*a2
t3 = t2-t1

tn(1,:)=t1
tn(2,:)=t2
tn(3,:)=t3
tn(4,:)=-t1
tn(5,:)=-t2
tn(6,:)=-t3

call WignerSeitz(coord,t1,t2,t3,cs,sn)


!!!!!! Rotate cell to symmetrize around x-axis


 do n=1,ndim
  coord(n,:) = matmul(ang_mat,coord(n,1:2))
 end do


coord2(:,1:2)=coord(:,1:2)

 t1 = matmul(ang_mat,t1)
 t2 = matmul(ang_mat,t2)
 t3 = t2-t1
 
tn(1,:)=t1
tn(2,:)=t2
tn(3,:)=t3
tn(4,:)=-t1
tn(5,:)=-t2
tn(6,:)=-t3
 
!!!!!!!!!!!!!!!!!!!!!!!!!!

if(nrelax.EQ.1)then
call Relaxation(coord,Gn)
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!! Set Points

ndimEV=numb
il=(ndim-2*ndimEV+2*ncb)/2+1
iu=(ndim+2*ncb)/2
nNP=(ndim/2-il)+1

ncount=numk/3+numk/2+numk/6+1
 
allocate(npointsBZ(ncount,8))
 
allocate(npointsBack(numk,numk))

allocate(bandsTB(ndimEV,ncount))

call samplePointsBands(npointsBZ,npointsBack,numk,ncount,vq1,vq12)


 write(parameters,'(A2,I0,A5,I0,A6,I0,A4)') '-I',ntheta,'-numk',numk,'-relax',nrelax,'.dat'
 
 
!$omp parallel do &
!$omp private(vk,icount,zh,evals) &
!$omp shared(coord,tn,vq1,vq12,npointsBZ,ncount,il,iu) &
!$omp shared(bandsTB) 
 do icount = 1,ncount

    vk(:)=(npointsBZ(icount,1)*vq1(:)+npointsBZ(icount,2)*vq12(:))/real(numk,dp)	 

   call ComplexH0(zh,coord,ndim,vk,reshape([t1,t2,t3],[2,3]))
   
   
   call diagonalize(zh,evals,'N',il,iu) ! Calculate only bands from il to iu
   !call diagonalize(zh,evals,'N') ! Calculate all bands

   bandsTB(:,icount)=evals(:)

 end do
!$omp end parallel do

fmu=(bandsTB(nNP,npointsBack(numk/3+1,numk/3+1))+bandsTB(nNP+1,npointsBack(numk/3+1,numk/3+1)))/2._dp

write(*,*) 'Neutrality point',fmu

write(filename,'(A5,A120)') 'Bands',parameters

open(99,file=filename,status='replace')

call plotBandsGamma(bandsTB,fmu,npointsBack,ndimEV,numk,ncount,vq1,vq12)


close(99)



      
end program
