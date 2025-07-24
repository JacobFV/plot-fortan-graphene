module lapack_routines

implicit none

 
interface diagonalize

 module procedure zdiagRed, cdiagRed, zdiag, cdiag

end interface


contains


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!% Start Diagonalization for complex(8)
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine zdiagRed(A,eigenvalues,ev,il,iu)

  implicit none

  complex(8), intent(inout) :: A(:,:)
  real(8)   , intent(inout) :: eigenvalues(:)
  character(1) , intent(in) :: ev
  integer(8) , intent(in)      :: il , iu
 
  integer(4) :: info, m, lwork, lrwork, liwork,il4,iu4
  real(8)    :: vl, vu, abstol, dlamch
  integer(4), allocatable :: iwork(:), isuppz(:)
  complex(8), allocatable :: work(:),z(:,:)
  real(8)   , allocatable :: rwork(:)

  if( (ev /= 'V') .and. (ev /= 'N') ) then
   write(*,*) 'ev has to be either "V" or "N" '
   stop
  end if

  il4=il

  iu4=iu

  allocate( isuppz(2*(iu-il+1)) )
  allocate( z(size(A,1),iu-il+1) )
  allocate(  work(1) )
  allocate( rwork(1) )
  allocate( iwork(1) )

  lwork = -1
  lrwork= -1
  liwork= -1
  call zheevr(ev,'I','L',size(A,1),A,size(A,1),vl,vu,il4,iu4,dlamch('S'),m,eigenvalues,&
              z,size(A,1),isuppz,work,lwork,rwork,lrwork,iwork,liwork,info)

  lwork  = int( work(1))
  lrwork = int(rwork(1))
  liwork = int(iwork(1))
  deallocate(work,rwork,iwork)
  allocate(work(lwork))
  allocate(rwork(lrwork))
  allocate(iwork(liwork))

  call zheevr(ev,'I','L',size(A,1),A,size(A,1),vl,vu,il4,iu4,dlamch('S'),m,eigenvalues,&
              z,size(A,1),isuppz,work,lwork,rwork,lrwork,iwork,liwork,info)

  A(:,1:iu-il+1) = z

  if(info.ne.0) write(*,*) 'zheev failed' , info

 
end subroutine zdiagRed

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!% Start Diagonalization for complex(4)
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subroutine cdiagRed(A,eigenvalues,ev,il,iu)

  implicit none

  complex, intent(inout) :: A(:,:)
  real   , intent(inout) :: eigenvalues(:)
  character(1) , intent(in) :: ev
  integer , intent(in)      :: il , iu
 
  integer :: info, m, lwork, lrwork, liwork
  real    :: vl, vu, abstol, dlamch
  integer, allocatable :: iwork(:), isuppz(:)
  complex, allocatable :: work(:),z(:,:)
  real   , allocatable :: rwork(:)

  if( (ev /= 'V') .and. (ev /= 'N') ) then
   write(*,*) 'ev has to be either "V" or "N" '
   stop
  end if

  allocate( isuppz(2*(iu-il+1)) )
  allocate( z(size(A,1),iu-il+1) )
  allocate(  work(1) )
  allocate( rwork(1) )
  allocate( iwork(1) )

  lwork = -1
  lrwork= -1
  liwork= -1
  call cheevr(ev,'I','L',size(A,1),A,size(A,1),vl,vu,il,iu,dlamch('S'),m,eigenvalues,&
              z,size(A,1),isuppz,work,lwork,rwork,lrwork,iwork,liwork,info)

  lwork  = int( work(1))
  lrwork = int(rwork(1))
  liwork = int(iwork(1))
  deallocate(work,rwork,iwork)
  allocate(work(lwork))
  allocate(rwork(lrwork))
  allocate(iwork(liwork))

  call cheevr(ev,'I','L',size(A,1),A,size(A,1),vl,vu,il,iu,dlamch('S'),m,eigenvalues,&
              z,size(A,1),isuppz,work,lwork,rwork,lrwork,iwork,liwork,info)

  A(:,1:iu-il+1) = z

  if(info.ne.0) write(*,*) 'cheev failed' , info

 
end subroutine cdiagRed

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!% End Diagonalization for complex(8)
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine zdiagFull(A,eigenvalues,ev)

  implicit none

  complex(8), intent(inout) :: A(:,:)
  real(8)   , intent(inout) :: eigenvalues(:)
  character(1) , intent(in) :: ev
 
  integer(4) :: info, lwork
  complex(8), allocatable :: work(:)
  real(8)   , allocatable :: rwork(:)


  if( (ev /= 'V') .and. (ev /= 'N') ) then
   write(*,*) 'ev has to be either "V" or "N" '
   stop
  end if


  allocate( work(1) )
  allocate( rwork(3*size(A,1)))


  lwork = -1

  call zheev(ev,'L',size(A,1),A,size(A,1),eigenvalues,work,lwork,rwork,info)

  lwork  = int( work(1))
  deallocate(work)
  allocate(work(lwork))

  call zheev(ev,'L',size(A,1),A,size(A,1),eigenvalues,work,lwork,rwork,info)
  if(info.ne.0) write(*,*) 'zheev failed' , info

 
end subroutine zdiagFull

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!% Start Full Diagonalization for complex(4)
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subroutine cdiagFull(A,eigenvalues,ev)

  implicit none

  complex, intent(inout) :: A(:,:)
  real   , intent(inout) :: eigenvalues(:)
  character , intent(in) :: ev
 
  integer :: ndim,info, lwork
  complex, allocatable :: work(:)
  real   , allocatable :: rwork(:)

  if( (ev /= 'V') .and. (ev /= 'N') ) then
   write(*,*) 'ev has to be either "V" or "N" '
   stop
  end if

  allocate( work(1) )
  allocate( rwork(3*size(A,1)))

  lwork = -1

  call cheev(ev,'L',size(A,1),A,size(A,1),eigenvalues,work,lwork,rwork,info)

  lwork  = int( work(1))
  deallocate(work)
  allocate(work(lwork))


  call cheev(ev,'L',size(A,1),A,size(A,1),eigenvalues,work,lwork,rwork,info)
  if(info.ne.0) write(*,*) 'cheev failed' , info

 
end subroutine cdiagFull

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!% End Full Diagonalization for complex(4)
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!% Diagonalization for complex(4)
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine cdiag(A,eigenvalues,ev,d)

  implicit none

  complex(4), intent(inout) :: A(:,:)
  real(4)   , intent(inout):: eigenvalues(:)
  character(1) , intent(in) :: ev
  character(1) , optional , intent(in) :: d
  
  integer(4) :: info
  integer(4), allocatable :: iwork(:)
  complex(4), allocatable :: work(:)
  real(4)   , allocatable :: rwork(:)

 
  if( (ev /= 'V') .and. (ev /= 'N') ) then
   write(*,*) 'ev has to be either "V" or "N" '
   stop
  end if

  if(present(d)) then

    if(ev=='V') then
      allocate( rwork(1+5*size(A,1)+2*size(A,1)**2) )
      allocate(  work(2*size(A,1)+size(A,1)**2) )
      allocate( iwork(3+5*size(A,1)) )
    else
      allocate( rwork(size(A,1)) )
      allocate(  work(size(A,1)+1) )
      allocate( iwork(1) )
    end if

    call cheevd(ev, 'L', size(A,1), A, size(A,1), eigenvalues, work, size(work), rwork, size(rwork), iwork ,size(iwork), info)

    if(info.ne.0) write(*,*) 'cheevd failed' , info

  else

    allocate( rwork(3*size(A,1)) )
    allocate(  work(2*size(A,1)) )

    call cheev(ev, 'L', size(A,1), A, size(A,1), eigenvalues, work, 2*size(A,1), rwork, info)

    if(info.ne.0) write(*,*) 'cheev failed' , info

  end if
  
end subroutine cdiag


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!% Diagonalization for complex(8)
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine zdiag(A,eigenvalues,ev,d)

  implicit none

  complex(8), intent(inout) :: A(:,:)
  real(8)   , intent(inout) :: eigenvalues(:)
  character(1) , intent(in) :: ev
  character(1) , optional , intent(in) :: d
  
  integer(4) :: info
  integer(4), allocatable :: iwork(:)
  complex(8), allocatable :: work(:)
  real(8)   , allocatable :: rwork(:)

  if( (ev /= 'V') .and. (ev /= 'N') ) then
   write(*,*) 'ev has to be either "V" or "N" '
   stop
  end if

  if(present(d)) then

    if(ev=='V') then
      allocate( rwork(1+5*size(A,1)+2*size(A,1)**2) )
      allocate(  work(2*size(A,1)+size(A,1)**2) )
      allocate( iwork(3+5*size(A,1)) )
    else
      allocate( rwork(size(A,1)) )
      allocate(  work(size(A,1)+1) )
      allocate( iwork(1) )
    end if

       call zheevd(ev, 'L', size(A,1), A, size(A,1), eigenvalues, &
                   work, size(work), rwork, size(rwork), iwork ,size(iwork), info)

    if(info.ne.0) write(*,*) 'zheevd failed' , info

  else

    allocate( rwork(3*size(A,1)) )
    allocate(  work(2*size(A,1)) )

    call zheev(ev, 'L', size(A,1), A, size(A,1), eigenvalues, work, 2*size(A,1), rwork, info)

    if(info.ne.0) write(*,*) 'zheev failed' , info

  end if
  
end subroutine zdiag


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!% End Diagonalization for complex(8)
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 


end module
