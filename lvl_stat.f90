!------------------------------------!
!            Trim file names         !
!------------------------------------!

function trim_name(file_name) !Function to trim file names
    implicit none
    character*200 :: trim_name
    character(LEN=*):: file_name
    integer :: ls1, ls2, i

    trim_name=''
    ls1 = len_trim(file_name)
    ls2 = 0
    do i = 1, ls1
        if(file_name(i:i).ne.' ') then
           ls2 = ls2 + 1
           trim_name(ls2:ls2) = file_name(i:i)
        end if
    end do
    return
end function



program lvl_stat

    implicit none
    !!!!!!!!!!!!!!!!!!
    !   Parameters   !
    !!!!!!!!!!!!!!!!!!

    integer, parameter :: v1list = 1
    integer, parameter :: v2list = 0
    integer, parameter :: sites = 10
    integer, parameter :: pts = 5
    integer, parameter :: v2min = 100
    integer, parameter :: v2max = 100
    integer, parameter :: frac = 3
    integer, parameter :: nev = 25
    integer, parameter :: ndis = 1

    character, parameter :: bc*1 = 'p'

    !!!!!!!!!!!!!!!
    !   Scalars   !
    !!!!!!!!!!!!!!!

    integer :: i = 0
    integer :: j = 0
    integer :: k = 0
    integer :: l = 0
    integer :: m = 0
    integer :: temp = 0
    integer :: mom = 0
    integer :: v2i = 0

    double precision :: t = 1.0d0
    double precision :: v1 = 1.0d0
    double precision :: v2 = 1.0d0
    double precision :: eps = 1.0d0
    double precision :: dis = 0.0d0
    double precision :: rmean = 0.0d0

    !!!!!!!!!!!!!!
    !   Arrays   !
    !!!!!!!!!!!!!!

    integer, allocatable :: mom_list(:)

    double precision, allocatable :: v_list(:)
    double precision, allocatable :: dpenergy(:)

    double complex, allocatable :: energy(:,:)

    character :: parameters*200, file_name*200
    character :: trim_name*200

    !!!!!!!!!!!!!!!
    !   Program   !
    !!!!!!!!!!!!!!!

    allocate(energy(2*nev, ndis))
    allocate(dpenergy(2*nev))
    allocate(v_list(376))
    allocate(mom_list(2))

    mom_list = (/ 1, -4 /)

    v_list = 0.d0
    do i = 1, size(v_list)
        v_list(i) = 0.1 * 1.5**(real(i-1)/10.)
!        v_list(i) = 0.1 * 1.5**((i-1)/10.)
!        print*, dble(i-1)/10., 'exp'
!        print*, 1.5**(dble(i-1)/10.), '1.5**0'
!        print*, 0.1 * 1.5**(dble(i-1)/10.), '0.1*'
!        print*, v_list(i), 'v_list(i)'
!        pause
    end do

    do i = v2min, v2max
        file_name = ""
        write (file_name,"('lvl_stat_L=',i0,'N=',i0,'t=',f12.4,'k1=',i0,'k2=',i0,'V2=',i0,'.dat')") sites, pts, t, mom_list(1), mom_list(2), i
        file_name = trim_name(file_name)
        open(99, file=file_name)
        do j = 0, size(v_list) - 1
            energy = 0.0d0
            v1 = v_list(j + 1)

!            v2 = v_list(i)
            v2 = 5.0d0
            do k = 1, size(mom_list)
                mom = mom_list(k)
!                write (parameters,"('L=',i0,'N=',i0,'k=',i0,'t=',f12.4,'V=',i0, &
!                    'V2=',i0,'W=',f12.4,'eps=',f12.4,'ndis=',i0,'BC=',a,'.dat')") sites,pts,mom,t,j, &
!                    i,dis,eps,ndis,bc
                write (parameters,"('L=',i0,'N=',i0,'k=',i0,'t=',f12.4,'V=',f12.4, &
                    'V2=',f12.4,'W=',f12.4,'eps=',f12.4,'ndis=',i0,'BC=',a,'.dat')") sites,pts,mom,t,v1, &
                    v2,dis,eps,ndis,bc
                parameters = trim_name(parameters)
                file_name = "ti_energies_"//parameters

                100 format(1000(F30.20))

                open(11, file=file_name)
                read(11, *) temp
                do l = 1, nev
                    do m = 1, ndis
                        read(11, 100) energy(l + (k-1)*nev, m)

                        dpenergy(l + (k-1)*nev) = dble(real(energy(l + (k-1)*nev, 1)))

                    end do
                end do

                close(11)
            end do
            rmean = 0.d0
            call quicksort_nr(size(dpenergy),dpenergy)
!print*, dpenergy, 'dpenergy'
!print*, size(dpenergy), 'size(dpenergy)'
!print*, v1, 'v1'
!print*, v2, 'v2'
!pause
            call lvlpar (frac, size(dpenergy), dpenergy, rmean) !No loop over disorder realizations
            write(99,*) v1, v2, rmean
        end do
        close(99)
    end do

end program




!----------------------------------------------!
!            Level spacing parameter           !
!----------------------------------------------!

subroutine lvlpar (frac, nev, evals, rmean)

    implicit none

    integer, intent(in) :: frac, nev
    double precision, intent(in) :: evals(nev)
    double precision, intent(out) :: rmean
    double precision :: diff(nev - 2*(int(nev/frac) -1) - 1)
    double precision, allocatable :: r(:)
    integer :: i = 0, cntr = 0

    allocate(r(size(diff) - 1))

    diff = 0
    r = 0
    cntr = 0


    do i = int(nev/frac), nev - int(nev/frac)
        cntr = cntr + 1
        diff(cntr) = evals(i+1) - evals(i)
    end do

    do i = 2, size(diff)
        if (max(diff(i), diff(i - 1)) > 0) then
            r(i - 1) = min(diff(i), diff(i - 1)) / max(diff(i), diff(i - 1))
        else if (i > 1 .and. max(diff(i), diff(i - 1)) == 0) then
            r(i - 1) = 0
        end if
    end do

    rmean = sum(r) / size(r)


end subroutine lvlpar

!!----------------------------------------------!
!!            Level spacing parameter           !
!!----------------------------------------------!
!
!subroutine lvlpar (frac, nev, evals, rmean)
!
!    implicit none
!
!    integer, intent(in) :: frac, nev
!    double precision, intent(in) :: evals(nev)
!    double precision, intent(out) :: rmean
!    double precision :: diff(nev - 2*(int(nev/frac) -1)), r(nev - 2*(int(nev/frac) -1))
!    integer :: i = 0, cntr = 0
!
!    diff = 0
!    r = 0
!    cntr = 0
!
!    do i = int(nev/frac), nev - int(nev/frac) - 1
!        cntr = cntr + 1
!        diff(cntr) = evals(i+1) - evals(i)
!        if (cntr > 1 .and. max(diff(cntr), diff(cntr - 1)) > 0) then
!            r(cntr - 1) = min(diff(cntr), diff(cntr - 1)) / max(diff(cntr), diff(cntr - 1))
!        else if (cntr > 1 .and. max(diff(cntr), diff(cntr - 1)) == 0) then
!            r(cntr - 1) = 0
!        end if
!    end do
!    rmean = sum(r) / size(r)
!
!end subroutine lvlpar



  !---------------------!
  !   Sorting routine   !
  !---------------------!


  ! This version maintains its own stack, to avoid needing to call
  ! itself recursively. By always pushing the larger "half" to the
  ! stack, and moving directly to calculate the smaller "half",
  ! it can guarantee that the stack needs no more than log_2(N)
  ! entries
  subroutine quicksort_nr(dim,array)

    implicit none

    integer, intent(in) :: dim
    double precision, intent(inout) :: array(dim)
    double precision :: temp,pivot
    integer :: i,j,left,right,low,high
    ! If your compiler lacks storage_size(), replace
    ! storage_size(i) by 64
    integer :: stack(2,storage_size(i)),stack_ptr

    low=1
    high=size(array)
    stack_ptr=1


    do

       if (high-low.lt.50) then ! use insertion sort on small arrays
          do i=low+1,high
             temp=array(i)
             do j=i-1,low,-1
                if (array(j).le.temp) exit
                array(j+1)=array(j)
             enddo
             array(j+1)=temp
          enddo
          ! now pop from stack
          if (stack_ptr.eq.1) return
          stack_ptr=stack_ptr-1
          low=stack(1,stack_ptr)
          high=stack(2,stack_ptr)
          cycle
       endif

       ! find median of three pivot
       ! and place sentinels at first and last elements
       temp=array((low+high)/2)
       array((low+high)/2)=array(low+1)
       if (temp.gt.array(high)) then
          array(low+1)=array(high)
          array(high)=temp
       else
          array(low+1)=temp
       endif
       if (array(low).gt.array(high)) then
          temp=array(low)
          array(low)=array(high)
          array(high)=temp
       endif
       if (array(low).gt.array(low+1)) then
          temp=array(low)
          array(low)=array(low+1)
          array(low+1)=temp
       endif
       pivot=array(low+1)

       left=low+2
       right=high-1
       do
          do while(array(left).lt.pivot)
             left=left+1
          enddo
          do while(array(right).gt.pivot)
             right=right-1
          enddo
          if (left.ge.right) exit
          temp=array(left)
          array(left)=array(right)
          array(right)=temp
          left=left+1
          right=right-1
       enddo
       if (left.eq.right) left=left+1
       !          call quicksort(array(1:left-1))
       !          call quicksort(array(left:))
       if (left.lt.(low+high)/2) then
          stack(1,stack_ptr)=left
          stack(2,stack_ptr)=high
          stack_ptr=stack_ptr+1
          high=left-1
       else
          stack(1,stack_ptr)=low
          stack(2,stack_ptr)=left-1
          stack_ptr=stack_ptr+1
          low=left
       endif

    enddo
  end subroutine quicksort_nr
