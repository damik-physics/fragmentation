program testfile

    implicit none
    include "omp_lib.h"
    integer :: threads, k, sites, j
    integer(kind=8):: dim, dim2, lr, ll, base(dim2), periods(dim2), transl(dim2), locs(dim2)
    double complex:: state(dim)
    double precision:: eps
    double precision :: entropy
    double precision, allocatable :: singval(:)

    !$omp parallel do default(firstprivate) shared(allstates, dimtot, dim_hs, totaventropy, units, ham_off, ham_di, occupation, basiss, rc_off, rc_di, cstate, cenergy) num_threads(2)


    do j = 1, nest

        if (allocated(weights)) deallocate(weights)
        allocate(weights(ncomp))
        weights = 0.d0

        do i = 1, dimtot
            if(allstates(i,2) < 0) cycle
            weights(comp(allstates(i, 5))) = weights(comp(allstates(i, 5))) +  abs(cstate(allstates(i,2),j,nv+1,nv2+1,ww))**2 * (dble(allstates(i,4)))**(-1)
        end do
        ipr(j) = sum((weights)**2)

        do la = 1, eecut
!                                !$omp critical
            !call ti_centent(num_threads, dim_hs, dimtot, la, sites - la , bjlr, eigstate(1:dim_hs,j,ww), thresh, mom, sites, allstates(1:dimtot,4), allstates(1:dimtot,3), allstates(1:dimtot,2), entropy(la), singval)
            call ti_centent(num_threads, dim_hs, dimtot, la, sites - la , allstates(1:dimtot,1), cstate(1:dim_hs,j,nv+1,nv2+1,ww), thresh, mom, sites, allstates(1:dimtot,4), allstates(1:dimtot,3), allstates(1:dimtot,2), entropy(la), singval)
!                                !$omp end critical
            write(11 + units(thread_num + 1, thread_num2 + 1), *) entropy(la)
            !write(81 + units(thread_num + 1, thread_num2 + 1), *) singval
            aventropy(la) = aventropy(la) + entropy(la)
            if (allocated(singval)) deallocate(singval)
            if (allocated(bjlr)) deallocate(bjlr)
        end do
    end do

end program



subroutine ti_centent(threads, dim, dim2, lr, ll, base, state, eps, k, sites, periods, transl, locs, entropy, singval)

    implicit none

    integer, intent(in) :: threads, k, sites
    integer(kind=8), intent(in) :: dim, dim2, lr, ll, base(dim2), periods(dim2), transl(dim2), locs(dim2)
    double complex, intent(in) :: state(dim)
    double precision, intent(in) :: eps
    double precision, intent(out) :: entropy
    double precision, allocatable , intent(out):: singval(:)

    double complex :: c(2**ll, 2**lr)
    double complex :: cd(2**lr, 2**ll)
    double complex :: ccd(2**lr, 2**lr)
    double complex :: phase
    integer :: loop = 0, loop2 = 0
    logical :: check
    double precision :: delta  = 0.0000001d0


    if (allocated(singval)) deallocate(singval)
    c = (0.d0, 0.d0)
    cd = (0.d0, 0.d0)
    ccd = (0.d0, 0.d0)


    loop = 0

!    !$omp critical
    do loop = 1, 252!size(base)
        if (periods(loop) < 0 .or. locs(loop) < 1) cycle
        c(base(loop)/(2**lr) + 1, mod(base(loop),(2**lr)) + 1) = state(locs(loop)) * sqrt(dble(periods(loop)))/dble(sites) * exp(-2*pi*ii*k*transl(loop)/sites)
        cd(mod(base(loop),(2**lr)) + 1, base(loop)/(2**lr) + 1) = dconjg(c(base(loop)/(2**lr) + 1, mod(base(loop),(2**lr)) + 1))
    end do
!    !$omp end critical


    ccd = 0.d0
    ccd = matmul(cd,c)


    call testhermitian (check, 2**lr, ccd, delta)

    !!$omp critical
    print*, 'Diagonalizing...'
    call cfulldiag(.False., 2**lr, ccd, singval)
    print*, 'Finished diagonalization.'
    !!$omp end critical
    entropy = 0.0d0

    !!$omp parallel do num_threads(threads)
    do loop2 = 1, 2**lr
        if (singval(loop2) > eps) then
            !$omp atomic
            entropy = entropy - singval(loop2) * dlog(singval(loop2))
        end if
    end do
    !!$omp end parallel do

end subroutine ti_centent
