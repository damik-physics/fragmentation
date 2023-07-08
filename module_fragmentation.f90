module module_fragmentation

implicit none
include "omp_lib.h"

double precision, parameter :: pi = 4*atan(1.d0)
double complex, parameter   :: ii = (0, 1)

contains


!----------------------------------------------------!
!            Spinless: Set initial variables         !
!----------------------------------------------------!

subroutine slsetvars(dim, dynamic, nested, nthreads, n, pts, threads, filling)

    implicit none

    integer, intent(in) :: dim
    integer, intent(in) :: nthreads
    double precision, intent(in) :: filling
    logical, intent(in) :: dynamic
    logical, intent(in) :: nested

    integer, intent(out) :: n
    integer, intent(out) :: pts
    integer, intent(out) :: threads


    !$ call omp_set_dynamic(dynamic)
    !$ call omp_set_nested(nested)

    call KMP_SET_STACKSIZE_S(1000000000)

    n = nex

    pts = int(filling*n)
    threads = nthreads

end subroutine




!------------------------------------------!
!            Spinless basis states         !
!------------------------------------------!

subroutine slbasis (io, sites, pts, dim, permuts)

    implicit none

    integer, intent(in) :: sites, pts, io
    integer(kind=8), intent(out) :: dim
    integer(kind=8), allocatable, intent(out) :: permuts(:)

    integer :: a(sites)
    integer :: temp = 0, I_in = 0, I_tail = 0
    integer :: i = 0, j = 0, l = 0, q = 0
    character :: file_name*100

    dim = int(fact(sites) / (fact(pts) * fact(max(1,sites-pts))),8) !Number of basis states

    if (dim == 0) then
        print*, 'No basis states available.'
        go to 112
    end if

    !Permutations contains integer values I of basis states, perm_up/dn contain the integer values I_up/dn and the corresponding J_up/dn
    if (allocated(permuts)) deallocate(permuts)
    allocate(permuts(dim))
    permuts = 0

    a(1 : sites-pts) = 0 !'a' contains the initial basis state with all '1's to the right
    a(sites-pts+1 : sites) = 1
    I_in = 0
    do  l = 1, sites !Binary representation of 'a'
        I_in = I_in + a(l) * 2**(sites-l)
    end do


    do j = 1, dim !Generates all possible configurations in 'permuts'
        I_tail = 0
        permuts(j) = I_in
        temp = I_in
        do i = 0, 64
            if (btest(temp,i)) then
                temp = ibclr(temp,i)
                if (not(btest(temp,i+1))) then   ! asks if pos i+1 is zero
                    temp = ibset(temp,i+1)
                    exit
                end if
            I_tail = ibset(ishft(I_tail,1),0) ! only generated if first loop finds no '01' pair and has to be repeated
            end if
        end do
        I_in = temp + I_tail
    end do

    print*, 'Basis successfully generated'
    print*, ''

    file_name=''
    write(file_name,"('basis_L=',i0,'N=',i0,'.dat')") sites, pts
    open(11 + io, file = trim_name(file_name))
    write(11 + io,*) dim
    write(11 + io,*) permuts(1:dim)
    close(11 + io)

    112 continue

end subroutine slbasis




!-------------------------------------!
!         Test hermiticity            !
!-------------------------------------!

subroutine testhermitian (hermit, dim, mat, eps)

    implicit none

    integer(kind=8), intent(in) :: dim
    double complex, intent(in) :: mat(dim, dim)
    double precision :: eps
    logical, intent(out) :: hermit

    integer :: i = 0 , j = 0

    !$OMP THREADPRIVATE (i, j)
    i = 0
    j = 0

    hermit=.true.
    do i = 1, dim
        do j = 1, dim

            if(real(conjg(mat(i,j))) - real(mat(j,i)) > eps .or. aimag(conjg(mat(i,j))) - aimag(mat(j,i)) > eps) then
                if(real(conjg(mat(i,j))) - real(mat(j,i)) > eps) then
                    print*, real(conjg(mat(i,j))), 'real(conjg(mat(i,j)))'
                    print*, real(mat(j,i)), 'real(mat(j,i))'
                    print*,  ''
                end if

                if(aimag(conjg(mat(i,j))) - aimag(mat(j,i)) > eps) then
                    print*, aimag(conjg(mat(i,j))), 'aimag(conjg(mat(i,j)))'
                    print*, aimag(mat(j,i)), 'aimag(mat(j,i))'
                    print*,  ''
                end if


                print*, i, 'i'
                print*, j, 'j'
                print*, mat(j,i), 'mat(j,i)'
                print*, mat(i,j), 'mat(i,j)'
                print*,  ''

                hermit=.false.
                exit
            endif
        end do
    end do
    if(.not. hermit) then
        print*, 'MATRIX IS NON-HERMITIAN!'
        print*, ''
        pause
    end if
    !else
    !    print*, 'MATRIX IS HERMITIAN'
    !    print*, ''
    !end if
    return



end subroutine testhermitian



!-----------------------------------------------------------------------!
!            1D: Translation representatives and periodicity            !
!-----------------------------------------------------------------------!

subroutine checkstate (s, n, k, r)
    !Checks whether state |s> is a new valid representative basis state with momentum compatibility. 'n' is the number of sites, 'k' the momentum and 'r' the periodicity of the state.
    implicit none
    integer, intent(in) :: n, k
    integer(kind=8), intent(in) :: s
    integer(kind=8), intent(out) :: r

    integer :: t, i

    r = -1
    t = s
    do i = 1, n
        t = ishftc(t, 1, n) !Translate by one to the left
        if (t < s) then !Representative is already in the list
            return
        else if (t == s) then !New potential representative found
            if (modulo(k, n/i) .ne. 0.d0) return !Check compatibility with 'momentum' k \in {-L/2 + 1, ..., L/2}
            r = i !Periodicity of state |s>
            return
        end if
    end do

end subroutine checkstate


!---------------------------------------------!
!            Momentum state basis             !
!---------------------------------------------!

subroutine momentumbasis (dim, n, k, basis, momdim, mombasis, periods, nonreps)

    implicit none
    integer(kind=8), intent(in) :: dim, basis(dim)
    integer, intent(in) :: n, k
    integer(kind=8), intent(out) :: momdim
    integer(kind=8), allocatable, intent(out) :: mombasis(:), periods(:), nonreps(:)

    integer :: j = 0, cntr = 0, cntr2 = 0
    integer(kind=8) :: r = 0
    integer :: mombasis_temp(dim), periods_temp(dim), nonreps_temp(dim)


    do j = 1, dim
        call checkstate(basis(j), n, k, r)
        if (r >= 0) then
            cntr = cntr + 1
            mombasis_temp(cntr) = basis(j)
            periods_temp(cntr) = r
        !else if (r == -1) then !Save non-representative states for calculation of entanglement entropy
        else !Save non-representative states for calculation of entanglement entropy
            cntr2 = cntr2 + 1
            nonreps_temp(cntr2) = basis(j)
        end if
    end do
    momdim = cntr
    allocate(mombasis(momdim), periods(momdim), nonreps(cntr2))
    mombasis(1:momdim) = mombasis_temp(1:momdim)
    periods(1:momdim) = periods_temp(1:momdim)
    nonreps(1:cntr2) = nonreps_temp(1:cntr2)


end subroutine momentumbasis





!---------------------------------------------!
!            Find representative              !
!---------------------------------------------!

subroutine representative (s, n, r, l)
    !Finds the representative 'r' state for state 's'. 'n' is the number of sites and 'l' the number of shifts needed to translate 's' to 'r'.
    implicit none
    integer, intent(in) :: n
    integer(kind=8), intent(in) :: s
!    integer, intent(out) :: r
    integer(kind=8), intent(out) :: l, r

    integer :: t, i

    r = s
    t = s
    l = 0
    do i = 1, n-1
        t = ishftc(t, 1, n)
        if (t < r) then
            r = t
            l = i
        end if
    end do

end subroutine representative





!------------------------------------------------!
!            Lookup of basis states              !
!------------------------------------------------!


subroutine findstate (dim, s, basis, loc)

    implicit none

    integer(kind=8), intent(in) :: dim, basis(dim)
    integer(kind=8), intent(in) :: s
    integer(kind=8), intent(out) :: loc

    integer :: left, right, mean

    left = 1
    right = dim
    do while (left <= right)
        mean = floor((left + right)/2.d0)
        if (basis(mean) < s) then
            left = mean + 1
        else if (basis(mean) > s) then
            right = mean - 1
        else
            loc = mean
            return
        end if
    end do
    loc = -1 !If no state is found
    return


end subroutine findstate














!------------------------------------------------------!
!            Spinless: Hopping Hamiltonian             !
!------------------------------------------------------!

subroutine slhopping_serial(dim, sites, bc, ti, ip, eps, t, basis, rc, ham, nnz, mom, periods, adj)

    implicit none

    integer(kind=8), intent(in) :: dim, basis(dim)
    integer, intent(in) :: sites, ti, ip
    integer, intent(in), optional :: mom
    integer(kind=8), intent(in), optional :: periods(dim)
    double precision, intent(in) :: t, eps
    character, intent(in) :: bc*1

    integer, intent(out) :: adj(dim, dim)
    integer, allocatable, intent(out) :: rc(:,:)
    double complex, allocatable, intent(out) :: ham(:)
    integer, intent(inout) :: nnz

    double complex, parameter :: ii = (0, 1)

    integer :: pntr = 0
    integer :: mask = 0
    integer(kind=8) :: loc = 0, loc2 = 0
    integer(kind=8) :: ntrans = 0
    integer(kind=8) :: rep = 0
    integer :: n_temp
    integer :: chunk = 0
    integer :: remainder = 0
    integer :: arrsize = 0
    integer :: n_threads = 0
    integer :: my_stat = 0
    integer :: i = 0 , j = 0 , m = 0 , k = 0 , l = 0 , q = 0
    integer :: l1 = 0
    integer :: cntr = 0 , cntrj = 0 , dbl = 0
    integer :: vbonds(2)
    integer, allocatable :: rc_temp(:,:)
    integer, allocatable :: temp_rc(:,:)
    double complex, allocatable :: ham_temp(:)
    double complex, allocatable:: temp(:)


    print*, 'Building hopping Hamiltonian ...'
    print*, ''



    adj = 0
    n_temp = 0
    arrsize = 2*dim*sites
    if (allocated(temp)) deallocate(temp)
    if (allocated(temp_rc)) deallocate(temp_rc)
    allocate(temp(arrsize))
    allocate(temp_rc(arrsize,2))

    do j = 1, size(basis)
        cntrj = n_temp      !Starting point in temp_rc for matrix elements of row 'j' to look for double entries with the same column 'loc'
        cntr  = 0           !Counts the number of already calculated matrix elements of each row 'j'
        do m = 0, sites - 1 ! m goes through all digits of basis states
            if (m == sites - 1 .and. bc /= 'p' ) cycle
            mask = ibset( ibset(0, m), mod(m + 1, sites ) )  !Procedure suggested in Lin-paper. Sets a mask with only non-zero components on sites i=m+1 and i=m+2
            k = iand( mask, basis( j ) ) !K records the occupancy of those two sites with up-spins
            l = ieor( mask, k ) !L records whether hopping is allowed or not. If it's 0 or the same as the mask, hopping is not allowed. If it is allowed the occupations of both sites (01 or 10) are swapped (10 or 01)
            if (l == 0 .or. l == mask ) then
                cycle
            end if

            if( ti == 1 .and. bc == 'p' ) then
                call representative(basis(j) - k + l, sites, rep, ntrans) !Finds the representative of scattered state in momentum orbit and determines the number of translations 'ntrans' needed to map to representative.
            else
                rep = basis(j) - k + l
                ntrans = 0
            end if
            call findstate(dim, rep, basis, loc) !Finds the location of representative in basis

            if (loc > 0) then
                dbl = 0 !Flag for whether a new matrix element should be created (0) or an old matrix element was updated (1)
                cntr = cntr + 1
                if (cntr > 1 .and. ti == 1) then
                    do i = max(cntrj,1), cntrj + cntr !Loop for updating existing elements: Goes through all matrix elements already calculated for row 'j' and checks whether column 'loc' already exists.
                        if(temp_rc(i,2) == loc .and. temp_rc(i,1) == j) then !If yes, update existing element.
                            if (m == sites-1 .and. bc == 'p' ) then !PBC hopping across boundary.
                                temp(i) = temp(i) + (-1)**(popcnt(basis(j))-1)*(-1)*t*eps* (sqrt(real(periods(j))/real(periods(loc))) * exp(-ii*2*pi*mom*ntrans/sites))**ti
                                !temp(i) = temp(i) + (-1)*t*eps* (sqrt(real(periods(j))/real(periods(loc))) * exp(-ii*2*pi*mom*ntrans/sites))**ti
                            else
                                temp(i) = temp(i) + (-1)*t*eps* (sqrt(real(periods(j))/real(periods(loc))) * exp(-ii*2*pi*mom*ntrans/sites))**ti
                            end if
                            dbl = 1 !Found and updated existing element.
                        end if
                    end do
                end if
                if (dbl == 0) then !If no existing element was found, create new matrix element.
                    n_temp = n_temp + 1 !Counter for non-zero matrix elements.
                    if (m == sites-1 .and. bc == 'p' ) then !PBC hopping
                        temp(n_temp) = temp(n_temp) + (-1)**(popcnt(basis(j))-1)*(-1)*t*eps* (sqrt(real(periods(j))/real(periods(loc))) * exp(-ii*2*pi*mom*ntrans/sites))**ti
                    else
                        temp(n_temp) = temp(n_temp) + (-1)*t*eps* (sqrt(real(periods(j))/real(periods(loc))) * exp(-ii*2*pi*mom*ntrans/sites))**ti
                    end if
                    temp_rc(n_temp,1) = j !Row of matrix element
                    temp_rc(n_temp,2) = loc !Column of matrix element
                    if (ti == 0 .and. ip == 1 .and. btest(basis(j), modulo(m - 1, sites) ) == btest(basis(j), modulo(m + 2, sites))) then
                        adj(j,loc) = 1
                    end if
                    if (ti == 0 .and. ip == 2 .and. btest(basis(j), modulo(m - 2, sites) ) + btest(basis(j), modulo(m + 2, sites) ) == btest(basis(j), modulo(m - 1, sites)) + btest(basis(j), modulo(m + 3, sites)) ) then
                        adj(j,loc) = 1
                    end if
                    !if (ti == 0 .and. ip == 1) then
                    !    call findstate(dim, ishftc(basis(j), 1, sites), basis, loc2)
                    !    adj(j,loc2) = 1
                    !end if
                end if
            end if
        end do
    end do

    nnz = n_temp
    if (allocated(ham_temp)) deallocate(ham_temp)
    if (allocated(rc_temp)) deallocate(rc_temp)
    allocate(rc_temp(nnz,2))
    allocate(ham_temp(nnz))

    ham_temp(1:nnz) = temp(1:nnz)
    rc_temp(1:nnz,1) = temp_rc(1:nnz,1)
    rc_temp(1:nnz,2) = temp_rc(1:nnz,2)

    call move_alloc(ham_temp,ham)
    call move_alloc(rc_temp,rc)

    if (allocated(ham_temp)) deallocate(ham_temp)
    if (allocated(rc_temp)) deallocate(rc_temp)
    if (allocated(temp)) deallocate(temp)
    if (allocated(temp_rc)) deallocate(temp_rc)

    print*, 'Finished hopping Hamiltonian.'
    print*, ''

end subroutine slhopping_serial





!-----------------------------------------!
!            Adjacency matrix (COO)       !
!-----------------------------------------!


subroutine adjacency_coo( ip, dim, sites, bc, basis, nnz, adj)

    implicit none

    integer(kind=8), intent(in) :: dim, basis(dim)
    integer, intent(in) :: ip, sites
    character, intent(in) :: bc*1

    integer, intent(out) :: nnz
    integer, allocatable, intent(out) :: adj(:,:)

    integer(kind=8) :: row     = 0
    integer(kind=8) :: col     = 0
    integer(kind=8) :: scat    = 0
    integer(kind=8) :: arrsize = 0
    integer :: mask = 0, cntr = 0
    integer :: i = 0, m = 0, k = 0, l = 0
    integer, allocatable :: adj_coo(:,:)

    cntr    = 0
    arrsize = 2 * dim * sites
    if (allocated( adj_coo ) ) deallocate( adj_coo )
    allocate( adj_coo( arrsize, 2 ) )

    do row = 1, dim
        do m = 0, sites - 1 ! m goes through all digits of basis states
            if ( m == sites - 1 .and. bc /= 'p' ) cycle
            mask = ibset( ibset( 0, m ), mod( m + 1, sites ) )  !Procedure suggested in Lin-paper. Sets a mask with only non-zero components on sites i=m+1 and i=m+2
            k = iand( mask, basis( row ) ) !K records the occupancy of those two sites with up-spins
            l = ieor( mask, k ) !L records whether hopping is allowed or not. If it's 0 or the same as the mask, hopping is not allowed. If it is allowed the occupations of both sites (01 or 10) are swapped (10 or 01)
            if (l == 0 .or. l == mask ) then
                cycle
            end if
            scat = basis(row) - k + l
            call findstate( dim, scat, basis, col ) !Finds the location of representative in basis

            if (ip == 1 .and. btest( basis(row), modulo(m - 1, sites) ) == btest( basis(row), modulo( m + 2, sites ) ) ) then
                cntr = cntr + 1
                adj_coo(cntr,1) = row
                adj_coo(cntr,2) = col
            end if
            if (ip == 2 .and. btest( basis(row), modulo( m - 2, sites ) ) + btest( basis(row), modulo( m + 2, sites ) ) &
                                         == btest(basis(row), modulo(m - 1, sites) ) + btest(basis(row), modulo(m + 3, sites) ) ) then
                cntr = cntr + 1
                adj_coo( cntr, 1 ) = row
                adj_coo( cntr, 2 ) = col
            end if
        end do
    end do

    nnz = cntr
    if( allocated( adj ) ) deallocate( adj )
    allocate( adj( nnz, 2 ) )
    do i = 1, nnz
        adj( i, 1 ) = adj_coo( i, 1 )
        adj( i, 2 ) = adj_coo( i, 2 )
    end do
    if( allocated( adj_coo ) ) deallocate( adj_coo )

    return

end subroutine adjacency_coo






!-----------------------------------------!
!            Adjacency matrix  (CSR)      !
!-----------------------------------------!


subroutine adjacency_csr(ti, ip, dim, sites, bc, basis, nnz, adjpntr, adjcol)

    implicit none

    integer(kind=8), intent(in) :: dim, basis(dim)
    integer, intent(in) :: ti, ip, sites
    character, intent(in) :: bc*1

    integer, intent(out) :: nnz
    integer, allocatable, intent(out) :: adjpntr(:), adjcol(:)

    integer(kind=8) :: row     = 0
    integer(kind=8) :: col     = 0
    integer(kind=8) :: scat    = 0
    integer(kind=8) :: arrsize = 0
    integer :: mask = 0, cntr  = 0
    integer :: i = 0, m = 0, k = 0, l = 0
    integer, allocatable :: adj_coo(:,:)
    double precision, allocatable :: dummy(:)

    cntr = 0
    arrsize = 2*dim*sites
    if (allocated(adj_coo)) deallocate(adj_coo)
    allocate(adj_coo(arrsize,2))


    do row = 1, dim
        do m = 0, sites - 1 ! m goes through all digits of basis states
            if (m == sites - 1 .and. bc /= 'p' ) cycle
            mask = ibset( ibset(0, m), mod(m + 1, sites ) )  !Procedure suggested in Lin-paper. Sets a mask with only non-zero components on sites i=m+1 and i=m+2
            k = iand( mask, basis( row ) ) !K records the occupancy of those two sites with up-spins
            l = ieor( mask, k ) !L records whether hopping is allowed or not. If it's 0 or the same as the mask, hopping is not allowed. If it is allowed the occupations of both sites (01 or 10) are swapped (10 or 01)
            if (l == 0 .or. l == mask ) then
                cycle
            end if
            scat = basis(row) - k + l
            call findstate(dim, scat, basis, col) !Finds the location of representative in basis

            if (ti == 0 .and. ip == 1 .and. btest(basis(row), modulo(m - 1, sites) ) == btest(basis(row), modulo(m + 2, sites))) then
                cntr = cntr + 1
                adj_coo(cntr,1) = row
                adj_coo(cntr,2) = col
            end if
            if (ti == 0 .and. ip == 2 .and. btest(basis(row), modulo(m - 2, sites) ) + btest(basis(row), modulo(m + 2, sites) ) &
                                         == btest(basis(row), modulo(m - 1, sites)) + btest(basis(row), modulo(m + 3, sites)) ) then
                cntr = cntr + 1
                adj_coo(cntr,1) = row
                adj_coo(cntr,2) = col
            end if
        end do
    end do

    nnz = cntr
    if( allocated( adjpntr ) ) deallocate( adjpntr )
    if( allocated( adjcol  ) ) deallocate( adjcol  )
    if( allocated( dummy   ) ) deallocate( dummy   )
    allocate( adjpntr( dim + 1 ) )
    allocate( adjcol( nnz ) )
    allocate( dummy( nnz ) )

    call coocsr( dim, nnz, dummy, adj_coo(1:nnz,1), adj_coo(1:nnz,2), dummy, adjcol, adjpntr ) !Creates CSR-format (ao,jao,iao) sparse matrix from coordinate format (a,ir,jc)

    if( allocated( adj_coo ) ) deallocate( adj_coo )
    if( allocated( dummy ) ) deallocate( dummy )

    return

end subroutine adjacency_csr



!---------------------------------------------------------!
!            Find connected components (sparse)           !
!---------------------------------------------------------!

!Eric's code
subroutine connected_components(NNZ, N, indices, indptr, labels)
    ! Label the connected components of a CSR or CSC matrix.
    !
    ! Input parameters
    ! ----------------
    ! NNZ : integer
    !     Number of non-zero entries.
    ! N : integer
    !     Number of rows or columns.
    ! indices : integer(NNZ)
    !     Row/column index of non-zero values.
    ! indptr : integer(N+1)
    !     Indexes where a given row/column starts.
    !
    ! Output parameters
    ! -----------------
    ! labels : integer(N)
    !     Array of labels of the connected components.
    !     Nodes with the same label are connected.
    integer, parameter :: VOID = -1
    integer, intent(in) :: NNZ
    integer(kind=8), intent(in) :: N
    integer, intent(in) :: indices(NNZ), indptr(N+1)
    integer, intent(out) :: labels(N)
    if (indptr(1) == 0) then
        call connected_components_0(NNZ, N, indices, indptr, labels)
    else if (indptr(1) == 1) then
        call connected_components_1(NNZ, N, indices, indptr, labels)
    else
        labels = VOID
    end if

end subroutine connected_components

subroutine connected_components_0(NNZ, N, indices, indptr, labels)
    ! Label the connected components of a CSR or CSC matrix (0-based indexing).
    integer, parameter :: VOID = -1, LAST = -2
    integer, intent(in) :: NNZ
    integer(kind=8), intent(in) :: N
    integer, intent(in) :: indices(0:NNZ-1), indptr(0:N)
    integer, intent(out) :: labels(0:N-1)
    integer :: label, row, col, ptr, head, v
    labels = VOID
    label = 0
    do row = 0, N-1
        if (labels(row) == VOID) then
            ! stack.push(row)
            head = row
            labels(row) = LAST
            do while (head /= LAST)
                ! v = stack.pop()
                v = head
                head = labels(v)
                ! Traverse edges
                labels(v) = label
                do ptr = indptr(v), indptr(v+1)-1
                    col = indices(ptr)
                    if (labels(col) == VOID) then
                        ! stack.push(col)
                        labels(col) = head
                        head = col
                    end if
                end do
            end do
            label = label + 1
        end if
    end do
end subroutine connected_components_0

subroutine connected_components_1(NNZ, N, indices, indptr, labels)
    ! Label the connected components of a CSR or CSC matrix (1-based indexing).
    integer, parameter :: VOID = -1, LAST = -2
    integer, intent(in) :: NNZ
    integer(kind=8), intent(in) :: N
    integer, intent(in) :: indices(NNZ), indptr(N+1)
    integer, intent(out) :: labels(N)
    integer :: label, row, col, ptr, head, v
    labels = VOID
    label = 1
    do row = 1, N
        if (labels(row) == VOID) then
            ! stack.push(row)
            head = row
            labels(row) = LAST
            do while (head /= LAST)
                ! v = stack.pop()
                v = head
                head = labels(v)
                ! Traverse edges
                labels(v) = label
                do ptr = indptr(v), indptr(v+1)-1
                    col = indices(ptr)
                    if (labels(col) == VOID) then
                        ! stack.push(col)
                        labels(col) = head
                        head = col
                    end if
                end do
            end do
            label = label + 1
        end if
    end do
end subroutine connected_components_1




!--------------------------------------------------------------!
!            Calculate connected components (sparse)           !
!--------------------------------------------------------------!

subroutine scc(dim, sites, bc, basis, ncompv1, ncompv2, compv1, compv2)

    integer(kind=8), intent(in) :: dim, basis(dim)
    integer, intent(in) :: sites
    character, intent(in) :: bc*1
    integer, intent(out) :: ncompv1, ncompv2
    integer, allocatable, intent(out) :: compv1(:), compv2(:)

    integer :: nnz   = 0
    integer :: ncomp = 0
    integer :: k = 0, i = 0, j = 0, cntr = 0
    integer, allocatable :: comp(:), adj(:,:), adj_temp(:,:), adjpntr(:), adjcol(:), adj_coo(:,:), compsize(:,:)
    double precision, allocatable :: dummy(:)

    do k = 1, 2
        if ( not( allocated( comp ) ) )  allocate( comp( dim ) )
        comp = 0
        !Generate first adjacency matrix
        !call adjacency_csr( dim, sites, bc, basis, adjpntr, adjcol )
        call adjacency_coo(k, dim, sites, bc, basis, nnz, adj_temp)

        if( allocated( adjpntr ) ) deallocate( adjpntr )
        if( allocated( adjcol  ) ) deallocate( adjcol  )
        if( allocated( dummy   ) ) deallocate( dummy   )
        allocate( adjpntr( dim + 1 ) )
        allocate( adjcol( nnz ) )
        allocate( dummy( nnz ) )

        call coocsr( dim, nnz, dummy, adj_temp(1:nnz,1), adj_temp(1:nnz,2), dummy, adjcol, adjpntr ) !Creates CSR-format (ao,jao,iao) sparse matrix from coordinate format (a,ir,jc)

        if( allocated( dummy ) ) deallocate( dummy )

        call connected_components( nnz, dim, adjcol, adjpntr, comp )
        if( allocated( adjpntr ) ) deallocate( adjpntr )
        if( allocated( adjcol  ) ) deallocate( adjcol  )

        ncomp = maxval( comp )
        print*, 'Number of connected components =', ncomp
        print*, ''

        !Generate list 'compsize' assigning size of CC, bond number of CC and a representative state
        if( allocated( compsize ) ) deallocate( compsize )
        if( allocated( adj_coo ) ) deallocate( adj_coo )
        allocate( compsize( ncomp, 3 ) ) !1st col: Size of CC; 2nd col: Bond number of CC; 3rd col: Representative of CC; Row # = CC #
        allocate( adj_coo( 2 * dim * sites, 2 ) )
        compsize = 0
        do i = 1, dim
            if ( compsize( comp(i), 1 ) == 0 ) then
                compsize( comp(i), 3 ) = i
                do j = 0, sites - 1
                    if ( k == 1 .and. btest( basis(i), j ) .and. btest( basis(i), modulo(j + 1, sites))) compsize(comp(i), 2) = compsize(comp(i), 2) + 1
                    if ( k == 2 .and. btest( basis(i), j ) .and. btest( basis(i), modulo(j + 2, sites))) compsize(comp(i), 2) = compsize(comp(i), 2) + 1
                end do
            end if
            compsize( comp(i), 1 ) = compsize( comp(i), 1 ) + 1
        end do
        cntr = 0
        do i = 1, ncomp - 1
            do j = i + 1, ncomp
                if ( compsize( i, 1 ) == compsize( j, 1 ) .and. compsize( i, 2 ) == compsize( j, 2 ) ) then
                    cntr = cntr + 1
                    adj_coo( cntr, 1) = compsize( i, 3 )
                    adj_coo( cntr, 2) = compsize( j, 3 )
                    cntr = cntr + 1
                    adj_coo( cntr, 1) = compsize( j, 3 )
                    adj_coo( cntr, 2) = compsize( i, 3 )
                end if
            end do
        end do

        if( allocated( adj ) ) deallocate( adj )
        allocate( adj( nnz + cntr, 2 ) )
        do i = 1, nnz
            adj( i, 1 ) = adj_temp( i, 1 )
            adj( i, 2 ) = adj_temp( i, 2 )
        end do
        do i = nnz + 1, nnz + cntr
            adj( i, 1 ) = adj_coo( i - nnz, 1 )
            adj( i, 2 ) = adj_coo( i - nnz, 2 )
        end do
        nnz = nnz + cntr

        if( allocated( adj_temp ) ) deallocate( adj_temp )
        if( allocated( adj_coo ) ) deallocate( adj_coo )
        if( allocated( adjpntr ) ) deallocate( adjpntr )
        if( allocated( adjcol  ) ) deallocate( adjcol  )
        if( allocated( dummy   ) ) deallocate( dummy   )

        allocate( adjpntr( dim + 1 ) )
        allocate( adjcol( nnz ) )
        allocate( dummy( nnz ) )

        call coocsr( dim, nnz, dummy, adj(1:nnz,1), adj(1:nnz,2), dummy, adjcol, adjpntr ) !Creates CSR-format (ao,jao,iao) sparse matrix from coordinate format (a,ir,jc)

        if( allocated( dummy   ) ) deallocate( dummy   )

        call connected_components( nnz, dim, adjcol, adjpntr, comp )
        if( allocated( adjpntr ) ) deallocate( adjpntr )
        if( allocated( adjcol  ) ) deallocate( adjcol  )

        ncomp = maxval( comp )

        if ( k == 1 ) then
            ncompv1 = ncomp
            compv1  = comp
        else if ( k == 2 ) then
            ncompv2 = ncomp
            compv2  = comp
        end if
    end do
    print*, 'V1: Number of connected components =', ncompv1
    print*, ''
    print*, 'V2: Number of connected components =', ncompv2
    print*, ''
    return
end subroutine scc



!-------------------------------------------------------------!
!            Calculate connected components (dense)           !
!-------------------------------------------------------------!

    !if (ip == 1) then
    !    do k = 1, 2
    !        if ( not( allocated( adj ) ) )   allocate( adj( dim_hs,dim_hs ) )
    !        if ( not( allocated( comp ) ) )  allocate( comp( dim_hs ) )
    !        if ( not( allocated( dad ) ) )   allocate( dad( dim_hs ) )
    !        if ( not( allocated( order ) ) ) allocate( order( dim_hs ) )
    !        if ( not( allocated( periods ) ) ) allocate( periods( dim_hs ) )
    !        adj     = 0
    !        comp    = 0
    !        dad     = 0
    !        order   = 0
    !        periods = 0
    !        !Generate first adjacency matrix
    !        call slhopping_serial ( dim_hs, sites, bc, 0, k, eps, t, permutations, pseudo_rc, pseudo_ham, pseudoi, mom, periods, adj )
    !        call digraph_adj_components ( adj, dim_hs, dim_hs, ncomp, comp, dad, order )
    !        deallocate( periods )
    !        print*, 'Number of connected components =', ncomp
    !        print*, ''
    !
    !        !Generate list 'compsize' assigning size of CC, bond number of CC and a representative state
    !        if( allocated( compsize ) ) deallocate( compsize )
    !        allocate( compsize( ncomp, 3 ) ) !1st col: Size of CC; 2nd col: Bond number of CC; 3rd col: Representative of CC; Row # = CC #
    !        compsize = 0
    !        do i = 1, dim_hs
    !            if ( compsize( comp(i), 1 ) == 0 ) then
    !                compsize( comp(i), 3 ) = i
    !                do j = 0, sites - 1
    !                    if ( k == 1 .and. btest( permutations(i), j ) .and. btest( permutations(i), modulo(j + 1, sites))) compsize(comp(i), 2) = compsize(comp(i), 2) + 1
    !                    if ( k == 2 .and. btest( permutations(i), j ) .and. btest( permutations(i), modulo(j + 2, sites))) compsize(comp(i), 2) = compsize(comp(i), 2) + 1
    !                end do
    !            end if
    !            compsize( comp(i), 1 ) = compsize( comp(i), 1 ) + 1
    !        end do
    !        do i = 1, ncomp - 1
    !            do j = i + 1, ncomp
    !                if ( compsize( i, 1 ) == compsize( j, 1 ) .and. compsize( i, 2 ) == compsize( j, 2 ) ) then
    !                    adj( compsize( i, 3 ), compsize( j, 3 ) ) = 1
    !                    adj( compsize( j, 3 ), compsize( i, 3 ) ) = 1
    !                end if
    !            end do
    !        end do
    !        ncomp = 0
    !        comp  = 0
    !        dad   = 0
    !        order = 0
    !        !Calculate second adjacency matrix
    !        call digraph_adj_components ( adj, dim_hs, dim_hs, ncomp, comp, dad, order )
    !        deallocate( dad, order, pseudo_ham, pseudo_rc  )
    !        if ( k == 1 ) then
    !            ncompv1 = ncomp
    !            compv1 = comp
    !        else if ( k == 2 ) then
    !            ncompv2 = ncomp
    !            compv2 = comp
    !        end if
    !    end do
    !end if






!------------------------------------------------------!
!            Spinless: Hopping Hamiltonian             !
!------------------------------------------------------!

subroutine slhopping(dim, sites, ti, bc, eps, t, basis, rc, ham, nnz, mom, periods)

    implicit none

    integer(kind=8), intent(in) :: dim, basis(dim)
    integer, intent(in) :: sites, ti
    integer, intent(in), optional :: mom
    integer(kind=8), intent(in), optional :: periods(dim)
    double precision, intent(in) :: t, eps
    character, intent(in) :: bc

    integer, allocatable, intent(out) :: rc(:,:)
    double complex, allocatable, intent(out) :: ham(:)
    integer, intent(inout) :: nnz


    double complex, parameter :: ii = (0, 1)
    integer :: pntr = 0
    integer :: mask = 0
    integer(kind=8) :: loc = 0
    integer(kind=8) :: ntrans = 0
    integer(kind=8) :: rep = 0
    integer, allocatable :: n_temp(:)
    integer :: chunk = 0
    integer :: remainder = 0
    integer :: arrsize = 0
    integer :: n_threads = 0
    integer :: my_stat = 0
    integer :: thread_num = 0
    integer :: level = 0
    integer :: team_size = 0
    integer :: i = 0 , j = 0 , m = 0 , k = 0 , l = 0 , q = 0
    integer :: cntr = 0 , cntrj = 0 , dbl = 0
    integer, allocatable :: rc_temp(:,:)
    integer, allocatable :: temp_rc(:,:,:)
    double complex, allocatable :: ham_temp(:)
    double complex, allocatable:: temp(:,:)

    print*, 'Building hopping Hamiltonian ...'
    print*, ''

    !$omp parallel
        !$ n_threads= omp_get_num_threads()
    !$omp end parallel

!    if(.not. dynamic) n_threads = max(int((n_threads - 2)/2),1)
    n_threads = max(int((n_threads - 2)/2),1)

    n_temp = 0
    !$omp parallel default(firstprivate) shared(basis, periods, temp, temp_rc, n_temp, n_threads, chunk, remainder, arrsize) num_threads(n_threads)
    !$ thread_num = omp_get_thread_num()


    !$omp single
        !$ n_threads = omp_get_num_threads()
        chunk = int((dim)/n_threads)
        remainder = mod((dim),n_threads)
        arrsize = (dim - (n_threads-1)*chunk + remainder) * sites
        if (allocated(temp)) deallocate(temp)
        if (allocated(temp_rc)) deallocate(temp_rc)
        if (allocated(n_temp)) deallocate(n_temp)
        allocate(temp(arrsize,n_threads))
        allocate(temp_rc(arrsize,n_threads,2))
        allocate(n_temp(n_threads))
        !$ level = omp_get_level()
        n_temp = 0

    !$omp end single

    !$omp do schedule(static, chunk)
    do j = 1, size(basis)
        !$ thread_num = omp_get_thread_num()

        cntrj = n_temp(thread_num + 1)
        cntr = 0
        do m = 0, sites - 1 ! m goes through all digits of basis states
            if (m == sites-1 .and. bc /= 'p') cycle

            mask=ibset(ibset(0,m),mod(m+1,sites))  !Procedure suggested in Lin-paper. Sets a mask with only non-zero components on sites i=m+1 and i=m+2
            k=iand(mask,basis(j)) !K records the occupancy of those two sites with up-spins
            l=ieor(mask,k) !L records whether hopping is allowed or not. If it's 0 or the same as the mask, hopping is not allowed. If it is allowed the occupations of both sites (01 or 10) are swapped (10 or 01)
            if (l == 0 .or. l == mask) then
                cycle
            end if

            if(ti == 1 .and. bc == 'p') then
                call representative(basis(j) - k + l, sites, rep, ntrans)
            else
                rep = basis(j) - k + l
                ntrans = 0
            end if
            call findstate(dim, rep, basis, loc)
            if (loc > 0) then
                dbl = 0
                cntr = cntr + 1
                if (cntrj > 0 .and. ti == 1) then
                    do i = cntrj + 1, cntrj + 1 + cntr
                        if(temp_rc(i,thread_num + 1,2) == loc) then
                            if (m == sites-1 .and. bc == 'p' ) then
                                temp(i,thread_num + 1) = temp(i,thread_num + 1) + (-1)**(popcnt(basis(j))-1)*(-1)*t*eps* (sqrt(real(periods(j))/real(periods(loc))) * exp(ii*2*pi*mom*ntrans/sites))**ti
                            else
                                temp(i,thread_num + 1) = temp(i,thread_num + 1) + (-1)*t*eps* (sqrt(real(periods(j))/real(periods(loc))) * exp(ii*2*pi*mom*ntrans/sites))**ti
                            end if
                            dbl = 1
                        end if
                    end do
                end if
                if (dbl == 0) then
                    n_temp(thread_num + 1) = n_temp(thread_num + 1) + 1
                    if (m == sites-1 .and. bc == 'p' ) then
                        temp(n_temp(thread_num + 1),thread_num + 1) = temp(n_temp(thread_num + 1),thread_num + 1) + (-1)**(popcnt(basis(j))-1)*(-1)*t*eps* (sqrt(real(periods(j))/real(periods(loc))) * exp(ii*2*pi*mom*ntrans/sites))**ti
                    else
                        temp(n_temp(thread_num + 1),thread_num + 1) = temp(n_temp(thread_num + 1),thread_num + 1) + (-1)*t*eps* (sqrt(real(periods(j))/real(periods(loc))) * exp(ii*2*pi*mom*ntrans/sites))**ti
                    end if
                    temp_rc(n_temp(thread_num + 1),thread_num + 1,1) = j
                    temp_rc(n_temp(thread_num + 1),thread_num + 1,2) = loc
                end if
            end if
        end do
    end do
    !$omp end do
    !$omp end parallel



    nnz = sum(n_temp)
    if (allocated(ham_temp)) deallocate(ham_temp)
    if (allocated(rc_temp)) deallocate(rc_temp)
    allocate(rc_temp(nnz,2))
    allocate(ham_temp(nnz))

    pntr = 1
    do i = 1, n_threads
        if (i > n_threads) exit
        ham_temp(pntr:pntr + n_temp(i) - 1) = temp(1:n_temp(i),i)
        rc_temp(pntr:pntr + n_temp(i) - 1,1) = temp_rc(1:n_temp(i),i,1)
        rc_temp(pntr:pntr + n_temp(i) - 1,2) = temp_rc(1:n_temp(i),i,2)
        pntr = pntr + n_temp(i)
    end do



    call move_alloc(ham_temp,ham)
    call move_alloc(rc_temp,rc)


    if (allocated(ham_temp)) deallocate(ham_temp)
    if (allocated(rc_temp)) deallocate(rc_temp)
    if (allocated(temp)) deallocate(temp)
    if (allocated(temp_rc)) deallocate(temp_rc)
    if (allocated(n_temp)) deallocate(n_temp)

    print*, 'Finished hopping Hamiltonian.'
    print*, ''

end subroutine slhopping





!----------------------------------------------------------!
!               Spinless: Diagonal Hamiltonian             !
!----------------------------------------------------------!

subroutine slhamdi (sites, dim, rc, ham, basis, occ, ndi)

    implicit none

    integer(kind=8), intent(in) :: dim, basis(dim)
    integer, intent(in) :: sites
    integer, intent(out) :: ndi
    integer, allocatable, intent(out) :: occ(:,:), rc(:)
    double precision, allocatable, intent(out) :: ham(:,:)

    integer :: counter_m_odd = 0, counter_m_even = 0, counter_v(sites), counter_v2(sites)
    integer :: pntr = 0
    integer :: chunk = 0
    integer :: remainder = 0
    integer :: arrsize = 0
    integer :: n_threads = 0
    integer :: thread_num = 0
    integer :: level = 0
    integer :: team_size = 0
    integer :: j = 0, m = 0, s = 0

    integer, allocatable :: rc_temp(:)
    double precision, allocatable :: ham_temp(:,:)


    print*, 'Building diagonal Hamiltonian ...'
    print*, ''

    !$omp parallel
        !$ n_threads= omp_get_num_threads()
    !$omp end parallel
    !if(.not. dynamic) n_threads = max(int((n_threads - 2)/2),1)
    n_threads = max(int((n_threads - 2)/2),1)


    if(allocated(occ)) deallocate(occ)
    allocate(occ(sites, dim))

    !$omp parallel default(firstprivate) shared(occ, rc_temp, ham_temp, n_threads) num_threads(n_threads)
        !$omp single
            !$ n_threads = omp_get_num_threads()
            if (allocated(ham_temp)) deallocate(ham_temp)
            if (allocated(rc_temp)) deallocate(rc_temp)
            allocate(ham_temp(dim,2))
            allocate(rc_temp(dim))
            !$ level = omp_get_level()
            !$ team_size = omp_get_team_size(level)
        !$omp end single


        !$omp do schedule(static,chunk)
        do j= 1, dim
            !$ thread_num = omp_get_thread_num()
            counter_v = 0
            counter_v2 = 0
            do m = 0, sites - 1 ! m goes through all digits of each configuration

                if (btest(basis(j),m) .and. btest(basis(j), modulo(m + 1, sites))) counter_v(m + 1) = counter_v(m + 1) + 1
                if (btest(basis(j),m) .and. btest(basis(j), modulo(m + 2, sites))) counter_v2(m + 1) = counter_v2(m + 1) + 1

                !---------------------------------!
                !             Disorder            !
                !---------------------------------!

                if (btest(basis(j),m)) occ(m+1,j) = occ(m+1,j) + 1
            end do
            ham_temp(j,1) = real(sum(counter_v))
            ham_temp(j,2) = real(sum(counter_v2))
            rc_temp(j) = j
        end do
        !$omp end do
    !$omp end parallel

    ndi = dim
    call move_alloc(rc_temp, rc)
    call move_alloc(ham_temp, ham)
    if (allocated(ham_temp)) deallocate(ham_temp)
    if (allocated(rc_temp)) deallocate(rc_temp)

end subroutine slhamdi




!-----------------------------------------------------------!
!            Spinless: Unify Hamiltonian  (complex)         !
!-----------------------------------------------------------!

subroutine slchamunify (eps, v1, v2, w, sites, occ, noff, ndi, hamoff, rcoff, hamdi, rcdi, ham, rc, nnz)

    implicit none

    double precision, intent(in) :: eps, v1, v2, w
    integer, intent(in):: sites, noff, ndi, occ(sites,*)
    double complex, intent(in):: hamoff(noff)
    double precision, intent(in)::  hamdi(ndi,2)
    integer, intent(in):: rcoff(noff,2), rcdi(ndi)
    double complex, allocatable, intent(out):: ham(:)
    integer, allocatable, intent(out):: rc(:,:)
    integer, intent(out)::  nnz

    integer :: j = 0
    integer :: temprcdi(ndi)
    double precision :: tempdi(ndi), rand(sites)

    call random_number(rand)
    rand = 2*(rand-0.5)


    nnz = noff + ndi !Number of non-zero elements
    do j = 1, ndi
        tempdi(j) = v1 * hamdi(j,1) + v2 * hamdi(j,2) + eps * w * sum(rand * occ(1:sites, j)) ! - mu * pts
        temprcdi(j) = rcdi(j)
    end do

    if( allocated( ham ) )  deallocate( ham )
    if( allocated( rc ) )   deallocate( rc )
    allocate( ham( nnz ) )
    allocate( rc( nnz, 2 ) )

    ham = 0
    rc  = 0

    ham( 1:noff )         = hamoff( 1:noff )
    ham( noff + 1:nnz )   = tempdi( 1:ndi )
    rc( 1:noff, 1 )       = rcoff( 1:noff, 1 )
    rc( noff + 1:nnz, 1 ) = temprcdi( 1:ndi )
    rc( 1:noff, 2 )       = rcoff( 1:noff, 2 )
    rc( noff+1:nnz, 2 )   = temprcdi( 1:ndi )


end subroutine slchamunify






!------------------------------------------------------------------!
!            Spinless: Fill dense Hamiltonian (complex)            !
!------------------------------------------------------------------!

subroutine slchamdense (dim, eps, v1, v2, w, sites, occ, noff, ndi, hamoff, rcoff, hamdi, rcdi, ham, nnz)

    implicit none

    integer(kind=8), intent(in) :: dim
    double precision, intent(in) :: eps, v1, v2, w
    integer, intent(in):: sites, noff, ndi, occ(sites,*)
    double complex, intent(in):: hamoff(noff)
    double precision, intent(in)::  hamdi(ndi,2)
    integer, intent(in):: rcoff(noff,2), rcdi(ndi)

    double complex, allocatable, intent(out):: ham(:,:)
    integer, intent(out)::  nnz

    integer :: temprcdi(ndi)
    double precision :: tempdi(ndi), rand(sites)
    integer :: jj = 0

    !$OMP THREADPRIVATE (jj)
    jj = 0
    nnz = 0

    call random_number(rand)
    rand = 2*(rand-0.5)

    if(allocated(ham)) deallocate(ham)
    allocate(ham(dim,dim))
    ham = 0
    if(noff > 0) then
        do jj = 1, noff
            !if(rcoff(jj,2) > rcoff(jj,1))
            ham(rcoff(jj,1),rcoff(jj,2)) = hamoff(jj)
        end do
    end if

    jj = 0
    if(ndi > 0) then
        do jj = 1, ndi
            ham(rcdi(jj),rcdi(jj)) = v1 * hamdi(jj,1) + v2 * hamdi(jj,2) + eps * w * sum(rand * occ(1:sites, jj)) ! - mu * pts
        end do
    end if

    nnz = noff + ndi !Number of non-zero elements

return
end subroutine slchamdense





!---------------------------------------------!
!            Factorization of basis           !
!---------------------------------------------!

subroutine factorize(dim, lin, lr, b)

    implicit none
    integer(kind=8), intent(in) :: dim, lin(dim), lr

    integer(kind=8), allocatable, intent(out) :: b(:,:)

    integer :: j = 0

    if (allocated(b)) deallocate(b)
    allocate(b(dim, 2))
    b = 0


    !!$omp parallel do num_threads(num_threads)
    do j = 1, dim
        b(j, 1) = lin(j)/(2**lr) !left subsystem
        b(j, 2) = mod(lin(j),(2**lr)) !right subsystem
    end do
    !!$omp end parallel do



end subroutine factorize




!---------------------------------------------------!
!            Complex: Entanglement entropy          !
!---------------------------------------------------!

subroutine centent(threads, dim, lr, ll, b, state, eps, entropy, singval)

    implicit none

    integer, intent(in) :: threads
    integer(kind=8), intent(in) :: dim, lr, ll, b(dim,2)
    double complex, intent(in) :: state(dim)
    double precision, intent(in) :: eps
    double precision, intent(out) :: entropy
    double precision, allocatable , intent(out):: singval(:)


    double complex :: c(2**ll, 2**lr)
    double complex :: cd(2**lr, 2**ll)
    double complex :: ccd(2**lr, 2**lr)
    integer :: j = 0


    if (allocated(singval)) deallocate(singval)
    c = 0
    cd = 0
    ccd = 0

    !!$omp parallel do num_threads(threads)
    do j = 1, dim
        c(b(j,1) + 1, b(j,2) + 1) = state(j)
        cd(b(j,2) + 1, b(j,1) + 1) = conjg(state(j))
    end do
    !!$omp end parallel do


    ccd = matmul(cd,c)


    !!$omp critical
!    call cfulldiag(.False., 2**lr, ccd, singval)
    call cfulldiag2(.False., 2**lr, ccd, singval)
    !!$omp end critical
    entropy = 0.0d0

    !!$omp parallel do num_threads(threads)
    do j = 1, 2**lr
        if (singval(j) > eps) then
            entropy = entropy - singval(j) * dlog(singval(j))
        end if
    end do
    !!$omp end parallel do

end subroutine centent








!-----------------------------------------------------------------------!
!            Complex: Entanglement entropy for momentum states          !
!-----------------------------------------------------------------------!

subroutine ti_centent(threads, dim, dim2, lr, ll, base, state, eps, k, sites, periods, transl, locs, entropy)

    implicit none

    integer, intent(in) :: threads, k, sites
    integer(kind=8), intent(in) :: dim, dim2, lr, ll, base(dim2), periods(dim2), transl(dim2), locs(dim2)
    double complex, intent(in) :: state(dim)
    double precision, intent(in) :: eps
    double precision, intent(out) :: entropy

    double precision, allocatable, save :: singval(:)
    double complex, allocatable, save :: c(:,:)
    double complex, allocatable, save :: cd(:,:)
    double complex, allocatable, save :: ccd(:,:)
    double complex :: phase
    integer, save :: loop = 0, loop2 = 0
    logical, save :: check
    double precision, save :: delta  = 0.0000001d0

    !$OMP THREADPRIVATE (check, loop, loop2, singval, c, cd, ccd, delta)
    loop  = 0
    loop2 = 0

    if (allocated(singval)) deallocate(singval)
    if (allocated(c)) deallocate(c)
    if (allocated(cd)) deallocate(cd)
    if (allocated(ccd)) deallocate(ccd)

    allocate(singval(dim))
    allocate(c(2**ll, 2**lr))
    allocate(cd(2**lr, 2**ll))
    allocate(ccd(2**lr, 2**lr))

    c = (0.d0, 0.d0)
    cd = (0.d0, 0.d0)

    loop = 0

    !!$omp critical
    do loop = 1, size(base)
        if (periods(loop) < 0 .or. locs(loop) < 1) cycle
        c(base(loop)/(2**lr) + 1, mod(base(loop),(2**lr)) + 1) = state(locs(loop)) * sqrt(dble(periods(loop)))/dble(sites) * exp(-2*pi*ii*k*transl(loop)/sites)
        cd(mod(base(loop),(2**lr)) + 1, base(loop)/(2**lr) + 1) = dconjg(c(base(loop)/(2**lr) + 1, mod(base(loop),(2**lr)) + 1))
    end do
    !!$omp end critical

    ccd = (0.d0, 0.d0)
    ccd = matmul(cd,c)

    call testhermitian (check, 2**lr, ccd, delta)
    call cfulldiag(.False., 2**lr, ccd, singval)
    !call cfulldiag2(.False., 2**lr, ccd, singval)
    entropy = 0.0d0

    !!$omp parallel do num_threads(threads)
    do loop2 = 1, 2**lr
        if (singval(loop2) > eps) then
            !!$omp atomic
            entropy = entropy - singval(loop2) * dlog(singval(loop2))
        end if
    end do
    !!$omp end parallel do

    if (allocated(singval)) deallocate(singval)

return
end subroutine ti_centent



!----------------------------------------------!
!            Level spacing parameter           !
!----------------------------------------------!

subroutine lvlpar (frac, nev, evals, rmean)

    implicit none

    integer, intent(in) :: frac, nev
    double complex, intent(in) :: evals(nev)
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
        diff(cntr) = real(evals(i+1)) - real(evals(i))
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



!-------------------------------------------------------------------------!
!            Calculate number of eigenvalues and Lanczos vectors          !
!-------------------------------------------------------------------------!

subroutine nevncv (thresh, exact, nevext, nestext, ncv0, dim, full, nev, ncv, nest)

    implicit none
    integer, intent(in) :: thresh, exact, nevext, nestext, ncv0
    integer(kind=8), intent(in) :: dim
    integer, intent(out) :: full, nev, ncv, nest

    if (exact == 0) then
        if (dim < thresh) then
            full = 1
            nev = dim
            ncv = dim
            nest = dim
        else
            full = 0
            if (dim == 1) then
                nev = 1
                ncv = 1
                nest = 1
            else
                nev = min(nevext, dim - 10)
                ncv = max(ncv0, 2*nev + 10) !nev+1<=ncv<=dim and ncv maximal
                if (ncv > dim) ncv = dim
                nest = min(nev, nestext)
            end if
        end if
    else if (exact == 1) then
        full = 1
        nev = dim
        ncv = dim
        nest = min(dim, nestext)
    end if
    print*, 'Number of eigenvalues = ', nev
    print*, ''
    print*, 'Number of eigenvectors = ', nest
    print*, ''

end subroutine nevncv




!---------------------------------------------------------!
!            Complex: Full exact diagonalization          !
!---------------------------------------------------------!

subroutine cfulldiag(eigvec, dim, mat, evals)

    implicit none
    logical, intent(in) :: eigvec
    integer(kind=8), intent(in) :: dim
    double complex, intent(inout) :: mat(dim, dim)
    double precision, allocatable, intent(out) :: evals(:)

    double complex, allocatable, save:: matloc(:,:)
    double precision, allocatable, save :: evalsloc(:)
    !integer(kind=8) :: lda = 0
    integer, save :: lda = 0
    integer, parameter :: lwmax = 100000
    integer, save :: i = 0, info = 0, lwork = 0
    character, save :: jobz
    double precision, allocatable, save :: rwork(:)
    complex*16, allocatable, save :: work(:)

    external :: zheev, zheev_2stage

    !$OMP THREADPRIVATE (matloc, evalsloc, lda, i, info, lwork, jobz, rwork, work)
    lda   = 0
    info  = 0
    lwork = 0

    if (eigvec) then
        jobz = 'V'
    else
        jobz = 'N'
    end if

    if (allocated(evals)) deallocate(evals)
    allocate(evals(dim))
    evals = 0.0d0

    if (allocated(evalsloc)) deallocate(evalsloc)
    if (allocated(matloc)) deallocate(matloc)
    allocate(evalsloc(dim))
    allocate(matloc(dim, dim))
    evalsloc = 0.0d0
    matloc = mat

    if (dim > 1) then
        lda = dim
        if (allocated(rwork)) deallocate(rwork)
        allocate(rwork(max(1,3*dim - 2)))
        rwork = 0.0d0

        if (allocated(work)) deallocate(work)
        allocate(work(lwmax))
        lwork = -1
        if (jobz == 'V') then
            !$omp critical
            call zheev (jobz, 'U', dim, matloc, dim, evalsloc, work, lwork, rwork, info)
            !$omp end critical
        else if (jobz == 'N') then
            call zheev_2stage ('N', 'U', dim, matloc, lda, evalsloc, work, lwork, rwork, info)
        end if
        !lwork = min(int(work(1)), lwmax)
        lwork = int(work(1))
        if (lwork < max(1,2*dim - 1)) lwork = max(1,2*dim - 1)
        if (allocated(work)) deallocate(work)
        allocate(work(lwork))
        work = 0.0d0

        !print*, 'Running diagoalization routine ...'
        if (jobz == 'V') then
            !$omp critical
            call zheev (jobz, 'U', dim, matloc, dim, evalsloc, work, lwork, rwork, info)
            !$omp end critical
        else if (jobz == 'N') then
            call zheev_2stage ('N', 'U', dim, matloc, lda, evalsloc, work, lwork, rwork, info)
        end if

        !print*, 'Finished diagoalization routine.'
        if (info .gt. 0 ) THEN
            write(*,*)'The algorithm failed to compute eigenvalues.'
            stop
        end if
        !print*, 'Eigenvalues:'
        !do i = 1, size(evalsloc)
        !    print*, evalsloc(i)
        !end do
        if (allocated(work)) deallocate(work)
        if (allocated(rwork)) deallocate(rwork)
    else
        evalsloc = real(matloc(1,1))
    end if

    evals = evalsloc
    mat = matloc
    if(allocated(evalsloc)) deallocate(evalsloc)
    if(allocated(matloc)) deallocate(matloc)

    return
end subroutine cfulldiag




!--------------------------------------------------------------!
!            Complex: Full exact diagonalization (NEW)         !
!--------------------------------------------------------------!

subroutine cfulldiag2(eigvec, dim, mat, evals)

    implicit none
    logical, intent(in) :: eigvec
    integer(kind=8), intent(in) :: dim
    double complex, intent(inout) :: mat(dim, dim)
    double precision, allocatable, intent(out) :: evals(:)

    integer, parameter :: lwmax = 100000
    integer, save :: dimloc = 0
    integer, save :: ldz = 0
    integer, save :: lda = 0
    integer, save :: m = 0
    integer, save :: i = 0, info = 0, lwork = 0, lrwork = 0, liwork = 0, il = 0, iu = 0
    character, save :: jobz
    double precision, save :: abstol, vl, vu

    integer, allocatable, save:: isuppz(:)
    integer, allocatable, save:: iwork(:)
    double precision, allocatable, save :: rwork(:)
    double complex, allocatable, save :: work(:)
    double complex, allocatable, save :: z(:,:)
    double complex, allocatable, save :: matloc(:, :)
    double precision, allocatable, save :: evalsloc(:)
    !     .. External Subroutines ..
    external :: zheevr, zheevr_2stage
    !external :: PRINT_MATRIX, PRINT_RMATRIX

    !     .. Intrinsic Functions ..
    intrinsic :: int, min

    !$OMP THREADPRIVATE (i, dimloc, lda, ldz, info, lwork, m, lrwork, liwork, il, iu, &
    !$OMP                abstol, vl, vu, &
    !$OMP                jobz, rwork, work, isuppz, iwork, z, matloc, evalsloc)

    dimloc = dim
    info  = 0
    lwork = 0
    m     = 0
    lda   = dim
    ldz   = dim

    if (eigvec) then
        jobz = 'V'
    else
        jobz = 'N'
    end if

    if (allocated(evals)) deallocate(evals)
    allocate(evals(dimloc))
    if (allocated(evalsloc)) deallocate(evalsloc)
    if (allocated(matloc)) deallocate(matloc)
    allocate(evalsloc(dimloc))
    allocate(matloc(dimloc,dimloc))
    evals    = 0
    evalsloc = 0
    matloc   = mat

    if (dimloc > 1) then

        !Negative ABSTOL means using the default value
        abstol = -1.0

        !Set VL, VU to compute eigenvalues in half-open (VL,VU] interval
        vl     = -5.0
        vu     = 5.0

        !Query the optimal workspace.
        lwork  = -1
        lrwork = -1
        liwork = -1

        if (allocated(work)) deallocate(work)
        if (allocated(iwork)) deallocate(iwork)
        if (allocated(isuppz)) deallocate(isuppz)
        if (allocated(z)) deallocate(z)
        if (allocated(rwork)) deallocate(rwork)

        allocate(work(lwmax))
        allocate(rwork(lwmax))
        allocate(z(dimloc,dimloc))
        allocate(iwork(lwmax))
        allocate(isuppz(dimloc))

        if (jobz == 'V') then
            call zheevr( 'Vectors', 'A', 'Upper', dimloc, matloc, lda, vl, vu, il,&
                         iu, abstol, m, evalsloc, z, ldz, isuppz, work, lwork, rwork,&
                         lrwork, iwork, liwork, info )

        else if (jobz == 'N') then
            call zheevr_2stage( 'N', 'A', 'Upper', dimloc, matloc, lda, vl, vu, il,&
                         iu, abstol, m, evalsloc, z, ldz, isuppz, work, lwork, rwork,&
                         lrwork, iwork, liwork, info )

        end if

        !lwork = min( lwmax, int( work(1) ) )
        !lrwork = min( lwmax, int( rwork(1) ) )
        !liwork = min( lwmax, iwork(1) )
        lwork  = int( work(1)  )
        lrwork = int( rwork(1) )
        liwork = iwork(1)

        !     Solve eigenproblem.

        if (allocated(work)) deallocate(work)
        if (allocated(iwork)) deallocate(iwork)
        if (allocated(isuppz)) deallocate(isuppz)
        if (allocated(rwork)) deallocate(rwork)

        allocate(work(lwork)) !lwmax
        allocate(rwork(lrwork)) !lwmax
        allocate(iwork(liwork)) !lwmax
        allocate(isuppz(dimloc))



        if (jobz == 'V') then
            call zheevr( 'Vectors', 'A', 'Upper', dimloc, matloc, lda, vl, vu, il,&
                         iu, abstol, m, evalsloc, z, ldz, isuppz, work, lwork, rwork,&
                         lrwork, iwork, liwork, info )
        else if (jobz == 'N') then
            call zheevr_2stage( 'N', 'A', 'Upper', dimloc, matloc, lda, vl, vu, il,&
                         iu, abstol, m, evalsloc, z, ldz, isuppz, work, lwork, rwork,&
                         lrwork, iwork, liwork, info )
        end if
        !     Check for convergence.

        evals = evalsloc
        mat   = z

        if (allocated(work)) deallocate(work)
        if (allocated(iwork)) deallocate(iwork)
        if (allocated(isuppz)) deallocate(isuppz)
        if (allocated(z)) deallocate(z)
        if (allocated(rwork)) deallocate(rwork)

        if( info .gt. 0 ) then
            write(*,*)'The algorithm failed to compute eigenvalues.'
            stop
        end if

        !   Print the number of eigenvalues found.

        !write(*,'(/a,i2)')' The total number of eigenvalues found:', m
        !print*, m, 'Number of eigenvalues found'


        !print*, 'Eigenvalues:'
        !do i = 1, size(evals)
        !    print*, evals(i)
        !end do
    else
        evals = real(mat(1,1))
    end if


    if (allocated(evalsloc)) deallocate(evalsloc)
    if (allocated(matloc)) deallocate(matloc)

return
end subroutine cfulldiag2



!----------------------------------------------------!
!            Parallel complex Spmmv: y = A*x         !
!----------------------------------------------------!

subroutine cpamux (threads, n, x, y, a, ja, ia) !My parallelized sparse matrix-vector multiplication
implicit none
!include "omp_lib.h"
double complex  :: x(*), y(*), a(*)
integer(kind=8) :: n
integer :: ja(*), ia(*)
integer :: threads

!-----------------------------------------------------------------------
!         A times a vector
!-----------------------------------------------------------------------
! multiplies a matrix by a vector using the dot product form
! Matrix A is stored in compressed sparse row storage.
!
! on entry:
!----------
! n     = row dimension of A
! x     = real array of length equal to the column dimension of
!         the A matrix.
! a, ja,
!    ia = input matrix in compressed sparse row format.
!
! on return:
!-----------
! y     = real array of length n, containing the product y=Ax
!
!-----------------------------------------------------------------------
! local variables
!
      double complex :: t
      integer :: i, k
!-----------------------------------------------------------------------


      !$omp parallel do private(k,t) num_threads(threads)
      do 100 i = 1,n

!
!     compute the inner product of row i with vector x
!

         t = 0.0d0
         do 99 k = ia(i), ia(i+1)-1
            t = t + a(k)*x(ja(k))
 99      continue

!
!     store result in y(i)
!
         y(i) = t
 100  continue
      !$omp end parallel do


      return
end subroutine




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



!-------------------------------------------!
!            Number of basis states         !
!-------------------------------------------!

integer function  combs(nsites, part, magn) !Calculates the number of possible combinations for a given magnetization and particle number!

    implicit none

    integer, intent(in) :: nsites, part
    double precision, intent(in) :: magn
    integer :: part_up, part_dn


    part_up=int((part+2*magn)/2)
    part_dn=part-part_up

    if (part_up < 0 .or. part_dn < 0) then
        combs = 0
    else
        combs = int(fact(nsites)/(fact(part_up)*fact(max(1,nsites-part_up))),8)*int(fact(nsites)/(fact(part_dn)*fact(max(1,nsites-part_dn))),8)
    end if

end function combs



!----------------------------------!
!            Factorial: n!         !
!----------------------------------!

double precision function fact(nn) !Calculates the factorial n!

    implicit none

    integer, intent(in) :: nn
    integer :: i

    if (nn < 0) error stop 'factorial is singular for negative integers'

    if(nn==0) then
        fact=1
    else
        fact = 1
        do i = 2, nn
            fact = fact * i
        end do
    end if

end function fact


!---------------------------------------!
!            Binomial coefficients      !
!---------------------------------------!


double precision function binomial(n, k)

    implicit none
    integer, intent(in) :: n, k

    binomial = fact(n) / (fact(k) * fact(n - k))

end function binomial





!--------------------------------------------!
!            Complex diagonalization         !
!--------------------------------------------!

subroutine complex_diag(threads, dim, nev, ncv, nst, mode, rvec, nnz, ham, rc, evals, evecs)

    implicit none

    integer, intent(in) :: threads, nev, ncv, nst
    integer(kind=8), intent(in) :: dim
    character*2, intent(in) :: mode
    logical, intent(in) :: rvec
    integer, intent(in) :: nnz
    double complex, intent(in) :: ham(nnz)
    integer, intent(in) :: rc(nnz,2)
    double complex, intent(out) :: evals(nev), evecs(dim, nst)

    double complex, allocatable:: ao(:)
    integer, allocatable:: jao(:), iao(:)

    integer(kind=8) :: maxn, maxnev, maxncv, ldv
    integer(kind=8) :: n, nx
    intrinsic :: abs
    double precision :: tol
    character :: bmat*1, which*2
    double complex :: sigma
    integer :: j, iparam(11), ipntr(14)
    integer :: ido, ishfts, lworkl, info, maxitr, mode1, nconv, ierr    !VARIABLES FOR DIAGONALIZATION ROUTINE
    logical, allocatable:: selector(:)
    double complex, allocatable :: ax(:), d(:), v(:,:), workd(:), workev(:), resid(:), workl(:)!VARIABLES FOR DIAGONALIZATION ROUTINE
    double precision, allocatable :: rwork(:), rd(:,:)

    !c
    !c     %-----------------------------%
    !c     | BLAS & LAPACK routines used |
    !c     %-----------------------------%
    !c

    double precision :: dznrm2 , dlapy2
    external :: dznrm2 , zaxpy , dlapy2, znaupd, zneupd
    double precision, external :: dnrm2



!c
!c     %-----------------------%
!c     | Executable Statements |
!c     %-----------------------%
!c
!c     %--------------------------------------------------%
!c     | The number NX is the number of interior points   |
!c     | in the discretization of the 2-dimensional       |
!c     | convection-diffusion operator on the unit        |
!c     | square with zero Dirichlet boundary condition.   |
!c     | The number N(=NX*NX) is the dimension of the     |
!c     | matrix.  A standard eigenvalue problem is        |
!c     | solved (BMAT = 'I').  NEV is the number of       |
!c     | eigenvalues to be approximated.  The user can    |
!c     | modify NX, NEV, NCV, WHICH to solve problems of  |
!c     | different sizes, and to get different parts of   |
!c     | the spectrum.  However, The following            |
!c     | conditions must be satisfied:                    |
!c     |                   N <= MAXN                      |
!c     |                 NEV <= MAXNEV                    |
!c     |           NEV + 2 <= NCV <= MAXNCV               |
!c     %--------------------------------------------------%
!c


    if (dim == 1) then
        evecs(1,1) = 1
        evals(1) = ham(1)
        print*, ''
        write(*,"('Only eigenvalue is ',f8.4)") ham(1)
    else

!    maxn = 10 + dim * dim !debug
    maxn = 10 + dim !* dim !debug
    maxnev = 10 + nev
    maxncv = 10 + ncv
    ldv = maxn

    if(allocated(workd)) deallocate(workd)
    if(allocated(workev)) deallocate(workev)
    if(allocated(rwork)) deallocate(rwork)
    if(allocated(rd)) deallocate(rd)
    if(allocated(ax)) deallocate(ax)
    if(allocated(d)) deallocate(d)
    if(allocated(resid)) deallocate(resid)
    if(allocated(selector)) deallocate(selector)
    if(allocated(workl)) deallocate(workl)
    if(allocated(v)) deallocate(v)
    allocate(ax(maxn), d(maxncv), resid(maxn), selector(maxncv), v(ldv,maxncv), &
            workl(3*maxncv*maxncv+5*maxncv), workd(3*maxn), workev(3*maxncv), &
            rwork(maxncv), rd(maxncv,3))


      nx = dim !debug
!      n  = nx * nx !debug
      n  = nx !* nx !debug


      if ( n .gt. maxn ) then
         print *, ' ERROR with _NDRV1: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NDRV1: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NDRV1: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat  = 'I'
      which = mode
!c
!c     %---------------------------------------------------%
!c     | The work array WORKL is used in ZNAUPD  as         |
!c     | workspace.  Its dimension LWORKL is set as        |
!c     | illustrated below.  The parameter TOL determines  |
!c     | the stopping criterion. If TOL<=0, machine        |
!c     | precision is used.  The variable IDO is used for  |
!c     | reverse communication, and is initially set to 0. |
!c     | Setting INFO=0 indicates that a random vector is  |
!c     | generated to start the ARNOLDI iteration.         |
!c     %---------------------------------------------------%
!c
    lworkl = 3*ncv**2+5*ncv
    tol    = 0.0
    ido    = 0
    info   = 0
    nconv  = 0


!c
!c     %---------------------------------------------------%
!c     | This program uses exact shift with respect to     |
!c     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!c     | IPARAM(3) specifies the maximum number of Arnoldi |
!c     | iterations allowed.  Mode 1 of ZNAUPD  is used     |
!c     | (IPARAM(7) = 1). All these options can be changed |
!c     | by the user. For details see the documentation in |
!c     | ZNAUPD .                                           |
!c     %---------------------------------------------------%
!c
    ishfts = 1
    maxitr = 300
    mode1  = 1
!c
    iparam(1) = ishfts
    iparam(3) = maxitr
    iparam(7) = mode1


!c
!c Create sparse matrix in CSR format
!c
    allocate(jao(nnz))
    allocate(iao(dim+1))
    allocate(ao(nnz))

    !call ccoocsr(dim,nnz,ham,rc(1:nnz,1),rc(1:nnz,2),ao,jao,iao)
    call ccoocsr(dim, nnz, ham, rc(1:nnz,1), rc(1:nnz,2), ao, jao, iao)


!c
!c     %-------------------------------------------%
!c     | M A I N   L O O P (Reverse communication) |
!c     %-------------------------------------------%
!c
    10   continue
!c
!c        %---------------------------------------------%
!c        | Repeatedly call the routine ZNAUPD  and take |
!c        | actions indicated by parameter IDO until    |
!c        | either convergence is indicated or maxitr   |
!c        | has been exceeded.                          |
!c        %---------------------------------------------%
!c
    call znaupd ( ido, bmat, n, which, nev, tol, resid, ncv, &
                v, ldv, iparam, ipntr, workd, workl, lworkl, &
                rwork, info )
!c
    if (ido .eq. -1 .or. ido .eq. 1) then
!c
!c           %-------------------------------------------%
!c           | Perform matrix vector multiplication      |
!c           |                y <--- OP*x                |
!c           | The user should supply his/her own        |
!c           | matrix vector multiplication routine here |
!c           | that takes workd(ipntr(1)) as the input   |
!c           | vector, and return the matrix vector      |
!c           | product to workd(ipntr(2)).               |
!c           %-------------------------------------------%
!c

        call cpamux (threads, nx, workd(ipntr(1)), workd(ipntr(2)), ao, jao, iao)

!c
!c           %-----------------------------------------%
!c           | L O O P   B A C K to call ZNAUPD  again. |
!c           %-----------------------------------------%
!c
        go to 10

    end if


!c
!c     %----------------------------------------%
!c     | Either we have convergence or there is |
!c     | an error.                              |
!c     %----------------------------------------%
!c
    if ( info .lt. 0 ) then
!c
!c        %--------------------------%
!c        | Error message, check the |
!c        | documentation in ZNAUPD   |
!c        %--------------------------%
!c
        print *, ' '
        print *, ' Error with znaupd, info = ', info
        print *, ' Check the documentation of _naupd'
        print *, ' '
!c
    else
!c
!c        %-------------------------------------------%
!c        | No fatal errors occurred.                 |
!c        | Post-Process using ZNEUPD .                |
!c        |                                           |
!c        | Computed eigenvalues may be extracted.    |
!c        |                                           |
!c        | Eigenvectors may also be computed now if  |
!c        | desired.  (indicated by rvec = .true.)    |
!c        %-------------------------------------------%
!c

!c
        call zneupd (rvec, 'A', selector, d, v, ldv, sigma, &
                    workev, bmat, n, which, nev, tol, resid, ncv, &
                    v, ldv, iparam, ipntr, workd, workl, lworkl, &
                    rwork, ierr)
!c
!c        %----------------------------------------------%
!c        | Eigenvalues are returned in the one          |
!c        | dimensional array D.  The corresponding      |
!c        | eigenvectors are returned in the first NCONV |
!c        | (=IPARAM(5)) columns of the two dimensional  |
!c        | array V if requested.  Otherwise, an         |
!c        | orthogonal basis for the invariant subspace  |
!c        | corresponding to the eigenvalues in D is     |
!c        | returned in V.                               |
!c        %----------------------------------------------%
!c


        if ( ierr .ne. 0) then
!c
!c           %------------------------------------%
!c           | Error condition:                   |
!c           | Check the documentation of ZNEUPD . |
!c           %------------------------------------%
!c
            print *, ' '
            print *, ' Error with zneupd, info = ', ierr
            print *, ' Check the documentation of _neupd. '
            print *, ' '
!c
        else
!c
            nconv = iparam(5)
            do 20 j = 1, nconv
!c
!c               %---------------------------%
!c               | Compute the residual norm |
!c               |                           |
!c               |   ||  A*x - lambda*x ||   |
!c               |                           |
!c               | for the NCONV accurately  |
!c               | computed eigenvalues and  |
!c               | eigenvectors.  (iparam(5) |
!c               | indicates how many are    |
!c               | accurate to the requested |
!c               | tolerance)                |
!c               %---------------------------%
!c

                call cpamux (threads, nx, v(1,j), ax, ao, jao, iao)
                !call av (threads, nup, ndn, ndi, hamup, hamdn, hamdi, dim, v(1,j), ax)
                !call av(nx, v(1,j), ax)
                call zaxpy (n, -d(j), v(1,j), 1, ax, 1)
                rd(j,1) = dble (d(j))
                rd(j,2) = dimag (d(j))
                rd(j,3) = dznrm2 (n, ax, 1)
                rd(j,3) = rd(j,3) / dlapy2 (rd(j,1),rd(j,2))
            20 continue
!c
!c            %-----------------------------%
!c            | Display computed residuals. |
!c            %-----------------------------%
!c
            call dmout (6, nconv, 3, rd, maxncv, -6, &
                        'Ritz values (Real, Imag) and relative residuals')
        end if


!c
!c        %-------------------------------------------%
!c        | Print additional convergence information. |
!c        %-------------------------------------------%
!c
        if ( info .eq. 1) then
            print *, ' '
            print *, ' Maximum number of iterations reached.'
            print *, ' '
        else if ( info .eq. 3) then
            print *, ' '
            print *, ' No shifts could be applied during implicit', &
                     ' Arnoldi update, try increasing NCV.'
            print *, ' '
        end if
!c
        print *, ' '
        print *, '_NDRV1'
        print *, '====== '
        print *, ' '
        print *, ' Size of the matrix is ', n
        print *, ' The number of Ritz values requested is ', nev
        print *, ' The number of Arnoldi vectors generated', &
                 ' (NCV) is ', ncv
        print *, ' What portion of the spectrum: ', which
        print *, ' The number of converged Ritz values is ', &
                   nconv
        print *, ' The number of Implicit Arnoldi update', &
                 ' iterations taken is ', iparam(3)
        print *, ' The number of OP*x is ', iparam(9)
        print *, ' The convergence criterion is ', tol
        print *, ' '
!c
    end if




    evals(1:nev) = d(1:nev)


    if (rvec) then
        do j= 1, nst
            evecs(1:dim, j) = v(1:dim,j)
        end do
    end if

    end if


!c
!c     %---------------------------%
!c     | Done with program zndrv1 . |
!c     %---------------------------%
!c

!c          Error flag for ZNAUPD on output.
!c          =  0: Normal exit.
!c          =  1: Maximum number of iterations taken.
!c                All possible eigenvalues of OP has been found. IPARAM(5)
!c                returns the number of wanted converged Ritz values.
!c          =  2: No longer an informational error. Deprecated starting
!c                with release 2 of ARPACK.
!c          =  3: No shifts could be applied during a cycle of the
!c                Implicitly restarted Arnoldi iteration. One possibility
!c                is to increase the size of NCV relative to NEV.
!c                See remark 4 below.
!c          = -1: N must be positive.
!c          = -2: NEV must be positive.
!c          = -3: NCV-NEV >= 2 and less than or equal to N.
!c          = -4: The maximum number of Arnoldi update iteration
!c                must be greater than zero.
!c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
!c          = -6: BMAT must be one of 'I' or 'G'.
!c          = -7: Length of private work array is not sufficient.
!c          = -8: Error return from LAPACK eigenvalue calculation;
!c          = -9: Starting vector is zero.
!c          = -10: IPARAM(7) must be 1,2,3.
!c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
!c          = -12: IPARAM(1) must be equal to 0 or 1.
!c          = -9999: Could not build an Arnoldi factorization.
!c                   User input error highly likely.  Please
!c                   check actual array dimensions and layout.
!c                   IPARAM(5) returns the size of the current Arnoldi
!c                   factorization.


!    Error flag for ZNEUPD on output.
!    0: Normal exit.
!    1: The Schur form computed by LAPACK routine csheqr could not be reordered by LAPACK routine ztrsen.
!       Re-enter subroutine zneupd with IPARAM(5) = NCV and increase the size of the array D to have dimension at least dimension NCV and allocate at least NCV columns for Z.
!       NOTE: Not necessary if Z and V share the same space. Please notify the authors if this error occurs.
!    -1: N must be positive.
!    -2: NEV must be positive.
!    -3: NCV-NEV >= 1 and less than or equal to N.
!    -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'.
!    -6: BMAT must be one of 'I' or 'G'.
!    -7: Length of private work WORKL array is not sufficient.
!    -8: Error return from LAPACK eigenvalue calculation. This should never happened.
!    -9: Error return from calculation of eigenvectors. Informational error from LAPACK routine ztrevc.
!    -10: IPARAM(7) must be 1, 2, 3.
!    -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
!    -12: HOWMANY = 'S' not yet implemented.
!    -13: HOWMANY must be one of 'A' or 'P' if RVEC = .true.
!    -14: ZNAUPD did not find any eigenvalues to sufficient accuracy.
!    -15: ZNEUPD got a different count of the number of converged Ritz values than ZNAUPD got. This indicates the user probably made an error in passing data from ZNAUPD to ZNEUPD or that the data was modified before entering ZNEUPD.

9000 continue


    if(allocated(workd)) deallocate(workd)
    if(allocated(workev)) deallocate(workev)
    if(allocated(rwork)) deallocate(rwork)
    if(allocated(rd)) deallocate(rd)
    if(allocated(ax)) deallocate(ax)
    if(allocated(d)) deallocate(d)
    if(allocated(resid)) deallocate(resid)
    if(allocated(selector)) deallocate(selector)
    if(allocated(workl)) deallocate(workl)
    if(allocated(v)) deallocate(v)
    if(allocated(jao)) deallocate(jao)
    if(allocated(iao)) deallocate(iao)
    if(allocated(ao)) deallocate(ao)



end subroutine




!--------------------------------------------------------!
!            COO to CSR for complex sparse matrix        !
!--------------------------------------------------------!


subroutine ccoocsr (nrow,nnz,a,ir,jc,ao,jao,iao)
    implicit none

!    integer, intent(in) :: nrow, nnz
    integer(kind=8), intent(in) :: nrow
    integer, intent(in) :: nnz
    double complex :: a(*),ao(*),x
    integer :: ir(*),jc(*),jao(*),iao(*)

    integer :: k, k0, i, j, iad


    do 1 k = 1, nrow + 1
        iao(k) = 0
    1 continue

    do 2 k = 1, nnz
        iao(ir(k)) = iao(ir(k)) + 1
    2 continue

    k = 1
    do 3 j = 1, nrow + 1
        k0 = iao(j)
        iao(j) = k
        k = k + k0
    3 continue

    do 4 k = 1, nnz
        i = ir(k)
        j = jc(k)
        x = a(k)
        iad = iao(i)
        ao(iad) =  x
        jao(iad) = j
        iao(i) = iad + 1
    4 continue

    do 5 j = nrow, 1, -1
        iao(j+1) = iao(j)
    5 continue
    iao(1) = 1
    return

end subroutine ccoocsr




!---------------------------------!
!            Thread units         !
!---------------------------------!

subroutine threadunits(ndv, ndv2, units)

    implicit none
    integer, intent(in) :: ndv, ndv2
    integer, allocatable, intent(out) :: units(:, :)
    integer :: i, j, h
    integer :: n_threads_tot

    n_threads_tot = 1


    if (ndv > 0) n_threads_tot = n_threads_tot + ndv
    if (ndv2 > 0) n_threads_tot = n_threads_tot + ndv2




    if(allocated(units)) deallocate(units)
    allocate(units(n_threads_tot,n_threads_tot))
    h = 0
    do i = 1, n_threads_tot
        do j = 1, n_threads_tot
            h = h + 1
            units(i,j) = h
        end do
    end do


end subroutine



!-------------------------------!
!            Step units         !
!-------------------------------!

subroutine stepunits(ndu, ndv, ndv2, units)
    !Creates array with units for parallel I/O. Each entry is a different unit number.
    implicit none
    integer, intent(in) :: ndu, ndv, ndv2
    integer, allocatable, intent(out) ::  units(:,:,:,:)
    integer :: i, j, k, l, q

    if (allocated(units)) deallocate(units)
    allocate(units(ndv2,ndv,ndu,2))
    q = 0
    do i = 1, ndv2
        do j = 1, ndv
            do k = 1, ndu
                do l = 1, 2
                    q = q + 1
                    units(i,j,k,l) = q
                end do
            end do
        end do
    end do

end subroutine



!-------------------------------------!
!            Print parameters         !
!-------------------------------------!

subroutine printparams(dim, sites, pts, ti, mom, nthreads, dim_hs, nev, ndis, v1list, v2list, &
    t, dis, dv, vmax, vmin, dv2, v2min, v2max)

    implicit none
    integer, intent(in)          :: dim, sites, pts, ti, mom, nthreads, nev, ndis, v1list, v2list
    integer(kind=8), intent(in)  :: dim_hs
    double precision, intent(in) :: t, dis, dv, vmax, vmin, dv2, v2min, v2max

    print*, '!-------------------------------------------------------------!'
    print*, '!                                                             !'
    print*, '!                          Parameters                         !'
    print*, '!                                                             !'
    print*, '!-------------------------------------------------------------!'
    print*, ''
    print('(1x, a,i0)'), 'Dimension = ', dim
    print*, ''
    print('(1x, a,i0)'), 'Sites = ', sites
    print*, ''
    print('(1x, a,i0)'), 'Particles = ', pts
    print*, ''
    if (ti == 1) print('(1x, a, i0)'), 'Momentum k = ', mom
    if (ti == 1) print*, ''
    print('(1x, a,i0)'), 'Number of Threads = ', nthreads
    print*, ''
    print('(1x, a,f10.4)'), 'Hopping amplitude t = ', t
    print*, ''
    print('(1x, a,i0)'), 'Hilbert space dimension = ', dim_hs
    print*, ''
    print('(1x, a,i0)'), 'Number of eigenvalues = ', nev
    print*, ''
    print('(1x, a,f10.4)'), 'Disorder strength W = ', dis
    print*, ''
    print('(1x, a,i0)'), 'Number of disorder realizations = ', ndis
    print*, ''
    if(v1list == 1) then
        print('(1x, a,f10.4)'), 'Stepsize dv = ',dv
        print*, ''
        print('(1x, a,f10.4)'), 'Vmax = ', vmax
        print*, ''
        print('(1x, a,f10.4)'), 'Vmin = ', vmin
    end if
    if(v2list == 1) then
        print*, ''
        print('(1x, a,f10.4)'), 'Stepsize dv2 = ',dv2
        print*, ''
        print('(1x, a,f10.4)'), 'V2max = ', v2max
        print*, ''
        print('(1x, a,f10.4)'), 'V2min = ', v2min
        print*, ''
    end if
    print*, '!-------------------------------------------------------------!'
    print*, '!-------------------------------------------------------------!'
    print*, ''

return
end subroutine



!-------------------------------------------!
!            Date and Time                  !
!-------------------------------------------!

subroutine datetime(startend, threads, params)
    implicit none

    integer, intent(in) :: startend
    integer, intent(in), optional :: threads
    character*(*), intent(in), optional :: params
    character(8), save  :: datei, datef
    character(10), save :: timei, timef
    integer,dimension(8), save :: valuesi, valuesf
    real, save :: start, finish
    character :: tempchar*100, file_name*150



    if(startend == 0) then
        call cpu_time(start)
        call date_and_time(date=datei,time=timei,values=valuesi)
        write(*,"('Calculation started at',x,a,' h',x,a,' min',x,a,' sec')") timei(1:2), timei(3:4), timei(5:6)
        print*, ''
        print*, 'Start date: ',datei(7:8), '.',datei(5:6), '.',datei(1:4)
        print*, ''
    else
        call cpu_time(finish)
        print*, 'Start = ', start
        print*, 'Finish = ', finish
        if (finish - start < 60) then
            write(*,"(' Elapsed CPU time = ',f12.3,' seconds.')") finish-start
            print*, ''
        else if (finish - start < 3600) then
            write(*,"(' Elapsed CPU time = ',f12.3,' minutes.')") (finish-start)/60
            print*, ''
        else
            write(*,"(' Elapsed CPU time = ',f12.3,' hours.')") (finish-start)/3600
            print*, ''
        end if
        call date_and_time(date=datef,time=timef,values=valuesf)
        write(*,"(' Calculation started at',x,a,'h',x,a,'min',x,a,'sec')") timei(1:2), timei(3:4), timei(5:6)
        print*, ''
        write(*,"(' Calculation ended at',x,a,'h',x,a,'min',x,a,'sec')") timef(1:2), timef(3:4), timef(5:6)
        print*, ''
        print*, 'Start date: ',datei(7:8), '.',datei(5:6), '.',datei(1:4)
        print*, ''
        print*, 'End date: ',datef(7:8), '.',datef(5:6), '.',datef(1:4)
        file_name = "times_"//params
        file_name = trim_name(file_name)
        open(77,file = file_name)
        write(77,"(' Number of threads = ',i0)") threads
        write(77,"(' Elapsed CPU time = ',f20.10,' seconds.')") finish-start
        write(77,"(' Calculation started at',x,a,'h',x,a,'min',x,a,'sec')") timei(1:2), timei(3:4), timei(5:6)
        write(77,"(' Start date:',x,a,'.',x,a,'.',x,a)") ,datei(7:8), datei(5:6), datei(1:4)
        write(77,"(' Calculation ended at',x,a,'h',x,a,'min',x,a,'sec')") timef(1:2), timef(3:4), timef(5:6)
        write(77,"(' End date:',x,a,'.',x,a,'.',x,a)") ,datef(7:8), datef(5:6), datef(1:4)
        close(77)
    end if
end subroutine datetime





!---------------------------------------------------!
!            Number of discretization steps         !
!---------------------------------------------------!


subroutine nsteps(min, max, delta, steps)

    implicit none
    double precision, intent(in) :: min, max, delta
    integer, intent(out) :: steps

    if(delta <= 1.d0) then
        steps = int(abs(max-min)/delta + delta/2) + 1
    else
        steps = int(abs(max-min)/delta) + 1
    end if


end subroutine nsteps



!------------------------------------------!
!            Check parallelization         !
!------------------------------------------!

subroutine parallelcheck()

    implicit none
    integer :: num_threads_loc
    integer :: thread_num_loc
    integer :: max_threads_loc

    num_threads_loc = 0
    thread_num_loc = 0
    max_threads_loc = 0
    !Critical block is a lock, forcing the instructions to be run in serial
    !$omp parallel
        !$omp critical
            !$ thread_num_loc = omp_get_thread_num()
            print('(1x,100(a,i0))'), 'Thread number ', thread_num_loc, ' is online.'
            print*, ''
            !$ num_threads_loc = omp_get_num_threads()
            if (thread_num_loc == 0) then
                print('(1x,100(a,i0))'), 'Available number of threads = ', num_threads_loc
                print*, ''
                !$ max_threads_loc = omp_get_max_threads()
                print('(1x,100(a,i0))'), 'Maximum number of threads = ', max_threads_loc
                print*, ''
            end if
        !$omp end critical
    !$omp end parallel

end subroutine




!--------------------------------------------------------!
!            Find connected components (Dense)           !
!--------------------------------------------------------!

subroutine digraph_adj_components ( adj, lda, nnode, ncomp, comp, dad, order )

!*****************************************************************************
!
!! DIGRAPH_ADJ_COMPONENTS finds the strongly connected components of a digraph.
!
!  Discussion:
!
!    A digraph is a directed graph.
!
!    A strongly connected component of a directed graph is the largest
!    set of nodes such that there is a directed path from any node to
!    any other node in the same component.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Reference:
!
!    K Thulasiraman, M Swamy,
!    Graph Theory and Algorithms,
!    John Wiley, New York, 1992.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ADJ.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) NCOMP, the number of strongly connected
!    components.
!
!    Output, integer ( kind = 4 ) COMP(NNODE), lists the connected component
!    to which each node belongs.
!
!    Output, integer ( kind = 4 ) DAD(NNODE), the father array for the depth
!    first search trees.  DAD(I) = 0 means that node I is the root of
!    one of the trees.  DAD(I) = J means that the search descended
!    from node J to node I.
!
!    Output, integer ( kind = 4 ) ORDER(NNODE), the order in which the nodes
!    were traversed, from 1 to NNODE.
!
  implicit none

!  integer ( kind = 4 ) lda
!  integer ( kind = 4 ) nnode
  integer ( kind = 8 ) lda
  integer ( kind = 8 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) comp(nnode)
  integer ( kind = 4 ) dad(nnode)
  integer ( kind = 4 ) iorder
  integer ( kind = 4 ) lowlink(nnode)
  integer ( kind = 4 ) mark(nnode)
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) point(nnode)
  integer ( kind = 4 ) stack(nnode)
  integer ( kind = 4 ) v
  integer ( kind = 4 ) w
  integer ( kind = 4 ) x
!
!  Initialization.
!
  comp(1:nnode) = 0
  dad(1:nnode) = 0
  order(1:nnode) = 0
  lowlink(1:nnode) = 0
  mark(1:nnode) = 0
  point(1:nnode) = 0

  iorder = 0
  nstack = 0
  ncomp = 0
!
!  Select any node V not stored in the stack, that is, with MARK(V) = 0.
!
  do

    v = 0

    do

      v = v + 1

      if ( nnode < v ) then
        adj(1:nnode,1:nnode) = abs( adj(1:nnode,1:nnode) )
        return
      end if

      if ( mark(v) /= 1 ) then
        exit
      end if

    end do

    iorder = iorder + 1

    order(v) = iorder
    lowlink(v) = iorder
    mark(v) = 1

    nstack = nstack + 1
    stack(nstack) = v
    point(v) = 1

30  continue
!
!  Consider each node W.
!
    do w = 1, nnode
!
!  Is there an edge (V,W) and has it not been examined yet?
!
      if ( 0 < adj(v,w) ) then

        adj(v,w) = - adj(v,w)
!
!  Is the node on the other end of the edge undiscovered yet?
!
        if ( mark(w) == 0 ) then

          iorder = iorder + 1
          order(w) = iorder
          lowlink(w) = iorder
          dad(w) = v
          mark(w) = 1

          nstack = nstack + 1
          stack(nstack) = w
          point(w) = 1

          v = w

        else if ( mark(w) == 1 ) then

          if ( order(w) < order(v) .and. point(w) == 1 ) then
            lowlink(v) = min ( lowlink(v), order(w) )
          end if

        end if

        go to 30

      end if

    end do

    if ( lowlink(v) == order(v) ) then

      ncomp = ncomp + 1

      do

        if ( nstack <= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DIGRAPH_ADJ_COMPONENTS - Fatal error!'
          write ( *, '(a)' ) '  Illegal stack reference.'
          stop
        end if

        x = stack(nstack)
        nstack = nstack - 1

        point(x) = 0
        comp(x) = ncomp

        if ( x == v ) then
          exit
        end if

      end do

    end if

    if ( dad(v) /= 0 ) then
      lowlink(dad(v)) = min ( lowlink(dad(v)), lowlink(v) )
      v = dad(v)
      go to 30
    end if

  end do

  return
end




end module module_fragmentation
