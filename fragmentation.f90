program fragmentation

    use mkl_pardiso
    use module_fragmentation
    use module_variables_fragmentation

    implicit none

    character, parameter :: mode*2 = 'SR' !'LM', 'LR', 'SR', 'LI', 'SI'
    double precision, parameter :: eps = 1.0d0

    integer(kind=8), allocatable :: mombasis(:), periods(:), basiss(:), nonreps(:), allstates(:,:)
    integer(kind=8), allocatable :: permutations(:)
    integer(kind=8), allocatable :: bjlr(:,:)
    integer, allocatable :: occupation(:,:)
    integer, allocatable :: rc_off(:,:), pseudo_rc(:,:)
    integer, allocatable :: rc_di(:)
    integer, allocatable :: rc(:,:)
    integer, allocatable :: adj(:,:), fakeadj(:,:)
    integer, allocatable :: comp(:), compv1(:), compv2(:)
    integer, allocatable :: compsize(:,:)
    integer, allocatable :: dad(:)
    integer, allocatable :: order(:)

    double complex, allocatable, save :: cenergy(:,:,:,:)
    double complex, allocatable :: work(:)
    double complex, allocatable :: eigstate(:,:,:)
    double complex, allocatable :: ham_off(:), pseudo_ham(:)
    double complex, allocatable :: ham(:), ham_d(:,:), evecs(:,:)

    double precision, allocatable :: rwork(:), ham_di(:,:)
    double precision, allocatable :: evals(:)
    double precision, allocatable :: weightsv1(:), weightsv2(:)
    double precision, allocatable :: singval(:), entropy(:), aventropy(:), totaventropy(:,:,:,:,:), totavtemp(:)
    double precision, allocatable :: iprv1(:), iprv2(:)
    double precision, allocatable :: v1_list(:)
    double precision, allocatable :: v2_list(:)

    integer(kind=8) :: la = 0
    integer(kind=8) :: dim_hs
    integer(kind=8) :: momdim = 0
    integer(kind=8) :: dimtot = 0
    integer(kind=8) :: rep = 0
    integer :: n = 0, nev = 0, ncv = 0, sites = 0
    integer :: pts = 0
    integer :: ndv = 0, ndv2 = 0
    integer :: h = 0, i = 0, j = 0, k = 0, l = 0, m = 0, q = 0, r = 0, s = 0, ss = 0, ne = 0, nu = 0, nb = 0, nv = 0, nv2 = 0 !LOOP PARAMETERS
    integer :: counter = 0, row = 0, col = 0, ls1 = 0, ls2 = 0
    integer :: full = 0
    integer :: temp = 0
    integer :: ti = 1
    integer :: n_off = 0
    integer :: n_di = 0
    integer :: nnz = 0
    integer :: nest = 0
    integer :: lwork = 0, info = 0
    integer :: ww = 0
    integer :: mom = 0, nmom = 0
    integer :: dimhs = 0, lwmax = 0
    integer :: npts = 0
    integer :: v1list = 0
    integer :: v2list = 0
    integer :: eecut = 0
    integer :: cntr = 0
    integer :: pseudoi = 0
    integer :: ncomp = 0, ncompv1 = 0, ncompv2 = 0
    integer :: vi  = 0 !SPARTAN: COMMENT OUT, EXTERNALLY PROVIDED
    integer :: v2i = 0 !SPARTAN: COMMENT OUT, EXTERNALLY PROVIDED

    double precision :: delta    = 0.d0
    double precision :: rmean    = 0.d0
    double precision :: thresh   = 0.d0
    double precision :: aviprv1  = 0.d0, aviprv2  = 0.d0
    double precision :: vv = 0, v2 = 0, mu = 0, e0 = 0

    logical :: check
    logical :: ok

    character :: globalpar*200, param*200, parameters*200, &
                parameters_short*200, parameters_medium*200, file_name*200

    !Variables for date and time
    character(8)         :: datei, datef
    character(10)        :: timei, timef
    integer,dimension(8) :: valuesi, valuesf

    !OpenMP VARIABLES
    integer, allocatable :: units(:,:), units2(:,:,:,:)
    integer :: num_threads, max_threads,  hamthreads
    integer :: ompt = 0
    integer :: level
    integer :: team_size
    integer :: thread_num = 0, thread_num2 = 0, thread_num3 = 0
    integer :: stack_size_check
    logical :: in_parallel
    logical :: get_nested
    logical :: get_dynamic

    !$OMP THREADPRIVATE (nv, nv2)
    nv  = 0
    nv2 = 0

    v1list = 1
    v2list = 1
    allocate(v1_list(23))
    do i = 1, 22
        v1_list(i+1) = 0.1 * 1.5**(real(i)-1)
    end do
    allocate(v2_list(23))
    do i = 1, 22
        v2_list(i+1) = 0.1 * 1.5**(real(i)-1)
    end do
    allocate(v1_list(1))
    v1_list(1) = 0

    call datetime(0)

    write(*,*) '--------------------------------------------------------------------------------------------------'
    write(*,*) ' Exact diagonalization '
    write(*,*) '--------------------------------------------------------------------------------------------------'

    call slsetvars(1, dynamic, nested, nthreads, sites, pts, num_threads, filling)
    if ( modulo( int( sites )/2, 2 ) == 0 ) pts = int( sites/2 ) - 1

    !$ call omp_set_num_threads(nthreads)
    !$omp parallel
    !$ num_threads= omp_get_num_threads()
    !$omp end parallel

    call parallelcheck()
    if (v1list == 0) then
        call nsteps(vmin, vmax, dv, ndv)
    else
        ndv = size(v1_list)
    end if
    if (v2list == 0) then
        call nsteps(v2min, v2max, dv2, ndv2)
    else
        ndv2 = size(v2_list)
    end if

    call stepunits(1, ndv, ndv2, units2)

    if (ee == 1) eecut = int(sites/2) !sites - 1

    !if (ee == 1 .and. ti == 1 .and. bc == 'p') then
    !    if (allocated(totaventropy)) deallocate(totaventropy)
    !    allocate(totaventropy(eecut, sites, ndv, ndv2, ndis))
    !    totaventropy = 0
    !end if

    !do mom = -int(sites/2) + 1, int(sites/2) !Momentum loop
    nmom = nmom + 1

    call threadunits(ndv, ndv2, units)

    mom = momext
    print*, 'Momentum sector = 2pi/L*', mom
    print*, ''

    if (bc == 'o') ti = 0
    call slbasis(25, sites, pts, dim_hs, permutations)

    if (ip == 1) then
        do k = 1, 2
            if ( not( allocated( adj ) ) )   allocate( adj( dim_hs,dim_hs ) )
            if ( not( allocated( comp ) ) )  allocate( comp( dim_hs ) )
            if ( not( allocated( dad ) ) )   allocate( dad( dim_hs ) )
            if ( not( allocated( order ) ) ) allocate( order( dim_hs ) )
            if ( not( allocated( periods ) ) ) allocate( periods( dim_hs ) )
            adj     = 0
            comp    = 0
            dad     = 0
            order   = 0
            periods = 0
            !Generate first adjacency matrix
            call slhopping_serial ( dim_hs, sites, bc, 0, k, eps, t, permutations, pseudo_rc, pseudo_ham, pseudoi, mom, periods, adj )
            call digraph_adj_components ( adj, dim_hs, dim_hs, ncomp, comp, dad, order )
            deallocate( periods )
            !Generate list 'compsize' assigning size of CC, bond number of CC and a representative state
            if( allocated( compsize ) ) deallocate( compsize )
            allocate( compsize( ncomp, 3 ) ) !1st col: Size of CC; 2nd col: Bond number of CC; 3rd col: Representative of CC; Row # = CC #
            compsize = 0
            do i = 1, dim_hs
                if ( compsize( comp(i), 1 ) == 0 ) then
                    compsize( comp(i), 3 ) = i
                    do j = 0, sites - 1
                        if ( k == 1 .and. btest( permutations(i), j ) .and. btest( permutations(i), modulo(j + 1, sites))) compsize(comp(i), 2) = compsize(comp(i), 2) + 1
                        if ( k == 2 .and. btest( permutations(i), j ) .and. btest( permutations(i), modulo(j + 2, sites))) compsize(comp(i), 2) = compsize(comp(i), 2) + 1
                    end do
                end if
                compsize( comp(i), 1 ) = compsize( comp(i), 1 ) + 1
            end do
            do i = 1, ncomp - 1
                do j = i + 1, ncomp
                    if ( compsize( i, 1 ) == compsize( j, 1 ) .and. compsize( i, 2 ) == compsize( j, 2 ) ) then
                        adj( compsize( i, 3 ), compsize( j, 3 ) ) = 1
                        adj( compsize( j, 3 ), compsize( i, 3 ) ) = 1
                    end if
                end do
            end do
            ncomp = 0
            comp  = 0
            dad   = 0
            order = 0
            !Calculate second adjacency matrix
            call digraph_adj_components ( adj, dim_hs, dim_hs, ncomp, comp, dad, order )
            deallocate( dad, order, pseudo_ham, pseudo_rc  )
            if ( k == 1 ) then
                ncompv1 = ncomp
                compv1 = comp
            else if ( k == 2 ) then
                ncompv2 = ncomp
                compv2 = comp
            end if
        end do
    end if

    if (ti == 1 .and. bc == 'p') then
        call momentumbasis(dim_hs, sites, mom, permutations, momdim, mombasis, periods, nonreps)
        dim_hs = momdim
        call move_alloc(mombasis,basiss)
        if (ee == 1 .or. ip == 1) then
            dimtot = momdim + size(nonreps)
            allocate(allstates(dimtot, 5)) !allstates contains all momentum states: 1 = state, 2 = position of representative in momentum basis, 3 = # translations, 4 = period
            allstates = 0
            do i = 1, momdim
                allstates(i, 1) = basiss(i)
                allstates(i, 2) = i
                allstates(i, 3) = 0
                allstates(i, 4) = periods(i)
                call findstate(dimtot, basiss(i), permutations, allstates(i, 5))
            end do
            do i = momdim + 1, dimtot
                allstates(i,1) = nonreps(i - momdim)
                call representative(nonreps(i - momdim), sites, rep, allstates(i,3))
                call checkstate(rep, sites, mom, allstates(i,4))
                call findstate(momdim, rep, basiss, allstates(i,2))
                call findstate(dimtot, nonreps(i - momdim), permutations, allstates(i, 5))
            end do
        end if
        deallocate(nonreps)
    else
        call move_alloc(permutations,basiss)
        if (allocated(periods)) deallocate(periods)
        allocate(periods(dim_hs))
        periods = 1
    end if

    print*, 'Creating Hamiltonian matrix ...'
    print*, ''
    !!$omp parallel private(thread_num) num_threads(min(2,num_threads))
    !    !$omp sections
    !        !$omp section
    !            call slhopping(dim_hs, sites, ti, bc, eps, t, basiss, rc_off, ham_off, n_off, mom, periods)
    !        !$omp section
    !            call slhamdi (sites, dim_hs, rc_di, ham_di, basiss, occupation, n_di)
    !    !$omp end sections
    !!$omp end parallel
    allocate(fakeadj(dim_hs, dim_hs))
    call slhopping_serial(dim_hs, sites, bc, ti, ip, eps, t, basiss, rc_off, ham_off, n_off, mom, periods, fakeadj)
    call slhamdi (sites, dim_hs, rc_di, ham_di, basiss, occupation, n_di)
    deallocate(fakeadj)

    if(allocated(permutations)) deallocate(permutations)
    if(allocated(mombasis)) deallocate(mombasis)
    if(allocated(periods)) deallocate(periods)

    call nevncv (dimthresh, exact, nevext, n_st, ncv0, dim_hs, full, nev, ncv, nest)
    if (ee == 1 .and. exact == 1) nest = dim_hs
    if (ip == 1 .and. exact == 1) nest = dim_hs

    if (allocated(cenergy)) deallocate(cenergy)
    allocate(cenergy(nev, ndv, ndv2, ndis))
    cenergy = 0

    if (dim_hs == 0) then
        nev = dim_hs
        ncv = dim_hs
        goto 116
    end if

    num_threads = max(int((nthreads - min(ndv2,nthreads) * max(min(ndv,int((nthreads-min(ndv2,nthreads))/min(ndv2,nthreads))),1)) &
                / (min(ndv2,nthreads) * max(min(ndv,int((nthreads-min(ndv2,nthreads))/min(ndv2,nthreads))),1))),1)

    print*, 'Number of V2 threads =',min(ndv2,nthreads)
    print*, ''
    print*, 'Number of V1 threads =',max(min(ndv,int((nthreads-min(ndv2,nthreads))/min(ndv2,nthreads))),1)
    print*, ''
    print*, 'Number of threads left = ', num_threads
    print*, ''
    if(num_threads < 1) error stop "NO THREADS LEFT! NUMBER OF THREADS"

    !--------------------------!
    !         V2LOOP           !
    !--------------------------!

    !$omp parallel do default(firstprivate) private(i, j ,la, ww, ham_d, eigstate) shared(allstates, units, ham_off, ham_di, occupation, basiss, rc_off, rc_di, cenergy) num_threads(min(ndv2,nthreads))
    do nv2 = 0, ndv2-1, 1
        if (v2list == 0) then
            !v2 = v2min + nv2*dv2
            v2  = 0! v2ext !SPARTAN
        else
            v2 = v2_list(nv2 + 1)
        end if

        print*, ''
        print('(1x,100(a,f12.4))'), '----------------------  V2 = ', v2, ' -------------------------'
        print*, ''

        !$ thread_num = omp_get_thread_num()

        !--------------------------!
        !         V1LOOP           !
        !--------------------------!

        !$omp parallel do default(firstprivate) private(i, j, la, ww, ham_d, eigstate) shared(allstates, units, ham_off, ham_di, occupation, basiss, rc_off, rc_di, cenergy) num_threads(max(min(ndv,int((nthreads-min(ndv2,nthreads))/min(ndv2,nthreads))),1))
        do nv = 0, ndv-1, 1

            if (v1list == 0) then
                !vv = vmin + nv*dv
                vv = 0!vext !SPARTAN
            else
                vv = v1_list(nv + 1)
            end if

            print*, ''
            print('(1x,100(a,f12.4))'), '----------------------  V = ', vv, ' -------------------------'
            print*, ''

            !$ thread_num2 = omp_get_thread_num()

            if (ti == 0 .and. bc == 'p') then
                if ( v1list == 1 .and. v2list == 0 ) then
                    write (parameters,"('L=',i0,'N=',i0,'t=',f12.4,'V=',i0, &
                        'V2=',i0,'W=',f12.4,'eps=',f12.4,'ndis=',i0,'BC=',a,'.dat')") sites,pts,t,nv, &
                        v2i,dis,eps,ndis,bc
                else if ( v1list == 0 .and. v2list == 1 ) then
                    write (parameters,"('L=',i0,'N=',i0,'t=',f12.4,'V=',i0, &
                        'V2=',i0,'W=',f12.4,'eps=',f12.4,'ndis=',i0,'BC=',a,'.dat')") sites,pts,t,vi, &
                        nv2,dis,eps,ndis,bc
                else if ( v1list == 1 .and. v2list == 1 ) then
                    write (parameters,"('L=',i0,'N=',i0,'t=',f12.4,'V=',i0, &
                        'V2=',i0,'W=',f12.4,'eps=',f12.4,'ndis=',i0,'BC=',a,'.dat')") sites,pts,t,nv, &
                        nv2,dis,eps,ndis,bc
                else
                    write (parameters,"('L=',i0,'N=',i0,'t=',f12.4,'V=',f12.4, &
                        'V2=',f12.4,'W=',f12.4,'eps=',f12.4,'ndis=',i0,'BC=',a,'.dat')") sites,pts,t,vv, &
                        v2,dis,eps,ndis,bc
                end if
            else
                if ( v1list == 1 .and. v2list == 0 ) then
                    write (parameters,"('L=',i0,'N=',i0,'k=',i0,'t=',f12.4,'V=',i0, &
                        'V2=',i0,'W=',f12.4,'eps=',f12.4,'ndis=',i0,'BC=',a,'.dat')") sites,pts,mom,t,nv, &
                        v2i,dis,eps,ndis,bc
                else if ( v1list == 0 .and. v2list == 1 ) then
                    write (parameters,"('L=',i0,'N=',i0,'k=',i0,'t=',f12.4,'V=',i0, &
                        'V2=',i0,'W=',f12.4,'eps=',f12.4,'ndis=',i0,'BC=',a,'.dat')") sites,pts,mom,t,vi, &
                        nv2,dis,eps,ndis,bc
                else if ( v1list == 1 .and. v2list == 1 ) then
                    write (parameters,"('L=',i0,'N=',i0,'k=',i0,'t=',f12.4,'V=',i0, &
                        'V2=',i0,'W=',f12.4,'eps=',f12.4,'ndis=',i0,'BC=',a,'.dat')") sites,pts,mom,t,nv, &
                        nv2,dis,eps,ndis,bc
                else
                    write (parameters,"('L=',i0,'N=',i0,'t=',f12.4,'V=',f12.4, &
                        'V2=',f12.4,'W=',f12.4,'eps=',f12.4,'ndis=',i0,'BC=',a,'.dat')") sites,pts,t,vv, &
                        v2,dis,eps,ndis,bc
                end if
            end if
            parameters=trim_name(parameters)
            print*, 'Parameters = ', parameters
            print*, ''

            !-----------------------------------------!
            !         Disorder realizations           !
            !-----------------------------------------!

            do ww = 1, ndis
                if (full == 0 ) then
                    call slchamunify (eps, vv, v2, dis, sites, occupation, n_off, n_di, ham_off, rc_off, ham_di, rc_di, ham, rc, nnz)
                else if (full == 1) then
                    call slchamdense (dim_hs, eps, vv, v2, dis, sites, occupation, n_off, n_di, ham_off, rc_off, ham_di, rc_di, ham_d, nnz)
                    delta  = 0.0000001
                    call testhermitian (check, dim_hs, ham_d, delta)
                end if

                if (allocated(eigstate)) deallocate(eigstate)
                allocate(eigstate(dim_hs, nest, ndis))

                !-----------------------------------!
                !         Diagonalization           !
                !-----------------------------------!

                print*, 'Diagonalizing Hamiltonian...'
                print*, ''
                if(full == 1) then !Full diagonalization
                    !$omp critical
                    call cfulldiag(rvec, dim_hs, ham_d, evals)
                    !call cfulldiag2(rvec, dim_hs, ham_d, evals)
                    !$omp end critical
                    cenergy(1:nev,nv+1,nv2+1,ww) = evals
                    do i = 1, nest
                        eigstate(1:dim_hs,i,ww) = ham_d(1:dim_hs,i)
                    end do
                else !Sparse diagonalization
                    call complex_diag(num_threads, dim_hs, nev, ncv, nest, mode, rvec, nnz, ham, rc, cenergy(1:nev,nv+1,nv2+1,ww), eigstate(1:dim_hs,1:nest,ww))
                end if
                print*, 'Finished Hamiltonian diagonalization.'
                print*, ''

                !-----------------------------------------------!
                !         Entanglement entropy and IPR          !
                !-----------------------------------------------!

                if (ee == 1 .or. ip == 1) then
                    if ( ee == 1 ) then
                        print*, 'Entanglement entropy ... '
                        print*, ''
                        file_name = "entanglement_entropy_"//parameters
                        open(11 + units(thread_num + 1, thread_num2 + 1), file=file_name)

                        if (allocated(entropy)) deallocate(entropy)
                        if (allocated(aventropy)) deallocate(aventropy)
                        allocate(entropy(eecut))
                        allocate(aventropy(eecut))

                        entropy   = 0.d0
                        aventropy = 0.d0
                    end if
                    if ( ip == 1 ) then
                        print*, 'IPR ... '
                        print*, ''
                        file_name = "ipr_"//parameters
                        open(61 + units(thread_num + 1, thread_num2 + 1), file=file_name)
                        if(allocated(iprv1)) deallocate(iprv1)
                        if(allocated(iprv2)) deallocate(iprv2)
                        allocate(iprv1(nest))
                        allocate(iprv2(nest))
                    end if

                    if (ti == 1) then
                        do j = 1, nest
                            if (ip == 1) then
                                if (allocated(weightsv1)) deallocate(weightsv1)
                                allocate(weightsv1(ncompv1))
                                if (allocated(weightsv2)) deallocate(weightsv2)
                                allocate(weightsv2(ncompv2))
                                weightsv1 = 0.d0
                                weightsv2 = 0.d0
                                do i = 1, dimtot
                                    if(allstates(i,2) < 0) cycle
                                    weightsv1(compv1(allstates(i, 5))) = weightsv1(compv1(allstates(i, 5))) +  abs(eigstate(allstates(i,2),j,ww))**2 * (dble(allstates(i,4)))**(-1)
                                    weightsv2(compv2(allstates(i, 5))) = weightsv2(compv2(allstates(i, 5))) +  abs(eigstate(allstates(i,2),j,ww))**2 * (dble(allstates(i,4)))**(-1)
                                end do
                                iprv1(j) = sum((weightsv1)**2)
                                iprv2(j) = sum((weightsv2)**2)
                                write( 61 + units( thread_num + 1, thread_num2 + 1 ), *) iprv1( j ), iprv2( j )
                            end if
                            if ( ee == 1 ) then
                                entropy = 0.d0
                                do la = 1, eecut
                                    thresh = 0.00001
                                    call ti_centent(num_threads, dim_hs, dimtot, la, sites - la , allstates(1:dimtot,1), eigstate(1:dim_hs,j,ww), thresh, mom, sites, allstates(1:dimtot,4), allstates(1:dimtot,3), allstates(1:dimtot,2), entropy(la))
                                    write(11 + units(thread_num + 1, thread_num2 + 1), *) entropy(la)
                                    aventropy(la) = aventropy(la) + entropy(la)
                                end do
                            end if
                        end do
                        if ( ip == 1 ) then
                            close( 61 + units( thread_num + 1, thread_num2 + 1 ) )
                            aviprv1 = sum(iprv1) / dble(size(iprv1))
                            aviprv2 = sum(iprv2) / dble(size(iprv2))
                            file_name = "avipr_"//parameters
                            open( 61 + units(thread_num + 1, thread_num2 + 1), file = file_name)
                            write( 61 + units(thread_num + 1, thread_num2 + 1), *) aviprv1, aviprv2
                            close( 61 + units(thread_num + 1, thread_num2 + 1))
                        end if
                    else if (ti == 0) then
                        do j = 1, nest
                            do la = 1, eecut
                                call factorize(dim_hs, basiss, la, bjlr) !Create factorization matrix B_jlr assigning bipartite factorization to each basis state
                                call centent(num_threads, dim_hs, la, sites - la , bjlr, eigstate(1:dim_hs,j,ww), thresh, entropy(la), singval)
                                write(11 + units(thread_num + 1, thread_num2 + 1), *) entropy(la)
                                write(81 + units(thread_num + 1, thread_num2 + 1), *) singval
                                aventropy(la) = aventropy(la) + entropy(la)
                                if (allocated(singval)) deallocate(singval)
                                if (allocated(bjlr)) deallocate(bjlr)
                            end do
                        end do
                    end if

                    if ( ee == 1 ) then
                        aventropy = aventropy / nest
                        close(11 + units(thread_num + 1, thread_num2 + 1))
                        file_name = "averaged_entanglement_entropy_"//parameters
                        open(11 + units(thread_num + 1, thread_num2 + 1), file = file_name)
                        do la = 1, eecut
                            write(11 + units(thread_num + 1, thread_num2 + 1), *) aventropy(la)
                        end do
                        close(11 + units(thread_num + 1, thread_num2 + 1))
                        print*, 'Finished calculation of entanglement entropy.'
                        print*, ''
                        if( allocated(entropy) ) deallocate(entropy)
                        if( allocated(aventropy) ) deallocate(aventropy)

                    end if

                    if (ip == 1 .and. allocated(iprv1) ) deallocate(iprv1)
                    if (ip == 1 .and. allocated(iprv2) ) deallocate(iprv2)
                end if

            end do !Disorder realization loop

            100 format(1000(F30.20))
            101 format(2000(F30.20))

            print*, 'Saving states and energies ... '
            print*, ''
            if (ti == 1 .and. bc == 'p') then
                file_name = "ti_energies_"//parameters
            else
                file_name = "energies_"//parameters
            end if

            open(11+units(thread_num+1,thread_num2+1),file=file_name)
            write(11+units(thread_num+1,thread_num2+1),*) ndis

            do j = 1, nev
                do i = 1, ndis
                    write(11+units(thread_num+1,thread_num2+1),100) real(cenergy(j,nv+1,nv2+1,i)), aimag(cenergy(j,nv+1,nv2+1,i))
                end do
            end do
            close(11+units(thread_num+1,thread_num2+1))

            if(states == 1 .and. rvec) then
                if (ti == 1 .and. bc == 'p') then
                    file_name = "ti_states_"//parameters
                else
                    file_name = "states_"//parameters
                end if
                open(11+units(thread_num+1,thread_num2+1),file=file_name)
                write(11+units(thread_num+1,thread_num2+1),*) dim_hs
                write(11+units(thread_num+1,thread_num2+1),*) ndis
                do j = 1, n_st
                    do ww = 1, ndis
                        do i = 1, dim_hs
                            write(11+units(thread_num+1,thread_num2+1),100) real(eigstate(i, j, ww)), aimag(eigstate(i, j, ww))
                        end do
                    end do
                end do
                close(11+units(thread_num+1,thread_num2+1))
            end if
            print*, 'Finished saving states and energies. '
            print*, ''
            if( allocated(eigstate) ) deallocate(eigstate)
            if( allocated(ham_d) ) deallocate(ham_d)
        end do !vloop
        !$omp end parallel do
    end do !v2loop
    !$omp end parallel do

    !--------------------------------------------!
    !         Level spacing statistics           !
    !--------------------------------------------!

    if(lvlp == 1) then
        print*, 'Level statistic parameter ... '
        print*, ''
        if(v1list == 1 .and. v2list == 1) then
            !rpar
            write (file_name,"('lvl_stat_L=',i0,'N=',i0,'t=',f12.4,'k=',i0,'.dat')") sites, pts, t, mom
            file_name = trim_name(file_name)
            open(99, file=file_name)
            do i = 1, ndv2
                do j = 1, ndv
                    call lvlpar (frac, nev, cenergy(1:nev,j,i,1), rmean)
                    write(99,*) v1_list(j), v2_list(i), rmean
                end do
            end do
            close(99)
        else if (v1list == 1 .and. v2list == 0) then
            !rpar
            v2 = 0!v2ext !SPARTAN
            write (file_name,"('lvl_stat_L=',i0,'N=',i0,'t=',f12.4,'k=',i0,'V2=',i0,'.dat')") sites, pts, t, mom, v2i
            file_name = trim_name(file_name)
            open(99, file=file_name)
            do j = 1, ndv
                call lvlpar (frac, nev, cenergy(1:nev,j,1,1), rmean)
                write(99,*) v1_list(j), rmean
            end do
            close(99)
        else if (v1list == 0 .and. v2list == 1) then
            !rpar
            vv = 0! vext !SPARTAN
            write (file_name,"('lvl_stat_L=',i0,'N=',i0,'t=',f12.4,'k=',i0,'V1=',i0,'.dat')") sites, pts, t, mom, vi
            file_name = trim_name(file_name)
            open(99, file=file_name)

            do j = 1, ndv2
                call lvlpar (frac, nev, cenergy(1:nev,j,1,1), rmean)
                write(99,*) v2_list(j), rmean
            end do
            close(99)
        end if
    end if

    if(allocated(ham)) deallocate(ham)
    if(allocated(ham_off)) deallocate(ham_off)
    if(allocated(ham_di)) deallocate(ham_di)
    if(allocated(ham_d)) deallocate(ham_d)
    if(allocated(rc)) deallocate(rc)
    if(allocated(rc_off)) deallocate(rc_off)
    if(allocated(rc_di)) deallocate(rc_di)
    if(allocated(eigstate)) deallocate(eigstate)
    if(allocated(cenergy)) deallocate(cenergy)
    if(allocated(evals)) deallocate(evals)
    if(allocated(occupation)) deallocate(occupation)
    if(allocated(work)) deallocate(work)
    if(allocated(units)) deallocate(units)
    if(allocated(units2)) deallocate(units2)
    if(allocated(basiss)) deallocate(basiss)
    if(allocated(entropy)) deallocate(entropy)
    if(allocated(aventropy)) deallocate(aventropy)

    print*, 'Finished calculation.'
    print*, ''
!    end do !Momentum loop

    116 continue

    !if (ee == 1 .and. ti == 1 .and. bc == 'p' .and. nmom > 1) then
    !    if (allocated(totavtemp)) deallocate(totavtemp)
    !    allocate (totavtemp(eecut))
    !    do nv2 = 0, ndv2-1, 1
    !        if (v2list == 0) then
    !            v2 = v2min + nv2*dv2
    !        else
    !            v2 = v2_list(nv2 + 1)
    !        end if
    !
    !        do nv = 0, ndv-1, 1
    !            if (v1list == 0) then
    !                vv = vmin + nv*dv
    !            else
    !                vv = v1_list(nv + 1)
    !            end if
    !
    !            write (file_name,"('L=',i0,'N=',i0,'t=',f12.4,'V=',f12.4, &
    !                        'V2=',f12.4,'W=',f12.4,'eps=',f12.4,'ndis=',i0,'BC=',a,'.dat')") sites,pts,t,vv, &
    !                        v2,dis,eps,ndis,bc
    !            file_name = trim_name(file_name)
    !            file_name = "tot_averaged_entanglement_entropy_"//file_name
    !            open(11, file = file_name)
    !            totavtemp = 0
    !            do nmom = 1, sites
    !                totavtemp = totavtemp + totaventropy(1:eecut, nmom, nv+1, nv2+1, 1)
    !            end do
    !            totavtemp = totavtemp / sites
    !            write(11, *) totavtemp
    !            close(11)
    !        end do
    !    end do
    !    if(allocated(totaventropy)) deallocate(totaventropy)
    !end if

    if (ti == 0 .and. bc == 'p') then
        write (parameters,"('L=',i0,'N=',i0,'t=',f12.4,'W=',f12.4,'eps=',f12.4,'ndis=',i0,'BC=',a,'.dat')") sites,pts,t,dis,eps,ndis,bc
    else
        write (parameters,"('L=',i0,'N=',i0,'k=',i0,'t=',f12.4,'W=',f12.4,'eps=',f12.4,'ndis=',i0,'BC=',a,'.dat')") sites,pts,mom,t,dis,eps,ndis,bc
    end if
    parameters=trim_name(parameters)
    call datetime(1, nthreads, parameters)
    call printparams(dim, sites, pts, ti, mom, nthreads, dim_hs, nev, ndis, v1list, v2list, &
                     t, dis, dv, vmax, vmin, dv2, v2min, v2max)

end program

!------------------------------!
!        END OF PROGRAM        !
!------------------------------!
