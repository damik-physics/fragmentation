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
    integer, intent(in) :: NNZ, N
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
    integer, intent(in) :: NNZ, N
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
    integer, intent(in) :: NNZ, N
    integer, intent(in) :: indices(NNZ), indptr(N+1)
    integer, intent(out) :: labels(N)
    integer :: label, row, col, ptr, head, v
    labels = VOID
    label = 0
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