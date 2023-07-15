module module_variables_fragmentation

implicit none

integer, parameter :: nthreads = 1 !SPARTAN (COMMENT OUT)

integer, parameter :: ee   = 1
integer, parameter :: ip   = 1
integer, parameter :: lvlp = 1

integer, parameter :: exact       = 1
integer, parameter :: printing    = 1
integer, parameter :: nn_int      = 1
integer, parameter :: nnn_int     = 1
integer, parameter :: disorder    = 1
integer, parameter :: spinless    = 1
integer, parameter :: states      = 1
integer, parameter :: shiftinvert = 0
integer, parameter :: dim         = 1
integer, parameter :: g_fact      = 2
integer, parameter :: ndis        = 1
integer, parameter :: frac        = 3
integer, parameter :: nevext      = 10
integer, parameter :: nevhf       = 10
integer, parameter :: ncv0        = 120
integer, parameter :: dimthresh   = 250
integer, parameter :: n_st        = 5
!integer, parameter :: ti = 1

double precision, parameter :: filling = 0.5

double precision, parameter :: dv    = 1.0d0
double precision, parameter :: dv2   = 1.0d0
double precision, parameter :: vmin  = 0.0d0
double precision, parameter :: vmax  = 0.0d0
double precision, parameter :: v2min = 0.0d0
double precision, parameter :: v2max = 0.0d0

double precision, parameter :: dis   = 0.0d0
double precision, parameter :: shift = nex/4 - 0.5
double precision, parameter :: t     = 1.0d0
double precision, parameter :: one   = 1.0D+0

logical, parameter :: dynamic = .False.
logical, parameter :: nested  = .True.
!logical, parameter :: rvec    = .True.

character, parameter :: bc*1     = 'p' !'p'
character, parameter :: switch*2 = 'SA' !'LA'

external :: dsyev, dsaupd, dseupd, dmout

!double precision, parameter :: pi = 4*atan(1.d0)

end module module_variables_fragmentation
