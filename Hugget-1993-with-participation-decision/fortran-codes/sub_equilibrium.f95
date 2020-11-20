subroutine sub_equilibrium
    use globals
    use mathutils
    implicit none

    ! Declare internal variables
    integer :: iter
    real :: r, rmin, rmax, error_k

    ! Initialize variables
    rmin    = -0.2
    rmax    = 0.2
    error_k = 1.0e4

    ! Iterate over the interest rate to clear the asset market
    do
        r = (rmax+rmin)/2            ! Update guess
        print *, 'Interest rate =', r
        print *, '------------------------------'
        iter = iter + 1              ! Update # iterations
    if ( error_k <= tolk ) exit      ! Convergence is achived
    if ( iter > 1e5 ) then           ! Protect against no convergence
        print *, "Too many iterations to solve the equilibrium"
        exit
    endif

    ! Value function iteration
    call sub_model(r)


    ! Converenge criterium
    error_k = abs(kmean)
    print *, '------------------------------'
    print *, 'Capital =', kmean

    ! Update boundaries
    if (kmean<0) then
        rmin=r
    else
        rmax=r
    endif

    enddo


end subroutine sub_equilibrium
