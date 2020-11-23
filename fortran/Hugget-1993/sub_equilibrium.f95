subroutine sub_equilibrium
    use globals
    use mathutils
    implicit none

    ! Declare internal variables
    integer :: iter
    real :: r, rmin, rmax, error_k

    ! Initialize variables
    iter    = 0     ! Number of iterations
    error_k = 1.0e4 ! Convergence criteria

    ! Bisection
    rmin    = -0.1  ! Lower bound interest rate
    rmax    = 0.2   ! Upper bound interest rate

    ! Iterate over the interest rate to clear the asset market
    do
        r = (rmax+rmin)/2            ! Update interest rate
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
    print *, 'Equilibrium error =', error_k

    ! Update boundaries
    if (kmean<0) then
        rmin=r
    else
        rmax=r
    endif

    enddo


end subroutine sub_equilibrium
