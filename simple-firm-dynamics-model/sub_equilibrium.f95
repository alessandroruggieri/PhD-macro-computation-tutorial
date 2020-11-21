subroutine sub_equilibrium
    use globals
    use mathutils
    implicit none

    ! Declare internal variables
    integer :: iter
    real :: p, pmin, pmax, error_k

    ! Initialize variables
    iter    = 0     ! Number of iterations
    error_k = 1.0e4 ! Convergence criteria

    ! Bisection
    pmin    = 0  ! Lower bound interest rate
    pmax    = 10 ! Upper bound interest rate

    ! Iterate over the interest rate to clear the asset market
    do
        p = (pmax+pmin)/2            ! Update interest rate
        iter = iter + 1              ! Update # iterations
    if ( error_k <= tolk ) exit      ! Convergence is achived
    if ( iter > 1e5 ) then           ! Protect against no convergence
        print *, "Too many iterations to solve the equilibrium"
        exit
    endif

    ! Value function iteration
    call sub_model(p)


    ! Converenge criterium
    error_k = abs(EV-centry)
    print *, 'Equilibrium error =', error_k
    print *, 'Equilibrium price =', p

    ! Update boundaries
    if (EV<centry) then
        pmin=p
    else
        pmax=p
    endif

    enddo


end subroutine sub_equilibrium
