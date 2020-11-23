subroutine sub_model(r)
    use globals
    use mathutils
    implicit none

! Declare input/output variables
    real, intent(in) :: r

! Declare internal variables
    integer :: iter, ind_z, ind_a, ind_a_tw, ind_z_tw, opt_a
    integer, dimension(1:2, 1:2, gp_a) :: pf_a
    real :: error_V, max_ut, c, val_ut
    real, dimension(gp_a) :: a_grid, exp_N, exp_J
    real, dimension(1:2, 1:2) :: gamma
    real, dimension(1:2, 1:2, gp_a) :: V, V_

! Define grid for assets
    a_grid = linspace(a_min, a_max, gp_a)

! Define gamma
    gamma(1,1) = rho
    gamma(1,2) = 1.0 - rho
    gamma(2,1) = 1.0 - rho
    gamma(2,2) = rho

! Define initial values
    iter = 0 ! Number of iterations

! Initial values for VF iteration (irrelevant)
    do ind_z = 1, 2
        do ind_a = 1, gp_a
            V_(1, ind_z, ind_a) = log((1.0/2.0)*z(ind_z)) ! Not working
            V_(2, ind_z, ind_a) = log(z(ind_z))           ! Working
        enddo
    enddo
    error_V = 1.0e4 ! V-V_

!--------------------------------------- Slide 1 -----------------------------------------
! Value function iteration
print *, "Start with Value Function Iteration..."

do
    V = V_ ! Update guess
    iter = iter + 1 ! Update # iterations
    if ( error_V <= tolV ) exit ! Convergence is achived
        if ( iter > 1e5 )  then ! Protect against no convergence
            print *, 'Too many iterations to solve value function'
            exit
        endif

    do ind_z = 1, 2

    ! Compute expected values
    exp_N = 0.0
    exp_J = 0.0

    do ind_a_tw = 1, gp_a

    ! No job offer tomorrow
    exp_N(ind_a_tw) = ( gamma(ind_z,1)*V(1,1,ind_a_tw) ) + ( gamma(ind_z,2)*V(1,2,ind_a_tw) )

    ! Job offer tomorrow
    do ind_z_tw = 1, 2
    if ( V(1,ind_z_tw,ind_a_tw) <= V(2,ind_z_tw,ind_a_tw)  ) then
    ! Work
        exp_J(ind_a_tw) = exp_J(ind_a_tw) + ( gamma(ind_z,ind_z_tw)*V(2,ind_z_tw,ind_a_tw) )
    else
    ! No work
        exp_J(ind_a_tw) = exp_J(ind_a_tw) + ( gamma(ind_z,ind_z_tw)*V(1,ind_z_tw,ind_a_tw) )
    endif ! Labor supply decision

    enddo ! Loop over z tomorrow
    enddo ! Loop over assets tomorrow

!--------------------------------------- Slide 2 -----------------------------------------
    do ind_a = 1, gp_a

    ! Non-employed agent
    max_ut = -1.0e6
    opt_a = 1

    do ind_a_tw = 1, gp_a
        c = (1.0+r)*a_grid(ind_a) - a_grid(ind_a_tw)
    if (c > 0.0) then
        val_ut = log(c) + beta*( (1.0-lambda)*exp_N(ind_a_tw) + lambda*exp_J(ind_a_tw) )
    else
        val_ut = -1.0e6
    endif
    if (val_ut >= max_ut) then
        max_ut = val_ut
        opt_a = ind_a_tw
    endif
    enddo ! Loop over a tomorrow

    V_(1,ind_z,ind_a) = max_ut
    pf_a(1,ind_z,ind_a) = opt_a
!--------------------------------------- Slide 3 -----------------------------------------
    ! Employed agent
    max_ut = -1.0e6
    opt_a = 1

    do ind_a_tw = 1, gp_a
        c = (1.0+r)*a_grid(ind_a) + z(ind_z) - a_grid(ind_a_tw)
    if (c > 0.0) then
        val_ut = log(c) - alpha + beta*( delta*exp_N(ind_a_tw) + (1.0-delta)*exp_J(ind_a_tw) )
    else
        val_ut = -1.0e6
    endif
    if (val_ut >= max_ut) then
        max_ut = val_ut
        opt_a = ind_a_tw
        endif
    enddo ! Loop over a tomorrow

    V_(2,ind_z,ind_a) = max_ut
    pf_a(2,ind_z,ind_a) = opt_a
    enddo ! Loop over a today
    enddo ! Loop over z today

    error_V = sum(abs(V-V_))
    if (mod(iter, 250) .eq. 0) then
            print *, 'VF error =', error_V
        endif
enddo ! Value function iteration
print *, "Done with Value Function Iteration"
print *, 'VF error =', error_V

! Call subroutine to simulate
call sub_sim(V,pf_a)

end subroutine sub_model
!--------------------------------------- Slide 4 -----------------------------------------
