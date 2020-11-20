subroutine sub_model(r)
    use globals
    use mathutils
    implicit none

! Declare input/output variables
    real, intent(in) :: r

! Declare internal variables
    integer :: iter, ind_z, ind_a, ind_a_tw, ind_z_tw, opt_a
    integer, dimension(1:2, gp_a) :: pf_a
    real :: error_V, max_ut, c, val_ut
    real, dimension(gp_a) :: a_grid, exp_V
    real, dimension(1:2, 1:2) :: gamma
    real, dimension(1:2, gp_a) :: V, V_

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
            V_(ind_z, ind_a) = log((1.0/2.0)*z(ind_z)) ! Not working
            V_(ind_z, ind_a) = log(z(ind_z))           ! Working
        enddo
    enddo
    error_V = 1.0e4 ! V-V_

! Value function iteration
print *, "Start with Value Function Iteration..."

do ! Start VFI loop

    V = V_   ! Update guess
    iter = iter + 1  ! Update # iterations
    if ( error_V <= tolV ) exit ! Convergence is achived
        if ( iter > 1e5 )  then ! Protect against no convergence
            print *, 'Too many iterations to solve value function'
            exit
        endif


    do ind_z = 1, 2  ! Loop over z today

        ! Compute expected values
        exp_V = 0.0
        do ind_a_tw = 1, gp_a  ! Loop over assets tomorrow
            do ind_z_tw = 1, 2   ! Loop over z tomorrow
                exp_V(ind_a_tw) = exp_V(ind_a_tw) + ( gamma(ind_z,ind_z_tw)*V(ind_z_tw,ind_a_tw) )
            enddo  ! Close loop over assets tomorrow
        enddo        ! Close loop over z tomorrow

    max_ut = -1.0e6
    opt_a = 1

        do ind_a = 1, gp_a   ! Loop over asset today
            do ind_a_tw = 1, gp_a    ! Loop over asset tomorrow
                c = (1.0+r)*a_grid(ind_a) + z(ind_z) - a_grid(ind_a_tw)
                if (c > 0.0) then
                    val_ut = log(c) + beta*exp_V(ind_a_tw)
                else
                    val_ut = -1.0e6
                endif
                if (val_ut >= max_ut) then
                    max_ut = val_ut
                    opt_a = ind_a_tw
                endif
            enddo ! Loop over a tomorrow

        V_(ind_z,ind_a) = max_ut    ! Store value function
        pf_a(ind_z,ind_a) = opt_a   ! Store policy function
        enddo  ! Close loop over a today
    enddo ! Close loop over z today

    ! Convergence flag
    error_V = sum(abs(V-V_))
    !if (mod(iter, 250) .eq. 0) then
    !       print *, 'VF error =', error_V
    !   endif

enddo ! Close VFI loop

! Print convergence error
print *, "Done with Value Function Iteration"
!print *, 'VF error =', error_V

! Call subroutine to simulate
call sub_sim(pf_a)

end subroutine sub_model
