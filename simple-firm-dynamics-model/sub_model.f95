subroutine sub_model(p)
    use globals
    use mathutils
    implicit none

! Declare input/output variables
    real, intent(in) :: p

! Declare internal variables
    integer :: iter, ind_z, ind_a, ind_a_tw, opt_a
    integer, dimension(1:2, gp_a) :: pf_a
    real :: error_V, max_ut, revenue, val_ut
    real, dimension(gp_a) :: a_grid
    real, dimension(1:2, 1) :: gamma
    real, dimension(1:2, gp_a) :: V, V_

! Define grid for assets
    a_grid = linspace(a_min, a_max, gp_a)

! Define gamma
    gamma(1,1) = rho
    gamma(2,1) = 1.0 - rho

! Define initial values
    iter = 0 ! Number of iterations

! Initial values for VF iteration (irrelevant)
    do ind_z = 1, 2
        do ind_a = 1, gp_a
            V_(ind_z, ind_a) = 1
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

    max_ut = -1.0e6
    opt_a = 1

        do ind_a = 1, gp_a   ! Loop over asset today
            do ind_a_tw = 1, gp_a    ! Loop over asset tomorrow

                revenue = (p*z(ind_z)*(1-alfa))**(1/alfa)*a_grid(ind_a)  - ch/chc*(a_grid(ind_a_tw)-a_grid(ind_a))**chc

                val_ut = revenue + beta*(1-delta)*V(ind_z,ind_a_tw)

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



        ! Expected value at entry
        EV = 0.0
        do ind_z = 1, 2   ! Loop over z tomorrow
            EV = EV + gamma(ind_z,1)*V(ind_z,1)
        enddo  ! Close loop over assets tomorrow



! Print convergence error
print *, "Done with Value Function Iteration"
!print *, 'VF error =', error_V


end subroutine sub_model
