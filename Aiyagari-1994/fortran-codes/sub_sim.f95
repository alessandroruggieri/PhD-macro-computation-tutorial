subroutine sub_sim(V,pf_a)
    use globals
    use mathutils
    implicit none
! Declare input/output variables
    integer, dimension(1:2, 1:2, gp_a), intent(in) :: pf_a
    real, dimension(1:2, 1:2, gp_a), intent(in) :: V
! Declare internal variables
    integer :: t, j, e_aux
    integer, dimension(periods,agents) :: employed, non_part, zed, assets, state
    real, dimension(periods,agents)  :: k
    real, dimension(1:3,1:3) :: aux_trans
    real, dimension(periods,agents) :: job, shock_z
    real, dimension(gp_a) :: a_grid
print *, "Starting Simulation..."

! Define grid for assets
    a_grid = linspace(a_min, a_max, gp_a)

! Define initial values
    state = 0
    state(1,:) = 3
    assets = 0
    assets(1,:) = 1
    zed = 0
    zed(1,:) = 1
    employed = 0
    non_part = 0
    non_part(1,:) = 1
    trans = 0.0
! Generate random shocks
    call random_number(job)
    call random_number(shock_z)
!--------------------------------------- Slide 1 -----------------------------------------
do t = 2, periods
do j = 1, agents
! Update productivity shock
    if (shock_z(t,j) >= rho) then ! Change productiviy
    if (zed(t-1,j) == 1) then
        zed(t,j) = 2
    elseif (zed(t-1,j) == 2) then
        zed(t,j) = 1
    endif
    else
        zed(t,j) = zed(t-1,j)
    endif
! Determine unconstrained labor supply
    if ( V(1,zed(t,j),assets(t-1,j)) <= V(2,zed(t,j),assets(t-1,j)) ) then
        e_aux = 1
    else
        e_aux = 0
    endif
!--------------------------------------- Slide 2 -----------------------------------------
! Simulation
if (employed(t-1,j) == 0 ) then ! Not employed
if (job(t,j) <= lambda ) then   ! Finds a job
    if (e_aux == 1) then        ! Wants to work
        employed(t,j) = 1
        state(t,j)    = 1
        trans(state(t-1,j),state(t,j)) = trans(state(t-1,j),state(t,j)) + 1.0
        assets(t,j)   = pf_a(2,zed(t,j),assets(t-1,j))
        k(t,j)        = a_grid(assets(t,j))
    else
        employed(t,j) = 0
        state(t,j)    = 3
        non_part(t,j) = 1
        trans(state(t-1,j),state(t,j)) = trans(state(t-1,j),state(t,j)) + 1.0
        assets(t,j)   = pf_a(1,zed(t,j),assets(t-1,j))
        k(t,j)        = a_grid(assets(t,j))
    endif
else ! Doesn't find job
    if (e_aux == 1) then ! Wants to work
        employed(t,j) = 0
        state(t,j)    = 2
        trans(state(t-1,j),state(t,j)) = trans(state(t-1,j),state(t,j)) + 1.0
        assets(t,j)   = pf_a(1,zed(t,j),assets(t-1,j))
        k(t,j)        = a_grid(assets(t,j))
    else
        employed(t,j) = 0
        state(t,j)    = 3
        non_part(t,j) = 1
        trans(state(t-1,j),state(t,j)) = trans(state(t-1,j),state(t,j)) + 1.0
        assets(t,j)   = pf_a(1,zed(t,j),assets(t-1,j))
        k(t,j)        = a_grid(assets(t,j))
    endif
endif
!--------------------------------------- Slide 3 -----------------------------------------
elseif (employed(t-1,j) == 1 ) then ! Employed



if (job(t,j) <= 1- delta ) then ! stay in the job
    if (e_aux == 1) then ! Wants to work
        employed(t,j)  = 1
        state(t,j)     = 1
        trans(state(t-1,j),state(t,j)) = trans(state(t-1,j),state(t,j)) + 1.0
        assets(t,j)    = pf_a(2,zed(t,j),assets(t-1,j))
        k(t,j)         = a_grid(assets(t,j))
    else
        employed(t,j)  = 0
        state(t,j)     = 3
        non_part(t,j)  = 1
        trans(state(t-1,j),state(t,j)) = trans(state(t-1,j),state(t,j)) + 1.0
        assets(t,j)    = pf_a(1,zed(t,j),assets(t-1,j))
        k(t,j)         = a_grid(assets(t,j))
    endif
else ! separation
    if (e_aux == 1) then ! Wants to work
        employed(t,j)  = 0
        state(t,j)     = 2
        trans(state(t-1,j),state(t,j)) = trans(state(t-1,j),state(t,j)) + 1.0
        assets(t,j)    = pf_a(1,zed(t,j),assets(t-1,j))
        k(t,j)         = a_grid(assets(t,j))
    else
        employed(t,j) = 0
        state(t,j)    = 3
        non_part(t,j) = 1
        trans(state(t-1,j),state(t,j)) = trans(state(t-1,j),state(t,j)) + 1.0
        assets(t,j)   = pf_a(1,zed(t,j),assets(t-1,j))
        k(t,j)        = a_grid(assets(t,j))
    endif
endif

endif

enddo ! Loop over agents
enddo ! Loop over periods
!--------------------------------------- Slide 3 -----------------------------------------
do j = 1, 3
do t = 1, 3
    aux_trans(j,t) = (trans(j,t) / (sum(trans(j,:))))*100.0
enddo
enddo
trans = aux_trans
Epop  = (real(sum(employed)) / real(periods*agents))*100.0
Opop  = (real(sum(non_part)) / real(periods*agents))*100.0
kmean = (real(sum(k))        / real(periods*agents))
print *, "Done with Simulation"
end subroutine sub_sim
