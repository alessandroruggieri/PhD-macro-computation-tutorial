subroutine sub_sim(pf_a)
    use globals
    use mathutils
    implicit none
! Declare input/output variables
    integer, dimension(1:2, gp_a), intent(in) :: pf_a
! Declare internal variables
    integer :: t, j
    integer, dimension(periods,agents) :: zed, assets
    real, dimension(periods,agents)  :: k, kprime, y, c, inv
    real, dimension(periods,agents) :: job, shock_z
    real, dimension(gp_a) :: a_grid
print *, "Starting Simulation..."

! Define grid for assets
    a_grid = linspace(a_min, a_max, gp_a)

! Define initial values
    assets = 0
    assets(1,:) = 4
    zed = 0
    zed(1,:) = 1


! Generate random shocks
    call random_number(job)
    call random_number(shock_z)
!--------------------------------------- Slide 1 -----------------------------------------
! Simulation
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

    assets(t,j)    = pf_a(zed(t,j),assets(t-1,j))
    k(t,j)         = a_grid(assets(t-1,j))
    kprime(t,j)    = a_grid(assets(t,j))
    c(t,j)         = z(zed(t,j))*k(t,j)**alfa -kprime(t,j) + (1-delta)*k(t,j)
    y(t,j)         = z(zed(t,j))*k(t,j)**alfa
    inv(t,j)       = kprime(t,j) - (1-delta)*k(t,j)

enddo ! Loop over agents
enddo ! Loop over periods

kmean = (real(sum(k))        / real(periods*agents))

cmean = (real(sum(c))        / real(periods*agents))

ymean = (real(sum(y))        / real(periods*agents))

invmean = (real(sum(inv))    / real(periods*agents))

print *, "Done with Simulation"
end subroutine sub_sim
