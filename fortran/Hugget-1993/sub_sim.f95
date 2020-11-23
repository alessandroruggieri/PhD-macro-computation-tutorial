subroutine sub_sim(pf_a)
    use globals
    use mathutils
    implicit none
! Declare input/output variables
    integer, dimension(1:2, gp_a), intent(in) :: pf_a
! Declare internal variables
    integer :: t, j
    integer, dimension(periods,agents) :: zed, assets, debtors
    real, dimension(periods,agents)  :: k
    real, dimension(periods,agents) :: job, shock_z
    real, dimension(gp_a) :: a_grid
print *, "Starting Simulation..."

! Define grid for assets
    a_grid = linspace(a_min, a_max, gp_a)

! Define initial values
    assets = 0
    assets(1,:) = 1
    zed = 0
    zed(1,:) = 1
    debtors = 0
    debtors(1,:) = 0

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
    k(t,j)         = a_grid(assets(t,j))

    if (k(t,j)<=0) then
        debtors(t,j)=1
    endif

enddo ! Loop over agents
enddo ! Loop over periods

!--------------------------------------- Slide 3 -----------------------------------------
kmean = (real(sum(k))        / real(periods*agents))
kstd  = sqrt((sum(k**2)-sum(k)**2/size(k))/(size(k)-1))
kneg  = (real(sum(debtors))        / real(periods*agents))

print *, "Done with Simulation"
end subroutine sub_sim
