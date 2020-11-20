subroutine sub_read
    use globals
    implicit none
    integer :: i, j ! Declare internal variables
    open(unit = 1, file = 'numerical.txt')
    read(1,*) tolk      ! Numerical tolerance asset market equilibrium
    read(1,*) tolV      ! Numerical tolerance VF iteration
    read(1,*) beta      ! Discount rate
    read(1,*) a_min     ! Minimum assets
    read(1,*) a_max     ! Maximum assets
    read(1,*) gp_a      ! Grid points for assets
    read(1,*) z(1)      ! Low productivity
    read(1,*) z(2)      ! High productivity
    read(1,*) rho       ! Persistency
    read(1,*) periods   ! Number of periods
    read(1,*) agents    ! Number of agents
    close(1)
end subroutine sub_read
