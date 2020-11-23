subroutine sub_read
    use globals
    implicit none
    !integer :: i, j ! Declare internal variables
    open(unit = 1, file = 'numerical.txt')
    read(1,*) tolk      ! Numerical tolerance asset market equilibrium
    read(1,*) tolV      ! Numerical tolerance VF iteration
    read(1,*) beta      ! Discount rate
    read(1,*) a_min     ! Minimum assets
    read(1,*) a_max     ! Maximum assets
    read(1,*) gp_a      ! Grid points for assets
    read(1,*) z(1)      ! Low productivity
    read(1,*) z(2)      ! High productivity
    read(1,*) periods   ! Number of periods
    read(1,*) agents    ! Number of agents
    close(1)
    open(unit = 2, file = 'parameters.txt')
    read(2,*) delta     ! Destruction probability
    read(2,*) lambda    ! Finding probability
    read(2,*) alpha     ! Working disutility
    read(2,*) rho       ! Persistence parameter
    close(2)
    !open(unit = 3, file = "data.txt")
    !read(3,*) dta_Epop
    !read(3,*) dta_Opop
    !do i = 1, 3
    !    do j = 1, 3
    !        read(3,*) dta_trans(i,j)
    !    enddo
    !enddo
    !close(3)
end subroutine sub_read
