module globals
! Declare integer variables
    integer :: gp_a, periods, agents
! Declare real variables
    real :: alfa, beta, lambda, delta, rho, tolk, tolV, a_min, a_max, alpha, dta_Epop, dta_Opop, Epop, Opop, kmean
! Declare real arrays
    real, dimension(1:2) :: z
    real, dimension(1:3, 1:3) :: dta_trans, trans
end module
