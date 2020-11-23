module mathutils
    implicit none
    contains
! Linspace
    function linspace(x_1, x_n, n)
        integer :: n, i
        real :: x_1, x_n, linspace(n)
        do i = 1, n ! Loop from i=1 to ind=n
            linspace(i) = x_1 + (((x_n - x_1)*(real(i-1)))/(real(n-1)))
        enddo ! close loop
    end function linspace
! Exponential utility
    function exp_ut(c, sigma)
        real :: exp_ut, c, sigma
        exp_ut = 1 - ( ( exp(-(sigma*c)) ) / sigma )
    end function exp_ut
! CRRA utility
    function crra_ut(c, sigma)
        real :: crra_ut, c, sigma
        if ( sigma == 1.0 ) then
            crra_ut = log(c)
        else
            crra_ut = (c**(1.0-sigma))/(1.0-sigma)
        endif
    end function crra_ut
end module mathutils
