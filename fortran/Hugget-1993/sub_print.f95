subroutine sub_print
use globals
    implicit none

! Print output on screen
    print *, '------------------------------'
    print *, "Average Wealth:", kmean
    print *, "Wealth Inequality:", kstd
    print *, "Share of debtors:", kneg
    print *, '------------------------------'



end subroutine
