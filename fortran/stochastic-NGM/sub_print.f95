subroutine sub_print
use globals
    implicit none

! Print output on screen
    print *, '------------------------------'
    print *, "Average Consumption", cmean
    print *, "Average Output", ymean
    print *, "Average Investment", invmean
    print *, "Average Capital Stock:", kmean
    print *, '------------------------------'



end subroutine
