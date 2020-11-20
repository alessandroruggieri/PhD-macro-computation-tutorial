subroutine sub_print
use globals
    implicit none

    real:: Upop, dta_Upop, d_E, d_U, d_O, avg_d, abs_d

    Upop=100.0-Epop-Opop
    dta_Upop=100.0-dta_Epop-dta_Opop


    d_E=(dta_Epop-Epop)/dta_Epop
    d_U=(dta_Upop-Upop)/dta_Upop
    d_O=(dta_Opop-Opop)/dta_Opop
    avg_d= (d_E+d_U+d_O)/3
    abs_d= abs(avg_d)

! Print output in txt
    open(unit=1, file="comparison.txt")
    write(1,*) d_E
    write(1,*) d_U
    write(1,*) d_O
    close(1)

! Print output in txt
    open(unit=2, file="moments.txt")
    write(2,*) Epop
    write(2,*) Upop
    write(2,*) Opop
    close(2)

! Print output on screen
    print *, "Employment", Epop
    print *, "Unemployment", Upop
    print *, "Inactive", Opop

! Print error
    open(unit=3, file="error.txt")
    write(3,*) abs_d
    close(3)

end subroutine
