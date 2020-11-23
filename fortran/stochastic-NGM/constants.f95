module constants
    implicit none
    save

    integer, parameter:: dbl = kind(1d0)
    real(dbl), parameter:: pi = 3.14159265358979d0
    character(1), parameter:: tab = char(9)
    character(1), parameter:: ret = char(13)

end module constants
