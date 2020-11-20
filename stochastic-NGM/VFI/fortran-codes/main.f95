program main
    use globals
    implicit none

! Declare variables
    integer :: t1, t2, clock_rate, clock_max
    real :: time

! Start CPU clock
    call system_clock(t1, clock_rate, clock_max)

! Read parmeters and data
    call sub_read

! Neo-classical growth model
    call sub_model

! Print model output
   call sub_print

! Display employed time
    call system_clock( t2, clock_rate, clock_max)
    time = real( t2 - t1 ) / real( clock_rate )
    print *, '------------------------------'
    print *, "Time to Solve (sec.) = "
    print "(f10.2)", time/60.0d0
end program main
