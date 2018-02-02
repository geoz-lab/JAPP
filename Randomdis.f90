    real*8 function Randomdis(Va)
        real*8 va
        real*8 r1, r2
                
        call random_number(r1)
        call random_number(r2)
        Randomdis = (va**0.5d0)* (( -2.d0*log(r1) )**(0.5d0)) * &
                    & cos(2.d0*3.141592654d0*r2)

    return
    end function