!****************************************************************************
!
!  Program: JAPP - Just A Poor Particle
!
!  Purpose: Not for any specific purpose.
!    
!       There is a poor particle names Brown. He borns in this code, 
!       and dies with the progress. Unfortunately, both you and he have no 
!       idea with his next movement which should be decided randomly.
!       Possibly, the Brown will travel over a distance if he runs enough time.
!       
!       Now, start the Brown's traveling.
! 
!       The author of this code is: Ge Zhang (E-mail: ge_zhang_phys@gmail.com)
!
!****************************************************************************

    program JAPP
    implicit none

    ! Variables
    ! The position, force, velocity, displacement
    ! Mass, Radi, Gravity, Time, Timestep

    integer start, loop, dims
    real*8 Poi(1:3), Poi0(1:3), Forc(1:3), Vel(1:3)
    real*8 Grav(1:3)
    real*8 Mass, Rho, Radi
    real*8 Time, TStep
    real*8 pi, kB
    real*8 FldRho, FldVis, FldTem, Fldfc, Varia
    real*8 FldByc(1:3)
    real*8 Randomdis
    character*256 filename01
    parameter (pi=3.141592654)
    parameter (kB=1.38064852E-23) ! Boltzmann constant, m2.kg.s-2.K-1
    external Randomdis

    ! /Particle Infor/
    common Mass, Rho, Radi, Grav, Time, TStep, &
            & Poi, Poi0, Forc, Vel, Fldfc


    
9   print*, '-----  Input number "0" to start -----'
    read*, start
    
    if (start /= 0) then 
        goto 9
    endif
    

    ! Initializing
    ! Length - meter, Time - second, Weight - kg
    Time    =   0.d0        ! s
    TStep   =   0.001d0    ! s
    ! The momenta relaxtion assumes: TStep  >>  m*D/(kB*T) = m/(6*pi*radi*miu), approximately 10E-9
    
    Radi    =   5E-6     ! m
    Rho     =   1.0E3     ! kg/m^3
    Mass    =   4.d0*pi*(Radi**3)*Rho/3.d0  ! kg
    
    Poi(1:3)    =   0.d0
    Poi0(1:3)   =   0.d0
    Forc(1:3)   =   0.d0


    Grav(1:2)   =   0.d0
    Grav(3)     =  -9.8d0       ! m/s^2

    ! Fluid properties
    FldRho      =  1.0E3        ! kg/m^3
    FldVis      =  8.9E-4       ! kg/(m*s)
    FldTem      =  25.0+273.15  ! Celsius to Kelvin, K

    FldByc(1)   =  0.d0
    FldByc(1)   =  0.d0
    FldByc(3)   = -4.d0*pi*(Radi**3)*FldRho*Grav(3)/3.d0

    Fldfc       =  (kB*FldTem)/(6*pi*FldVis*Radi)     ! Fluid coeff tensor for i=j

    Varia       =  2.d0*Fldfc*TStep

    call random_seed()

    write(filename01,*) Time
    open(1001,file=''//trim(adjustl(filename01))//'.DAT',status='unknown')

        write(1001,*)  'Title="Brownian Particle"'
        write(1001,*)  'variables="X","Y","Z","Vx","Vy","Vz"'
        write(1001,*)  'zone   ',' F=point'
        write(1001,19,advance="no") Poi(1:3), Vel(1:3)
            
     close(1001)

    ! Main loop
    do 99 loop = 1, 10000, 1

    do 99 dims = 1, 3, 1

        Forc(dims)     =    Grav(dims)*Mass + FldByc(dims)
        
        Poi0(dims)     =    Poi(dims)
       
        Poi(dims)      =    Poi0(dims)  +  (Fldfc*Forc(dims)*TStep)/(kB*FldTem)   +   Randomdis(Varia)

        Vel(dims)      =    (Poi(dims)   -   Poi0(dims))/TStep

        if ( mod(loop,50)==0 ) then
            !** The trajectory of the particle/position, velocity with time **
            Time    =   TStep*(loop*1.d0)

            write(filename01,*) Time
            open(1001,file=''//trim(adjustl(filename01))//'.DAT',status='unknown')

              write(1001,*)  'Title="Brownian Particle"'
              write(1001,*)  'variables="X","Y","Z","Vx","Vy","Vz"'
              write(1001,*)  'zone   ',' F=point'
              write(1001,19,advance="no") Poi(1:3), Vel(1:3)
            
            close(1001)

        19  format(6E20.10,/)

        endif

99  continue
    
    print *, 'The Travel Ends'
    
    end program JAPP
 
    