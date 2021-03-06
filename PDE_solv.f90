program main
use ziggurat
implicit none
character :: model !Update model A or B
integer :: i, i_interval, i_z, j_z, j, N_tot, N_cnf, i_cnf, j_cnf, i_time, fnsh_time, obs_time, dbg_count = 1
real(kind=8) :: mean_rho, dot_prod, s_new, L_box, eps, dz,dz2, str_len ! ,conv_eps !square of conv criteria
real(kind=8) , dimension(:), allocatable :: rho, rho_new, mu, s_par, dist_cnf
real(kind=8) , dimension(:,:), allocatable :: tang, poly_coeff, rho_tot, mu_tot, mu_perp
!rho_tot is the matrix with all rho as a function of the string
logical :: conver = .False., PBC = .True.! Converged to solution, Periodic
                                        !Boundary Condition True or False

!Read input
call read_input(L_box,eps,N_tot,fnsh_time,obs_time, model, dz, N_cnf, mean_rho)
dz2 = dz**2

allocate(rho(N_tot),rho_new(N_tot), mu(N_tot), rho_tot(N_tot,N_cnf), s_par(N_cnf),dist_cnf(N_cnf-1))
allocate(poly_coeff(4,N_cnf-1),mu_tot(N_tot,N_cnf),mu_perp(N_tot,N_cnf),tang(N_tot,N_cnf))

!open files to save data
open (unit = 73, file = "rho_vs_z.dat")
open (unit = 74, file = "F_vs_s.dat")
open (unit = 75, file = "N_part_vs_time.dat")
open (unit = 76, file = "Finterface_vs_time.dat")

call init_rand_seed()

!Make initial guess of rho for each point in string
do i_cnf=1,N_cnf
    if(i_cnf.le.N_cnf/2) then
        call init_guess(rho_tot(:,i_cnf),N_tot,2)   ! Last input is case:
    else                                            ! 1 is ramp.           
        call init_guess(rho_tot(:,i_cnf),N_tot,3)   ! 2 is random near mean value defined inside routine
                                                    ! 3 is near equilibrium
                                                    ! 4 is read config from file        
    end if
    !DEBUG                                   
    !print *,"Va Hasta aca!",i_cnf           
    !do i_z = 1, N_tot
    !    write(i_cnf+100,*) dz*float(i_z), rho_tot(i_z,i_cnf)   
 
    !end do
    !write(i_cnf+100,*)

end do                                   

print*, " Starting Iteratrive approach to Solution"

do i_time = 1, fnsh_time !while (.not.conver)
    do i_cnf = 1, N_cnf ! loop over the string
        rho(:) = rho_tot(:,i_cnf)
        mu(:) = mu_tot(:,i_cnf)
        !Update solution to diffretian equation
        call update(rho,rho_new,mu,N_tot,eps, dz2, PBC, model)

        !Meassure Free Energy, N_particles, F interface, etc .
        if ( mod(i_time,obs_time) .eq. 0 ) then
            call observation(rho_new,mu,N_tot, L_box, i_time,dz,PBC)
            !do j_z = 1, N_tot
            !    write(i_cnf+100,*) dz*float(j_z),rho_new(j_z)                   
            !end do
            !write(i_cnf+100,*) 
        end if
        
        !Evaluate convergence !OLD
        !call eval_conv(rho,rho_new,N_tot,conver,conv_eps)    

        ! move New to old variables
        rho_tot(:,i_cnf) = rho_new(:)
        mu_tot(:,i_cnf) = mu(:)
    end do
     
    !Calculate Distances between  rho(:,i) and  rho(:,i+1)
    str_len = 0.
    do i_cnf = 1, N_cnf-1 ! loop over the string
        dist_cnf(i_cnf) = sqrt( sum( ( rho_tot(:,i_cnf+1) - rho_tot(:,i_cnf) )**2 ) )
        if ( dist_cnf(i_cnf) .le. 10E-6) then
            dist_cnf(i_cnf) = 10E-6
        end if
        !write(199,*) i_cnf,dist_cnf(i_cnf)
        str_len = str_len + dist_cnf(i_cnf)
    end do
    !write(199,*) 
    
    !Calculate normalized distance between configurations
    s_par(1) = 0
    do i_cnf = 1, N_cnf-1 ! loop over the string
       s_par(i_cnf+1) = s_par(i_cnf) +  dist_cnf(i_cnf)  / str_len
       !write(200,*) i_cnf,s_par(i_cnf)
    end do    
    !write(200,*) i_cnf,s_par(i_cnf)
    !write(200,*)
    
     
    ! Get third order spline for each z0 f(s|z0) = rho(z0,s)
    ! from f(s|z0) get new rho(s*N_tot+1,z0)
    do i_z = 1, N_tot  ! loop over z positions
        call spline(s_par, rho_tot(i_z,:), N_cnf, poly_coeff) !splines for this z value
         
        !get rho_tot(i_z,:) from splines
        do i_cnf = 2, N_cnf-1
            s_new = float(i_cnf-1) / float(N_cnf-1)

            !find i_interval for this s_new value
            do j_cnf = 1, N_cnf-1
                if( ( s_new.ge.s_par(j_cnf) ) .and. (s_new.le.s_par(j_cnf+1))) then
                    i_interval = j_cnf
                    exit !get out of loop
                end if
            end do
             
            !!DEBUG     
            !print*,"s_new",s_new
            !print*,"s_par(j_cnf), i_interval, s_par(j_cnf+1)",s_par(j_cnf),s_new, s_par(j_cnf+1)
            !print*,"" 

            !evaluate rho_tot(i_z,i_cnf) = poly_coeff(1,i_interval) + poly_coeff(2,i_interval) * (s_new-s_par(i_interval)) + poly_coeff(3,i_interval) * (s_new-s_par(i_interval))**2 + poly_coeff(4,i_interval) * (s_new-s_par(i_interval))**3
            rho_tot(i_z,i_cnf) = poly_coeff(1,i_interval) + poly_coeff(2,i_interval) * (s_new-s_par(i_interval)) + poly_coeff(3,i_interval) * (s_new-s_par(i_interval))**2 + poly_coeff(4,i_interval) * (s_new-s_par(i_interval))**3

            ! Get tangent vector tang(i_z,i_cnf) = d (rho_tot) / ds
            tang(i_z,i_cnf) = poly_coeff(2,i_interval) + 2. * poly_coeff(3,i_interval) * (s_new-s_par(i_interval)) + 3. * poly_coeff(4,i_interval) * (s_new-s_par(i_interval))**2

            
        end do !s_loop
   
    end do!i_z loop

    !!Check that the number of particles remains constant i.e. sum ( m(z) ) = cte
    do i_cnf = 1, N_cnf   
        rho_tot(:,i_cnf) = rho_tot(:,i_cnf) - sum( rho_tot(:,i_cnf) ) / float(N_tot) + mean_rho
    end do   

    !Get mu_perp 
    do i_cnf=2, N_cnf-1
        dot_prod = sum( tang(:,i_cnf) * mu_tot(:,i_cnf)  ) / sum(tang(:,i_cnf)**2 )
        do i_z=1, N_tot
            mu_perp(i_z,i_cnf) = mu_tot(i_z,i_cnf) - dot_prod * tang(i_z,i_cnf)
        end do
    end do
    
    write(305,*) i_time, maxval(  mu_perp(:,:) )
end do !time loop


!Save result
call write_rho(rho, N_tot, L_box)

call write_results(rho_tot, N_tot, L_box, N_cnf)
!

!Close file with rho vs L
close (73)
close (74)
close (75)
close (76)

print*, " Done!"

end program

subroutine write_results(rho_tot, N_tot, L_box, N_cnf)
implicit none
integer, intent(in) :: N_tot, N_cnf
real(kind=8),intent(in) :: rho_tot(N_tot,N_cnf), L_box
integer :: i_z, i_cnf
real(kind=8) :: Free_Energy(N_cnf), D_rho(N_cnf), dz 

dz = L_box / float(N_tot )

!Save Free Energy
Free_Energy(:) = 0.
do i_cnf = 1, N_cnf
    do i_z =1 ,N_tot-1
        D_rho(i_cnf) = ( rho_tot(i_z+1,i_cnf) - rho_tot(i_z,i_cnf) ) / dz
        Free_Energy(i_cnf) =  Free_Energy(i_cnf) - rho_tot(i_z,i_cnf)**2 + .5*rho_tot(i_z,i_cnf)**4 + D_rho(i_cnf)**2
    end do

    !Account for PBC
!    if(PBC) then
        D_rho(i_cnf) = ( rho_tot(1,i_cnf) - rho_tot(N_tot,i_cnf) ) / dz
     
        Free_Energy(i_cnf) =  Free_Energy(i_cnf) - rho_tot(N_tot,i_cnf)**2 + .5*rho_tot(N_tot,i_cnf)**4 + D_rho(i_cnf)**2
!    end if
    Free_Energy(i_cnf) = .5 * Free_Energy(i_cnf) * dz
end do

open(unit = 60, file = "F_vs_s.dat", status= "unknown")

do i_cnf = 1, N_cnf
    write(60,*) float(i_cnf-1)/float(N_cnf-1),Free_Energy(i_cnf)
end do

!Save last configurations

do i_cnf = 1, N_cnf
    do i_z = 1, N_tot
    write(100+i_cnf,*) float(i_z-1)*dz,rho_tot(i_z,i_cnf) 
    end do
end do

end subroutine write_results

subroutine observation(rho,mu,N_tot, L_box, i_time,dz,PBC)
real(kind=8) , dimension(N_tot), intent(in) :: rho, mu
integer, intent(in) :: N_tot, i_time
real(kind=8), intent(in) :: dz, L_box
logical,intent(in) :: PBC
real(kind=8) :: Free_Energy=0., N_part=0., D_rho, F_inter=0.
integer :: i

N_part=0
Free_Energy=0
F_inter=0

do i=1, N_tot-1
    D_rho = ( rho(i+1) - rho(i) ) / dz !derivative of rho
    F_inter = F_inter + D_rho**2 
    N_part = N_part + rho(i) 
    Free_Energy =  Free_Energy - rho(i)**2 + .5*rho(i)**4 + D_rho**2
end do

N_part = N_part + rho(N_tot)

!Account for PBC
if(PBC) then
    D_rho = ( rho(1) - rho(N_tot) ) / dz
    F_inter = F_inter + D_rho**2
    Free_Energy =  Free_Energy - rho(N_tot)**2 + .5*rho(N_tot)**4 + D_rho**2
end if

!Multiply by dz because I'm doing an integral
N_part = N_part * dz
F_inter = .5 * F_inter * dz 
Free_Energy = .5 * Free_Energy * dz


!Write N_part and Free_Energy to file
write(74,*) i_time, Free_Energy
write(75,*) i_time, N_part
write(76,*) i_time, F_inter
!print*, "in obs",i_time

call write_rho(rho, N_tot, L_box)

end subroutine

subroutine write_rho(rho, N_tot, L_box)
real(kind=8) , dimension(N_tot), intent(in) :: rho
integer, intent(in) :: N_tot
real(kind=8) :: z
real(kind=8), intent(in) :: L_box
integer :: i

do i = 1, N_tot
    z = float(i-1)/float(N_tot-1)*L_box 
    write(73,*) z, rho(i)
end do

write(73,*) ""

end subroutine
 
subroutine eval_conv(rho, rho_new, N_tot, conver, conv_eps)
integer :: i 
real(kind=8), intent(in) :: conv_eps
real(kind=8) , dimension(N_tot), intent(in) :: rho, rho_new
real(kind=8) , dimension(N_tot) ::  diff
integer, intent(in) :: N_tot
logical, intent(inout) :: conver
real(kind=8) :: max_diff_2

diff(:) = (rho_new(:) - rho(:) ) ** 2

max_diff_2 = maxval(diff)

if(max_diff_2.le.conv_eps) then !Apply convergence criteria
    conver = .True.
    print*, " Converged to Solution"
    !print*,max_diff
end if

end subroutine

subroutine init_guess(rho,N_tot,ini_case)
use ziggurat 
integer :: i, N_int
real(kind=8) :: inv_Ntot_1, amp=.05, mean = -0.8, dummy
integer, intent(in) :: N_tot,ini_case
real(kind=8) , dimension(N_tot), intent(inout) :: rho
select case(ini_case)
case(1) ! rho initiates as linear function 
    inv_Ntot_1 = 1. / float((N_tot-1))
    
    do i=2, N_tot-1
        rho(i) = 2. * inv_Ntot_1 *(float(i)-1.) - 1.
    end do
    
    !Boundary conditions
    rho(1) = -1.
    rho(N_tot) = 1.
  
case(2) ! rho initiates randomly around mean
    amp=.05
        mean = -0.8

    do i=1, N_tot
         rho(i) = amp*rnor()
    end do
    
    rho(:) = rho(:) - sum(rho(:)) / float(N_tot) + mean !Set mean to mean

case(3) ! rho initiates with 2 interfaces
    !Get interface position acording to mean value selected.
    !Other interface is at the beginning of simulation box z=-L
    N_int = int( (1.-mean) / 2. * float(N_tot) )
    !print*, "N_int,N_tot",N_int,N_tot 
    do i=1, N_int-1
        rho(i) = -1.
    end do
       
    do i=N_int,N_tot
        rho(i) = 1.
    end do

case(4) ! read rho from file
    open(unit = 89, file = "input_rho_vs_z.dat", status= "old") 

    !Read data
    do i=1,N_tot
        read(89,*) dummy, rho(i)
    end do

    close(89)

end select
!print*, " Initial Guess done"

end subroutine

subroutine init_rand_seed()
use ziggurat
logical :: file_ex
integer :: rand_seed
inquire(file='random_seed.dat',exist=file_ex)
if (file_ex) then
    open(unit = 68, file = "random_seed.dat", status= "old")
    read(68,*) rand_seed
    close(68)
else
    rand_seed = 465178
end if

call zigset( rand_seed )

open(unit = 68, file = "random_seed.dat", status= "unknown")
write(68,*) shr3()   
close(68)
   
end subroutine 

subroutine update(rho, rho_new, mu, N_tot, eps, dz2, PBC, model) 
integer :: i
logical, intent(in) :: PBC
integer, intent(in) :: N_tot
character, intent(in) :: model
real(kind=8), intent(in) :: eps, dz2
real(kind=8) , dimension(N_tot), intent(inout) :: rho, rho_new, mu
real(kind=8) , dimension(N_tot) :: delta_rho
real(kind=8) ::  mu_exc

!Update Boundaries
select case(PBC)
case(.False.) ! Fixed Boundary conditions
    mu(1) = 0
    mu(N_tot) = 0
    rho_new(1) = rho(1)
    rho_new(N_tot) = rho(N_tot)
    
case(.True.) ! Perioduc Boundary Conditions
    mu(1) = -rho(1) + rho(1)**3 -( rho(2) + rho(N_tot) - 2 * rho(1)) / dz2
    mu(N_tot) = -rho(N_tot) + rho(N_tot)**3 -( rho(1) + rho(N_tot-1) - 2 * rho(N_tot)) / dz2

end select

!Update mu bulk
do i = 2, N_tot-1 ! Dont touch boundaries
    mu(i) = -rho(i) + rho(i)**3 -( rho(i+1) + rho(i-1) - 2 * rho(i)) / dz2
end do

!Normlize mu: set sum(mu(:)) = 0.
!This ensures that integral over space of rho ( # particles) remains constant
mu_exc = sum(mu(:)) / float(N_tot)
mu(:) = mu(:) - mu_exc
!write(55,*) sum(mu(:))

! Calculate new densities
select case(model)
case("A")
    rho_new(:) = rho(:) - eps * mu(:)

case("B")
    do i = 2, N_tot-1
        delta_rho(i) = ( mu(i+1) + mu(i-1) - 2.*mu(i) ) / dz2
    end do
    
    if(PBC) then!If PBC are ON
        delta_rho(1) = (mu(2) + mu(N_tot) - 2.*mu(1)) / dz2
        delta_rho(N_tot) = (mu(1) + mu(N_tot-1) - 2.*mu(N_tot)) / dz2
    end if
     
    rho_new(:) = rho(:) + eps * delta_rho(:)    

end select

end subroutine

subroutine read_input(L, eps, N_steps, fnsh_time, obs_time,model, dz, N_cnf, mean_rho)
real(kind=8), intent(out) :: L, eps, dz, mean_rho
character, intent(out) :: model
integer, intent(out) :: N_steps, fnsh_time, obs_time, N_cnf
open (unit = 53, file = "input.dat", status= "old")

read(53,*) L
read(53,*) N_steps
read(53,*) eps
read(53,*) fnsh_time
read(53,*) obs_time
read(53,*) model
read(53,*) N_cnf
read(53,*) mean_rho

close(53)

if((model.ne."A") .and. (model.ne."B")) then
    print*, "ERROR: integration model should be 'A' or 'B'"
    print*, "Edit input.dat file"
    stop
end if

print*, " Reading done"
dz = L / float(N_steps) ! Bin size

!N_steps = N_steps * 2 + 1! Total # of steps is double of input
!print*,"L,N_steps,eps,fnsh_time,obs_time"
!print*,L,N_steps,eps,fnsh_time,obs_time
end subroutine

subroutine splines(x_points, y_points, n_points, poly_coeff)
!This routine performs a 3rd degree spline of the form
! f_i(x) = poly_coeff(1,i) + poly_coeff(2,i) * (x-x_i) + poly_coeff(3,i) *
! (x-x_i)**2 + poly_coeff(4,i) * (x-x_i)**4
! for x_i < x < x_(i+1) 
implicit none
integer, intent(in) :: n_points
real(kind=8), dimension(n_points), intent(in) :: x_points,y_points
real(kind=8), dimension(4,n_points-1), intent(out) :: poly_coeff
integer :: i_point
real(kind=8) :: dx

dx = x_points(2) - x_points(1)
if( dx .le. 0.000001) dx = 0.000001

!For the first interval x_points(1) < x < x_points(2) take a linear spline
poly_coeff(1,1) = y_points(1)
poly_coeff(2,1) = ( y_points(2) - y_points(1) ) / ( x_points(2) - x_points(1) )
poly_coeff(3,1) = 0.
poly_coeff(4,1) = 0.

!Now define iteratively the polinomial coefficients for each segment
do i_point = 2, n_points - 1 
    dx = x_points(i_point + 1) - x_points(i_point)
    if( dx .le. 0.000001) dx = 0.000001    
    poly_coeff(1,i_point) = y_points(i_point)
    poly_coeff(2,i_point) = poly_coeff(2,i_point-1) + 2*poly_coeff(3,i_point-1)*dx + 3*poly_coeff(4,i_point-1)*dx**2
    poly_coeff(3,i_point) = poly_coeff(3,i_point-1) + 3*poly_coeff(4,i_point-1)*dx
    poly_coeff(4,i_point) = (y_points(i_point+1) - poly_coeff(1,i_point) - poly_coeff(2,i_point)*dx - poly_coeff(3,i_point)*dx**2)/dx**3
end do
end subroutine splines



subroutine spline (x, y, n, poly_coeff)
    !======================================================================
    !  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
    !  for cubic spline interpolation
    !  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
    !  for  x(i) <= x <= x(i+1)
    !  Alex G: January 2010
    !----------------------------------------------------------------------
    !  input..
    !  x = the arrays of data abscissas (in strictly increasing order)
    !  y = the arrays of data ordinates
    !  n = size of the arrays xi() and yi() (n>=2)
    !  output..
    !  b, c, d  = arrays of spline coefficients
    !  comments ...
    !  spline.f90 program is based on fortran version of program spline.f
    !  the accompanying function fspline can be used for interpolation
    !======================================================================
    implicit none
    integer, intent(in) :: n
    real(kind=8), intent(in) :: x(n), y(n)
    real(kind=8) :: b(n), c(n), d(n)
    integer :: i, j, gap
    real(kind=8) ::  h

    !Add by Kevo
    real(kind=8), dimension(4,n-1), intent(out) :: poly_coeff

    gap = n-1
    ! check input
    if ( n < 2 ) return
    if ( n < 3 ) then
    b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
    c(1) = 0.
    d(1) = 0.
b(2) = b(1)
    c(2) = 0.
    d(2) = 0.
    return
    end if
    !
    ! step 1: preparation
    !
    d(1) = x(2) - x(1)
c(2) = (y(2) - y(1))/d(1)
    do i = 2, gap
    d(i) = x(i+1) - x(i)
    b(i) = 2.0*(d(i-1) + d(i))
    c(i+1) = (y(i+1) - y(i))/d(i)
c(i) = c(i+1) - c(i)
    end do
    !
    ! step 2: end conditions 
    !
    b(1) = -d(1)
b(n) = -d(n-1)
    c(1) = 0.0
    c(n) = 0.0
    if(n /= 3) then
    c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
    c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
    c(1) = c(1)*d(1)**2/(x(4)-x(1))
c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
    end if
    !
    ! step 3: forward elimination 
    !
    do i = 2, n
    h = d(i-1)/b(i-1)
    b(i) = b(i) - h*d(i-1)
c(i) = c(i) - h*c(i-1)
    end do
    !
    ! step 4: back substitution
    !
c(n) = c(n)/b(n)
    do j = 1, gap
    i = n-j
c(i) = (c(i) - d(i)*c(i+1))/b(i)
    end do
    !
    ! step 5: compute spline coefficients
    !
b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
    do i = 1, gap
    b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
    d(i) = (c(i+1) - c(i))/d(i)
c(i) = 3.*c(i)
    end do
    c(n) = 3.0*c(n)
d(n) = d(n-1)

!Get all coefficients in one matrix! Add by Kevo
poly_coeff(1,1:n-1) = y(1:n-1)
poly_coeff(2,1:n-1) = b(1:n-1)
poly_coeff(3,1:n-1) = c(1:n-1)
poly_coeff(4,1:n-1) = d(1:n-1)
    end subroutine spline

