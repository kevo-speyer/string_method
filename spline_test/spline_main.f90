program spline_test
implicit none
real(kind=8), dimension(:), allocatable :: x, y
real(kind=8), dimension(:,:), allocatable :: poly_coeff
integer :: n_points = 10

allocate(x(n_points), y(n_points))

!DEBUG
print*, "va 0"

call read_data(x, y, n_points)

!DEBUG
print*, "va1"

allocate(poly_coeff(4,n_points-1))

call splines(x, y, n_points, poly_coeff) !My routine

!DEBUG
print*, "va1.5"

!call spline(x, y, n_points, poly_coeff) ! Internet Routine

!DEBUG
print*, "va2"

call test_splines(poly_coeff, n_points, x, y)

!DEBUG
print*, "run OK"
end program

subroutine test_splines(poly_coeff, n_points, x_data, y_data)
implicit none
integer, intent(in) :: n_points
real(kind=8), dimension(4,n_points-1), intent(in) :: poly_coeff
real(kind=8), dimension(n_points), intent(in) :: x_data, y_data
integer :: i, j, n_eval = 10
real(kind=8) :: x_eval, y_eval, dx

open (unit = 74, file = "eval_spline.dat", status= "unknown")

do i = 1, n_points-1 !lopp data points; i is the interval
    do j = 0, n_eval !inside interval i, evaluate n_eval points of splines
        !x_eval is the point to evaluate the pline
        x_eval = x_data(i) + float(j) / float(n_eval) * (x_data(i+1) - x_data(i))
        dx = x_eval - x_data(i)
        y_eval = poly_coeff(1,i) + poly_coeff(2,i)*dx + poly_coeff(3,i)*dx**2 + poly_coeff(4,i)*dx**3
        write(74,*) x_eval, y_eval
    end do
end do

end subroutine

subroutine read_data(x, y, n_points)
implicit none
real(kind=8), dimension(n_points), intent(out) :: x, y
integer, intent(in) :: n_points 
integer :: i, dummy
open (unit = 53, file = "data_points.dat", status= "old")

!DEBUG
print*, "read_data 1"

!allocate (x(n_points),y(n_points))

!DEBUG
print*, "read_data 2"


do i = 1, n_points

    read(53,*) x(i) , y(i)
    
!DEBUG
print*, "read_data 3",i,x(i) , y(i)


end do

end subroutine

subroutine splines(x_points, y_points, n_points, poly_coeff)
!This routine performs a 3rd degree spline of the form
! f_i(x) = poly_coeff(1,i) + poly_coeff(2,i) * (x-x_i) + poly_coeff(3,i) *
! (x-x_i)**2 + poly_coeff(4,i) * (x-x_i)**3
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



subroutine spline (x, y, n,poly_coeff)
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
