program ver_vel
implicit none
real,dimension(20,131,105)::u,v,du_dx,dv_dx,du_dy,dv_dy,div,omega,avg_div,cor_div,cor_avg_div,cor_omega,w,w_cor
real::t,f(131),dx(131),bias(131,105)
integer::i,j,k,p(20),dy

open(300,file='pres_levels.dat')
open(100,file='u_20_levels.dat')
open(200,file='v_20_levels.dat')
open(400,file='corrected_divergence.dat')
open(500,file='omega.dat')
open(600,file='corrected_omega.dat')
open(1000,file='vertical_velocity.dat')
open(700,file='vertical_velocity_925.dat')
open(900,file='vertical_velocity_850.dat')
open(950,file='vertical_velocity_600.dat')
open(970,file='vertical_velocity_208.dat')
open(800,file='divergence.dat')

do k=1,20
read(300,*)p(k)
end do

do k=1,20
do i=1,131
read(100,*)(u(k,i,j),j=1,105)
read(200,*)(v(k,i,j),j=1,105)
end do																							 
end do

dy=0.5*110000.0
i=1
do t=-20.0,45.0,0.5
dx(i)=0.625*cos(t*(3.14/180.0))*110000.0
if(i.le.131)i=i+1
end do
i=1
do t=-20.0,45.0,0.5
f(i)=2.0*7.292*(10.0**(-5.0))*sin(t*(3.14/180.0))
if(i.le.131)i=i+1
end do

call derivative(u,du_dx,du_dy)
call derivative(v,dv_dx,dv_dy)

omega(1,i,j)=0.0
!avg_div(1,i,j)=0.0
cor_omega(1,i,j)=0.0
!cor_avg_div(1,i,j)=0.0

do k=1,20
do i=1,131
do j=1,105
div(k,i,j)=du_dx(k,i,j)+dv_dy(k,i,j)
if(k.le.19)then
avg_div(k+1,i,j)=(div(k+1,i,j)+div(k,i,j))/2.0
omega(k+1,i,j)=omega(k,i,j)+(((p(k)-p(k+1))*100)*avg_div(k+1,i,j))
w(k,i,j)=(-omega(k,i,j))/(1.225*(9.81))
end if
end do 
end do
end do

do i=1,131
do j=1,105
bias(i,j)=omega(20,i,j)/((985-108)*100)
end do 
end do

do k=1,20
do i=1,131
do j=1,105
cor_div(k,i,j)=	div(k,i,j)-bias(i,j)
if(k.le.19)then
cor_avg_div(k+1,i,j)=(cor_div(k+1,i,j)+cor_div(k,i,j))/2.0
cor_omega(k+1,i,j)=	cor_omega(k,i,j)+(((p(k)-p(k+1))*100)*cor_avg_div(k+1,i,j))
w_cor(k,i,j)=(-cor_omega(k,i,j))/(1.225*(9.81))
end if
end do 
end do
end do


do k=1,20
do i=1,131
write(800,10)(div(k,i,j),j=1,105)
write(400,10)(cor_div(k,i,j),j=1,105)
write(500,10)(omega(k,i,j),j=1,105)
write(600,10)(cor_omega(k,i,j),j=1,105)
write(1000,10)(w_cor(k,i,j),j=1,105)
write(700,10)(w_cor(3,i,j),j=1,105)
write(900,10)(w_cor(5,i,j),j=1,105)
write(950,10)(w_cor(9,i,j),j=1,105)
write(970,10)(w_cor(16,i,j),j=1,105)
end do																					  
end do

10 format(105(f12.7,1x))

contains
subroutine derivative(p,dp_dx,dp_dy)
implicit none
integer::k,i,j
real,dimension(20,131,105)::p,dp_dx,dp_dy

do k=1,20
do i=1,131
do j=1,105

if(j.gt.1.and.j.lt.105)dp_dx(k,i,j)=(p(k,i,j+1)-p(k,i,j-1))/(2*dx(i))
if(i.gt.1.and.i.lt.131)dp_dy(k,i,j)=(p(k,i+1,j)-p(k,i-1,j))/(2*dy)

if(i.eq.1)dp_dy(k,i,j)=(p(k,i+1,j)-p(k,i,j))/(dy)
if(j.eq.1)dp_dx(k,i,j)=(p(k,i,j+1)-p(k,i,j))/(dx(i))

if(i.eq.131)dp_dy(k,i,j)=(p(k,i,j)-p(k,i-1,j))/(dy)
if(j.eq.105)dp_dx(k,i,j)=(p(k,i,j)-p(k,i,j-1))/(dx(i))
 
end do 
end do
end do

end subroutine derivative
end program ver_vel

