module globe
double precision pi, G, mu, m, au, pc, Msun, Mj, Rsun, Rj
double precision Ms, Mp, Mm, Rs, Rp, Rm, density_m, Bs, Bp, Bm, omega_s, omega_p, omega_m
double precision alpha_p, moment_p !alpha_p for planet's moment of inertia
double precision time, ap, wp, angmom_p, am, wm, angmom_m, ap_Roche, am_Roche !orbit
double precision Qs, Qp, Qm, tau_s, tau_p, tau_m !tidal Q and lag time
double precision R_obs, density, v, va, Mach, Cd, Lambda, Lambda1, Lambda2, sigma_m !planet's resistance and moon's conductivity
double precision torque_spt, torque_pst, torque_spm, torque_pmt, torque_pmm, torque_smm !torques for migration timescales and dynamics
end module globe

program main
use globe
implicit none
pi=acos(-1.d0)
G=6.67d-8
mu=4*pi
m=1.67d-24
au=1.5d13
pc=3.086d18
Msun=2.d33
Mj=2.d30
Rsun=7.d10
Rj=7.d9

Ms=Msun
Mp=Mj
Mm=1.d26
Rs=2*Rsun
Rp=2*Rj
Rm=3.d8
density_m=3.d0
alpha_p=0.1
moment_p=alpha_p*Mp*Rp**2
am_Roche=(Mp/Mm)**(1.d0/3.d0)*Rm
omega_s=2.d0*pi/(3.d0*24.d0*3.6d3)
omega_p=2.d0*pi/(40.d0*3.6d3)
tau_s=1.d0
tau_p=1.d0
Bs=1.5d3
Bp=1.d2
Cd=1.0
Lambda2=0.1d9
sigma_m=1.d-8

!call typical_torques
!call magneticQ
!call migration ! plt_migration
!call rTrH ! plt_rTrH.m
call moon ! plt_moon
!call dynamics ! plt_dynamics
!call radio
end program main

subroutine typical_torques
use globe
implicit none
ap=0.05*au
wp=(G*Ms/ap**3)**0.5
am=6.0*Rj
wm=(G*Mp/am**3)**0.5
write(6,*) 'Gamma_sp^t', G*Mp**2*Rs**5/ap**6*(wp-omega_s)*tau_s
write(6,*) 'Gamma_ps^t', G*Ms**2*Rp**5/ap**6*(wp-omega_p)*tau_p
R_obs=Rp*(Bp/Bs)**(1.d0/3.d0)*(Rs/ap)**(-1)
density=5.0*(5.d6/4.6d9)**(-2)*(ap/au)**(-2)
v=(wp-omega_s)*ap
va=Bs*(Rs/ap)**3/sqrt(density*m*mu)
Mach=abs(v)/va
Cd=Mach/sqrt(1.0+Mach**2)
Lambda1=mu*abs(v)/Cd
write(6,*) 'magnetosphere radius (Rp)  ', 'density  ', 'va  ', 'v  ', 'Mach  ', 'Cd'
write(6,'(6E15.6)') R_obs/Rp, density, va, v, Mach, Cd
write(6,*) 'Gamma_sp^m', 4.d0*Bs**2*ap*Rp**2*(Bp/Bs)**(2.d0/3.d0)*(Rs/ap)**4*v/(Lambda1+Lambda2)
v=(wp-omega_p)*ap
va=Bp*(Rp/ap)**3/sqrt(density*m*mu)
Mach=abs(v)/va
Cd=Mach/sqrt(1.0+Mach**2)
Lambda1=mu*abs(v)/Cd
write(6,*) 'Gamma_ps^m', 4.d0*Bp**2*ap*Rs**2*(Bs/Bp)**(2.d0/3.d0)*(Rp/ap)**4*v/(Lambda1+Lambda2)
write(6,*) 'Gamma_pm^t', G*Mm**2*Rp**5/am**6*(wm-omega_p)*tau_p
write(6,*) 'Gamma_mp^t', G*Mp**2*Rm**5/am**6/67.d0
write(6,*) 'Gamma_pm^m', 4.d-9*sigma_m*(wm-omega_p)*am**2*Rm**3*Bp**2*(Rp/am)**6
write(6,*) 'Gamma_sm^m', 4.d-9*sigma_m*(wm-omega_p)*am**2*Rm**3*Bs**2*(Rs/ap)**6
end subroutine typical_torques

subroutine magneticQ
use globe
implicit none
double precision k2, E0, x
Rs=Rsun
Rp=Rj
omega_s=2*pi/(10*24*3.6d3)
Bs=1.d3
Bp=1.d2
k2=0.035
ap=0.05*au
wp=(G*Ms/ap**3)**0.5
E0=G*Mp**2*Rs**5/ap**6
v=(wp-omega_s)*ap
Lambda1=mu*abs(v)/Cd
torque_spm=4.d0*Bs**2*ap*Rp**2*(Bp/Bs)**(2.d0/3.d0)*(Rs/ap)**4*v/(Lambda1+Lambda2)
Qs=E0/abs(torque_spm)/k2
write(6,'(4E15.6)') (wp-omega_s)/omega_s, E0, abs(torque_spm), Qs
open(1,file='Q.dat',form='formatted')
do x=-0.999,1.0,0.0001
 wp=(1.0+x)*omega_s
 ap=(G*Ms/wp**2)**(1.0/3.0)
 E0=G*Mp**2*Rs**5/ap**6
 v=(wp-omega_s)*ap
 Lambda1=mu*abs(v)/Cd
 torque_spm=4.d0*Bs**2*ap*Rp**2*(Bp/Bs)**(2.d0/3.d0)*(Rs/ap)**4*v/(Lambda1+Lambda2)
 Qs=E0/abs(torque_spm)/k2
 write(1,'(7E15.6)') x, ap/au, Lambda1+Lambda2, E0, abs(torque_spm), -log10(Qs), -log(Qs*0.01**(-4.0/3.0)*0.1**(2.0/3.0))
enddo
close(1)
end subroutine magneticQ

subroutine migration
use globe
implicit none
double precision Pp, Pp_min, Pp_max, am_min, am_max
integer i, n
n=10000
omega_s=2.d0*pi/(3.d0*24.d0*3.6d3)
omega_p=2.d0*pi/(40.d0*3.6d3) ! Omega_p=40 hours (fort.3 and fort.4), 3 days (fort.31 and fort.41)
write(6,*) "moon's corotation radius", (G*Mp/omega_p**2)**(1.d0/3.d0)/Rj

Pp_min=1.d0
Pp_max=1.d2
do i=0,n
 Pp=Pp_min+dble(i)*(Pp_max-Pp_min)/dble(n)
 wp=2.d0*pi/(Pp*24.d0*3.6d3)
 ap=(G*Ms/wp**2)**(1.d0/3.d0)
 angmom_p=Mp*(G*Ms*ap)**0.5d0
 torque_spt=G*Mp**2*Rs**5/ap**6*(wp-omega_s)*tau_s
 torque_pst=G*Ms**2*Rp**5/ap**6*(wp-omega_p)*tau_p
 v=(wp-omega_s)*ap
 Lambda1=mu*abs(v)/Cd
 torque_spm=4.d0*Bs**2*ap*Rp**2*(Bp/Bs)**(2.d0/3.d0)*(Rs/ap)**4*v/(Lambda1+Lambda2)
 write(3,'(4E15.6)') Pp, angmom_p/abs(torque_spt)/3.1536d7, angmom_p/abs(torque_pst)/3.1536d7, angmom_p/abs(torque_spm)/3.1536d7
enddo

am_min=6.d0*Rj
ap=1.0*au
am_max=(Mp/Ms)**(1.d0/3.d0)*ap
do i=0,n
 am=am_min+dble(i)*(am_max-am_min)/dble(n)
 wm=(G*Mp/am**3)**0.5d0
 angmom_m=Mm*(G*Mp*am)**0.5d0
 torque_pmt=G*Mm**2*Rp**5/am**6*(wm-omega_p)*tau_p
 torque_pmm=4.d-9*sigma_m*(wm-omega_p)*am**2*Rm**3*Bp**2*(Rp/am)**6
 torque_smm=4.d-9*sigma_m*(wm-omega_p)*am**2*Rm**3*Bs**2*(Rs/ap)**6
 write(4,'(4E15.6)') am/Rj, angmom_m/abs(torque_pmt)/3.1536d7, angmom_m/abs(torque_pmm)/3.1536d7, angmom_m/abs(torque_smm)/3.1536d7
enddo
end subroutine migration

subroutine rTrH 
use globe
implicit none
double precision dotMp, P, rH, rT
integer i,j,ni,nj
dotMp=1.d-7*Mp/3.1536d7
open(1,file='rTrH.dat',form='formatted')
ni=100
nj=100
do i=0, ni
 P=1.d0+(100.d0-1.d0)/dble(ni)*i
 ap=(G*Ms/(2.d0*pi/(P*24.d0*3.6d3))**2)**(1.d0/3.d0)
 rH=(Mp/Ms)**(1.d0/3.d0)*ap
 do j=0, nj
  Bp=50.d0+(500.d0-50.d0)/dble(nj)*j
  rT=(Bp**4*Rp**12/(G*Mp*dotMp**2))**(1.d0/7.d0)
  write(1,'(5E15.3)') P, Bp, rH/Rj, rT/Rj, rH/rT
 enddo
enddo
close(1)
end subroutine rTrH

subroutine moon
use globe
implicit none
double precision Pp, Pp_min, Pp_max, Pm, am_min, am_max
double precision tau1, tau2
integer i, ni, j, nj
ni=200
nj=200
Mm=1.d27 ! 1.d25 for moon1.dat, 1.d26 for moon2.dat, 1.d27 for moon3.dat
Rm=(3.d0/4.d0*Mm/density_m)**(1.d0/3.d0)
omega_p=2.0d0*pi/(3.6d3*40.0d0)
am_min=(G*Mp/omega_p**2)**(1.0d0/3.0d0)
Pp_min=2.0d0
Pp_max=200.0d0
open(1,file='moon3.dat',form='formatted')
do i=0, ni
 Pp=Pp_min+dble(i)*(Pp_max-Pp_min)/dble(ni)
 wp=2.0d0*pi/(3.6d3*24.0d0*Pp)
 if(wp.gt.omega_p) write(6,*) 'spin must be faster than orbit'
 ap=(G*Ms/wp**2)**(1.0d0/3.0d0)
 torque_pst=G*Ms**2*Rp**5/ap**6*(wp-omega_p)*tau_p
 tau1=1.5d0*alpha_p*Mp*Rp**2*omega_p/abs(torque_pst)
 am_max=(Mp/Ms)**(1.0d0/3.0d0)*ap
 do j=1, nj
  am=am_min+dble(j)*(am_max-am_min)/dble(nj)
  wm=(G*Mp/am**3)**0.5d0
  Pm=2.0d0*pi/wm/(3.6d3*24.0d0)
  angmom_m=Mm*(G*Mp*am)**0.5d0
  torque_pmt=G*Mm**2*Rp**5/am**6*(wm-omega_p)*tau_p
  torque_pmm=4.d-9*sigma_m*(wm-omega_p)*am**2*Rm**3*Bp**2*(Rp/am)**6
  torque_smm=4.d-9*sigma_m*(wm-omega_p)*am**2*Rm**3*Bs**2*(Rs/ap)**6
  tau2=angmom_m/(2.d0*abs(torque_pmt+torque_pmm+torque_smm))
  write(1,'(5E15.6)') Pp, Pm, log10(tau1/3.1536d7)/10, log10(tau2/3.1536d7)/10, log10(tau1/tau2)
 enddo
enddo
close(1)
end subroutine moon

subroutine dynamics
use globe
implicit none
double precision h, time0, y1, f1, f2, f3
external:: f1, f2, f3, rk4
double precision cor
integer i, n
n=100000000
h=3.1536d7
time0=3.1536d7*2.d6
time=time0
wp=2.d0*pi/(1.d2*24.d0*3.6d3) !planet's initial orbit, 4d for dyn1 and 4, 10d for dyn2 and 5, 100d for dyn3 and 6
ap=(G*Ms/wp**2)**(1.d0/3.d0)
omega_p=2.d0*pi/(40.d0*3.6d3) !planet's initial spin
cor=(G*Mp/omega_p**2)**(1.d0/3.d0)
am=10.d0*Rj !moon's initial orbit !moon's initial orbit, 6Rj for dyn1-3, 10Rj for dyn4-6
wm=(G*Mp/am**3)**0.5d0
write(6,'(6E15.6)') ap, cor, 6.d0*Rj, am, wm, omega_p
open(1,file='dyn6.dat',form='formatted')

do i=1,n
 angmom_p=Mp*(G*Ms*ap)**0.5d0
 torque_pst=G*Ms**2*Rp**5/ap**6*(wp-omega_p)*tau_p
 density=5.0*(time/(4.6d9*3.1536d7))**(-2)*(ap/au)**(-2)
 v=(wp-omega_s)*ap
 va=Bs*(Rs/ap)**3/sqrt(density*m*mu)
 Mach=abs(v)/va
 Cd=Mach/sqrt(1.0+Mach**2)
 Lambda1=mu*abs(v)/Cd
 torque_spm=4.d0*Bs**2*ap*Rp**2*(Bp*(time/time0)**(-0.267)/Bs)**(2.d0/3.d0)*(Rs/ap)**4*v/(Lambda1+Lambda2)
 angmom_m=Mm*(G*Mp*am)**0.5d0
 torque_pmt=G*Mm**2*Rp**5/am**6*(wm-omega_p)*tau_p
 torque_pmm=4.d-9*sigma_m*(wm-omega_p)*am**2*Rm**3*(Bp*(time/time0)**(-0.267))**2*(Rp/am)**6
 torque_smm=4.d-9*sigma_m*(wm-omega_p)*am**2*Rm**3* Bs                        **2*(Rs/ap)**6

 call rk4(time,ap,h,y1,f1)
 ap=y1
 wp=(G*Ms/ap**3)**0.5d0
 call rk4(time,omega_p,h,y1,f2)
 omega_p=y1
 cor=(G*Mp/omega_p**2)**(1.d0/3.d0)
 call rk4(time,am,h,y1,f3)
 am=y1
 wm=(G*Mp/am**3)**0.5d0
 time=time+h
 if(mod(i,n/10000000).eq.0) write(1,'(11E25.10)') (time-time0)/3.1536d7, ap/au, cor/Rj, am/Rj, omega_p, wm, &
 &                       torque_pst, torque_spm, torque_pmt, torque_pmm, Cd
 if(am.lt.am_Roche) then
  write(6,*) "moon already inside planet's Roche radius"
  stop
 endif
 ap_Roche=(Mp/Ms)**(1.d0/3.d0)*ap
 if(am.gt.ap_Roche) then
  write(6,*) "moon already outside circumplanetary disk"
  stop
 endif
enddo
close(1)
end subroutine dynamics

subroutine rk4(t,y,h,y1,f)
implicit none
double precision t, y, h, y1, f ! y1: y at next time step
external:: f
double precision k1, k2, k3, k4
k1=f(t,y)
k2=f(t+h/2.d0,y+h/2.d0*k1)
k3=f(t+h/2.d0,y+h/2.d0*k2)
k4=f(t+h     ,y+h     *k3)
y1=y+h/6.d0*(k1+2.d0*(k2+k3)+k4)
end subroutine rk4

function f1(t,y)
use globe
implicit none
double precision t, y, f1
f1=2.d0*ap*(-torque_pst-torque_spm)/angmom_p
end function f1

function f2(t,y)
use globe
implicit none
double precision t, y, f2
f2=torque_pst/moment_p
end function f2

function f3(t,y)
use globe
implicit none
double precision t, y, f3
f3=2.d0*am*(-torque_pmt-torque_pmm-torque_smm)/angmom_m
end function f3

subroutine radio
use globe
implicit none
double precision eta, B, W, Theta, d, delta_nu, F, d_limit
Bs=1.d0 !1.5d3 for young and 1.d0 for mature
Bp=1.d1 !1.d2 for young and 1.d1 for mature
eta=2.d-3
Theta=1.58
d=1.d3*pc
delta_nu=100.0*1d6

omega_s=2.d0*pi/(3.d0*24.d0*3.6d3)
ap=5.d-2*au
wp=(G*Ms/ap**3)**(1.d0/2.d0)
v=abs(wp-omega_s)*ap
B=Bs*(Rs/ap)**3
R_obs=Rp*(Bp/Bs)**(1.d0/3.d0)*(Rs/ap)**(-1)
write(6,*) 'R_obs/Rp'
write(6,*) R_obs/Rp
W=eta*(v*B**2/mu)*pi*R_obs**2
write(6,*) 'radio power of star-planet (erg/s)'
write(6,'(E10.3)') W
F=W/(Theta*d**2*delta_nu)/1.d-26
write(6,*) 'radio flux (mJr)'
write(6,'(E10.3)') F
d_limit=(1.d0/F)**(-0.5)*d/pc
write(6,*) 'detection limit (pc)'
write(6,'(E10.3)') d_limit

omega_p=2.d0*pi/(40.d0*3.6d3)
am=6.d0*Rj
wm=(G*Mp/am**3)**(1.d0/2.d0)
v=abs(wm-omega_p)*am
B=Bp
R_obs=Rm
W=eta*(v*B**2/mu)*pi*R_obs**2
write(6,*) 'radio power of planet-moon (erg/s)'
write(6,'(E10.3)') W
F=W/(Theta*d**2*delta_nu)/1.d-26
write(6,*) 'radio flux (mJr)'
write(6,'(E10.3)') F
d_limit=(1.d0/F)**(-0.5)*d/pc
write(6,*) 'detection limit (pc)'
write(6,'(E10.3)') d_limit
end subroutine radio
