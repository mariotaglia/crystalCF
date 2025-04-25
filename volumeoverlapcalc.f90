subroutine calcvolumeoverlap
! this subroutine calculate the solaping volume for a spherical NP superlattice.

use ellipsoid
use transform
use system
use chainsdat
use MPI, only : rank
implicit none
real*8 vect1T(3), vect1R(3),vect2T(3),vect2R(3)
integer j, k, ir,it,ip, ix, iy, iz, nn
real*8 x,y,z, vect2
real*8, external :: pos2x, pos2y, pos2z
real*8, parameter :: pi = acos(-1.0)
integer Nradius, Ntita, Nphi
real*8 dradius, dtita, dphi
real*8 sum, rlig
real*8 rA, rB, tita, phi
real*8 r2, radius
integer flag

nn = 2
Ntita = 100
dtita = 2.0*pi/float(Ntita)
Nphi=100
dphi = pi/float(Nphi)
Nradius = 200
rlig = 0.088*(long+1)+0.525
dradius = rlig/float(Nradius)
sum = 0.0

do j = 1, NNN
	do k = 1, NNN
		if (k.ge.j) then
		rA = Aell(1,j)
		rB = Aell(1,k)

		vect1T(1) = pos2x(j,0) 
		vect1T(2) = pos2y(j,0)
		vect1T(3) = pos2z(j,0)
		vect1R = MATMUL(IMAT,vect1T) ! coordinates of particle in real space

		do ir = 1, Nradius
		radius = float(ir)*dradius + rA
		do it = 1, Ntita
		tita = dtita*float(it-1)
		do ip = 1, Nphi
			phi = dphi*float(ip-1)

			x = radius*sin(phi)*cos(tita) + vect1R(1)
			y = radius*sin(phi)*sin(tita) + vect1R(2)
			z = radius*cos(phi) + vect1R(3)

			flag = 0

			do ix = -nn, nn
			do iy = -nn, nn
			do iz = -nn, nn

			vect2T(1) = pos2x(k,ix) 
			vect2T(2) = pos2y(k,iy)
			vect2T(3) = pos2z(k,iz)
			vect2R = MATMUL(IMAT,vect2T) ! coordinates of particle in real space
			
			if(j.eq.k) then
			if((ix.ne.0).or.(iy.ne.0).or.(iz.ne.0)) then
				vect2 = (vect2R(1) - x)**2 + (vect2R(2) - y)**2 + (vect2R(3) - z)**2
				r2 = (rlig+rA)**2
				if(vect2.lt.r2) then
	  				flag = 1
	  				exit
	  			endif
  			endif
  			else
	  			vect2 = (vect2R(1) - x)**2 + (vect2R(2) - y)**2 + (vect2R(3) - z)**2
				r2 = (rA+rlig)**2
				if(vect2.lt.r2) then
	  				flag = 1
	  				exit
  				endif
  			endif

			enddo
			if(flag.eq.1)exit
			enddo
			if(flag.eq.1)exit
			enddo

			if(flag.eq.0) sum = sum + radius*radius*sin(phi)*dtita*dphi*dradius
		enddo
		enddo
		enddo
	end if
	enddo
enddo

if(rank.eq.0) then
open(unit=46, file= 'volume_overlap.dat', access='APPEND')
write(46,*) sum
close(46)
endif

end subroutine

double precision function pos2x(j,jx)
use ellipsoid
use transform
use system
implicit none
integer j, jx
pos2x =  Rellf(1,j)*delta*dfloat(dimx)+dfloat(dimx*jx)*delta
end

double precision function pos2y(j,jy)
use ellipsoid
use transform
use system
implicit none
integer j, jy
pos2y =  Rellf(2,j)*delta*dfloat(dimy)+dfloat(dimy*jy)*delta
end

double precision function pos2z(j,jz)
use ellipsoid
use transform
use system
implicit none
integer j, jz
pos2z =  Rellf(3,j)*delta*dfloat(dimz)+dfloat(dimz*jz)*delta
end