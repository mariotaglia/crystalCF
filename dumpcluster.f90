subroutine dumpclusterinfo
! this subroutine dumps cluster information to disk

use ellipsoid
use transform
use system
use clusters
use MPI, only : rank
implicit none
integer, parameter :: Nlat = 7 ! number of lattices to scan around central cell
integer j,i
integer Nlatx, Nlaty, Nlatz ! periodic images to scan in each direction, may be different due to PBC
real*8 co
integer depth

co = cutoffcluster ! cutoff distance

allocate(listcluster(4,dumpcluster, maxcluster)) ! matrix saving the clusters  
                                                ! first index is type of particle, cellx, celly and cellz
                                                ! second index is the size of the cluster
                                                ! third index is the number of cluster

allocate(tmpcluster(4,dumpcluster))                                                

! check PBC

if(PBC(1).eq.1) then
 Nlatx=Nlat
else 
 Nlatx=0
endif

if(PBC(3).eq.1) then
 Nlaty=Nlat
else 
 Nlaty=0
endif

if(PBC(5).eq.1) then
 Nlatz=Nlat
else 
 Nlatz=0
endif

! loop over cell indexes 


indexcluster = 1

do j = 1, NNN ! loop over all particles in central cell 
   tmpcluster(1,1) = j
   tmpcluster(2,1) = 0
   tmpcluster(3,1) = 0
   tmpcluster(4,1) = 0
   depth = 2
   call findneighbors(j,0,0,0, Nlatx,Nlaty,Nlatz,co, depth)
enddo


if(rank.eq.0) then
print*,'Found ', indexcluster, ' clusters'
print*,'List:'
do j = 1, indexcluster
   print*, (i, listcluster(1,i,indexcluster), i=1,dumpcluster)
enddo
endif

end


recursive subroutine findneighbors(j,jx,jy,jz,Nlatx,Nlaty,Nlatz,co, depth)
use ellipsoid
use transform
use system
use clusters
use MPI, only : rank
implicit none

real*8 co
integer i,ix,iy,iz
integer j,jx,jy,jz
real*8 vect1T(3), vect1R(3),vect2T(3),vect2R(3),vectdiff(3)
integer Nlatx,Nlaty,Nlatz

integer depth
real*8 dist


! Find coordinates of first particle in real space
 vect1T(1) = Rellf(1,j)*delta*dfloat(dimx)+dfloat(dimx*jx)*delta
 vect1T(2) = Rellf(2,j)*delta*dfloat(dimy)+dfloat(dimy*jy)*delta
 vect1T(3) = Rellf(3,j)*delta*dfloat(dimz)+dfloat(dimz*jz)*delta

 vect1R = MATMUL(IMAT,vect1T) ! coordinates of particle in real space
 

! Loop over all other particles      

 do ix = -Nlatx,Nlatx 
 do iy = -Nlaty,Nlaty
 do iz = -Nlatz,Nlatz 
   do i = 1, NNN ! loop over all particles in cell ix,iy,iz

      vect2T(1) = Rellf(1,i)*delta*dfloat(dimx)+dfloat(dimx*ix)*delta
      vect2T(2) = Rellf(2,i)*delta*dfloat(dimy)+dfloat(dimy*iy)*delta
      vect2T(3) = Rellf(3,i)*delta*dfloat(dimz)+dfloat(dimz*iz)*delta
      vect2R = MATMUL(IMAT,vect2T) ! coordinates of particle in real space

      vectdiff = vect1R-vect2R
      dist = norm2(vectdiff) 

      if((dist.ne.0.0).and.(dist.le.co)) then ! found a particle within co distance and it's not the same particle

      tmpcluster(1,depth) = i
      tmpcluster(2,depth) = ix
      tmpcluster(3,depth) = iy
      tmpcluster(4,depth) = iz

      if (depth.eq.dumpcluster) then ! we are done

      listcluster(:,:,indexcluster)=tmpcluster(:,:)
      if(indexcluster.eq.maxcluster) then
              if(rank.eq.0)print*, 'maximum number of clusters reached, increase maxcluster'
             call endall
      endif  
      indexcluster = indexcluster + 1

      else

       depth = depth + 1 ! increase depth
       call findneighbors(i,ix,iy,iz, Nlatx,Nlaty,Nlatz,co, depth) ! find neighbors for this particle
      endif ! depth
        
      endif ! dist
 
    enddo ! i
  enddo ! ix
  enddo ! iy
  enddo ! iz
end





