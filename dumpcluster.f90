subroutine dumpclusterinfo
! this subroutine dumps cluster information to disk

use ellipsoid
use transform
use system
use clusters
use MPI, only : rank
implicit none
integer, parameter :: Nlat = 7 ! number of lattices to scan around central cell
integer j, i, jx, jy, jz, ix, iy ,iz, ii, jj
integer Nlatx, Nlaty, Nlatz ! periodic images to scan in each direction, may be different due to PBC
real*8 co
integer depth
integer combs
real*8, allocatable :: distlisttmp(:,:)
real*8, allocatable :: distlist(:,:,:)
integer d1, d2, dd
real*8 vect1T(3), vect1R(3),vect2T(3),vect2R(3),vectdiff(3)
real*8, external :: posx, posy, posz
integer, allocatable :: listcluster_in(:,:,:)
real*8, allocatable :: distlist_in(:,:,:)
real*8 weight(maxcluster)
real*8 tol
integer found, flag
integer indexcluster_in
character*15 filename


if(rank.eq.0) then
print*,"Enter dump-cluster"
print*, "Search level:", dumpcluster
if(cluster_same.eq.0)print*, 'Treat all particles as different'
if(cluster_same.eq.1)print*, 'Treat all particles as equal'
endif

co = cutoffcluster ! cutoff distance


allocate(listcluster(4,dumpcluster, maxcluster)) ! matrix saving the clusters  
                                                ! first index is type of particle, cellx, celly and cellz
                                                ! second index is the size of the cluster
                                                ! third index is the number of cluster

allocate(listcluster_in(4,dumpcluster, maxcluster))
allocate(tmpcluster(4,dumpcluster))                                                

combs = (dumpcluster-1)*dumpcluster/2 ! possible ways to arrange a cluster of dumpcluster particles
                                      ! used to make an ordered list of distances and then find unique clusters
allocate(distlisttmp(3,combs))  
allocate(distlist(3,combs, maxcluster))  ! distlist: index1: dist, i, j
allocate(distlist_in(3,combs, maxcluster))

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

indexcluster = 0 

do j = 1, NNN ! loop over all particles in central cell 
   tmpcluster(1,1) = j
   tmpcluster(2,1) = 0
   tmpcluster(3,1) = 0
   tmpcluster(4,1) = 0
   depth = 2
   call findneighbors(j,0,0,0, Nlatx,Nlaty,Nlatz,co, depth)
enddo



! Now, find only non-equivalent clusters
! To do that, we will make am ordered list of the N(N-1)/2 distances in the cluster

do ii = 1, indexcluster

dd = 0 

do d1 = 1, dumpcluster
 do d2 = d1+1, dumpcluster 

 dd = dd + 1

 j = listcluster(1,d1,ii)
 jx = listcluster(2,d1,ii)
 jy = listcluster(3,d1,ii)
 jz = listcluster(4,d1,ii)

 i = listcluster(1,d2,ii)
 ix = listcluster(2,d2,ii)
 iy = listcluster(3,d2,ii)
 iz = listcluster(4,d2,ii)

 vect1T(1) = posx(j,jx) 
 vect1T(2) = posy(j,jy)
 vect1T(3) = posz(j,jz)
 vect1R = MATMUL(IMAT,vect1T) ! coordinates of particle in real space
       
 vect2T(1) = posx(i,ix)
 vect2T(2) = posy(i,iy)
 vect2T(3) = posz(i,iz)
 vect2R = MATMUL(IMAT,vect2T) ! coordinates of particle in real space

 vectdiff = vect1R-vect2R
   
 distlisttmp(1, dd) = norm2(vectdiff) 

 if(i.lt.j) then
   distlisttmp(2, dd) = i
   distlisttmp(3, dd) = j
 else
   distlisttmp(2, dd) = j
   distlisttmp(3, dd) = i
 endif


 enddo ! d2
enddo ! d1

call bubble_sort(distlisttmp, combs)

distlist(:,:,ii) = distlisttmp(:,:)



enddo ! ii
! Now find unique clusters

indexcluster_in = 0 ! index of unique clusters
tol = 1e-3 ! tolerance for comparing distances
weight = 0.0

do ii = 1, indexcluster ! loop over all clusters

  found = 0 ! cluster ii not in the list

  do jj = 1, indexcluster_in ! loop over all cluster already in list 

    do dd = 1, combs ! loop over all distance pairs
    flag = 1 ! assume it is there

    if(cluster_same.eq.0) then ! cluser_same = 1 -> treat all particles as equal for cluster
         if((abs(distlist(1,dd,ii)-distlist_in(1,dd,jj)).gt.tol).or.  &
            (abs(distlist(2,dd,ii)-distlist_in(2,dd,jj)).gt.tol).or.  &
            (abs(distlist(3,dd,ii)-distlist_in(3,dd,jj)).gt.tol)) then
          flag = 0 ! not there
          exit
         endif
    else
         if((abs(distlist(1,dd,ii)-distlist_in(1,dd,jj)).gt.tol)) then
          flag = 0 ! not there
          exit
         endif
    endif
  
    enddo ! dd   

    if(flag.eq.1) then
            found = jj ! already in list
!            print*,'found', distlist(1,dd,ii), distlist_in(1,dd,jj), ii, jj
            exit
    endif
  enddo ! jj   
  
  if(found.ne.0) weight(found) = weight(found) + 1 ! increase weight 

  if(found.eq.0) then ! not in list, add to it
       indexcluster_in = indexcluster_in + 1
       distlist_in(:,:,indexcluster_in) = distlist(:,:,ii)
       listcluster_in(:,:,indexcluster_in) = listcluster(:,:,ii)
       weight(indexcluster_in) = 1.0
   endif

enddo ! ii

!!! PRINT ALL

if(rank.eq.0) then
print*,'Found ', indexcluster, ' clusters'
print*,'List:'
print*, 'Distance'
   print*, 'Distances between all particles'
   if(cluster_same.eq.0)print*, 'Followed by particle types involucrated:'
do j = 1, indexcluster
!   print*, listcluster(:,:,j)

   print*, distlist(1,:,j)
   if(cluster_same.eq.0)print*, int(distlist(2,:,j))
   if(cluster_same.eq.0)print*, int(distlist(3,:,j))
   if(cluster_same.eq.0)print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
enddo
endif

   if(cluster_same.eq.0)print*,""
   if(cluster_same.eq.0)print*,""
   if(cluster_same.eq.0)print*,""
   if(cluster_same.eq.0)print*,""

if(rank.eq.0) then
print*,'Found ', indexcluster_in, ' unique clusters'
   if(cluster_same.eq.0)print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'   
print*,'List of distances'
print*,'Number of occurences'
if(cluster_same.eq.0)print*, 'Particle types involucrated'
   if(cluster_same.eq.0)print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'   

    open(unit=45, file='distances.dat')
    write(45,*)indexcluster_in

do j = 1, indexcluster_in
!   print*, listcluster_in(:,:,j)
   print*, distlist_in(1,:,j)
   print*, int(weight(j))
   write(45,*)distlist_in(1,:,j)
   write(45,*)int(weight(j))


   if(cluster_same.eq.0)print*, int(distlist_in(2,:,j))
   if(cluster_same.eq.0)write(45,*)int(distlist_in(2,:,j))

   if(cluster_same.eq.0)print*, int(distlist_in(3,:,j))
   if(cluster_same.eq.0)write(45,*)int(distlist_in(3,:,j))

   if(cluster_same.eq.0)print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'   
enddo
endif

print*,'Saved to distances.dat'
if(rank.eq.0)print*, 'Saving coordinates of unique clusters to disk, files cluster.***.xyz'

do j = 1, indexcluster_in

    write(filename,'(A8, I3.3, A4)')'cluster.', j, '.xyz'
    open(unit=45, file=filename)

    do i = 1, dumpcluster 
    write(45,*) posx(listcluster_in(1,i,j), listcluster_in(2,i,j)), &
                posy(listcluster_in(1,i,j), listcluster_in(3,i,j)), &
                posz(listcluster_in(1,i,j), listcluster_in(4,i,j))
    enddo
    close(45)
enddo

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
real*8, external :: posx, posy, posz
integer depth
real*8 dist
integer flagin, dd

! Find coordinates of first particle in real space
 vect1T(1) = posx(j,jx) 
 vect1T(2) = posy(j,jy)
 vect1T(3) = posz(j,jz)

 vect1R = MATMUL(IMAT,vect1T) ! coordinates of particle in real space
 

! Loop over all other particles      

 do ix = -Nlatx,Nlatx 
 do iy = -Nlaty,Nlaty
 do iz = -Nlatz,Nlatz 
   do i = 1, NNN ! loop over all particles in cell ix,iy,iz

      vect2T(1) = posx(i,ix)
      vect2T(2) = posy(i,iy)
      vect2T(3) = posz(i,iz)
      vect2R = MATMUL(IMAT,vect2T) ! coordinates of particle in real space

      vectdiff = vect1R-vect2R
      dist = norm2(vectdiff) 

      if((dist.ne.0.0).and.(dist.le.co)) then ! found a particle within co distance and it's not the same particle

!      print*, 'found 1', ix,iy,iz, depth, dist

! Before adding it to the list, we must check that the new particle is not in the cluster already

      flagin = 0 ! not in the cluster

      do dd = 1, depth-1
        if ((tmpcluster(1,dd).eq.i).and. & 
        (tmpcluster(2,dd).eq.ix).and. &
        (tmpcluster(3,dd).eq.iy).and. &
        (tmpcluster(4,dd).eq.iz)) then
         
        flagin = 1 ! already in the cluster 
        exit
        endif
      enddo

     if (flagin.eq.0) then

      tmpcluster(1,depth) = i
      tmpcluster(2,depth) = ix
      tmpcluster(3,depth) = iy
      tmpcluster(4,depth) = iz

      if (depth.eq.dumpcluster) then ! we are done

      indexcluster = indexcluster + 1
      listcluster(:,:,indexcluster)=tmpcluster(:,:)



      if(indexcluster.eq.maxcluster) then
              if(rank.eq.0)print*, 'maximum number of clusters reached, increase maxcluster'
             call endall
      endif  

      else

       call findneighbors(i,ix,iy,iz, Nlatx,Nlaty,Nlatz,co, depth+1) ! find neighbors for this particle
      endif ! depth

      endif ! flagin      
      endif ! dist
 
    enddo ! i
  enddo ! ix
  enddo ! iy
  enddo ! iz
end


double precision function posx(j,jx)
use ellipsoid
use transform
use system
implicit none
integer j, jx
posx =  Rellf(1,j)*delta*dfloat(dimx)+dfloat(dimx*jx)*delta
end

double precision function posy(j,jy)
use ellipsoid
use transform
use system
implicit none
integer j, jy
posy =  Rellf(2,j)*delta*dfloat(dimy)+dfloat(dimy*jy)*delta
end

double precision function posz(j,jz)
use ellipsoid
use transform
use system
implicit none
integer j, jz
posz =  Rellf(3,j)*delta*dfloat(dimz)+dfloat(dimz*jz)*delta
end


! ChatGTP routine to order an array
subroutine bubble_sort(arr, ss)
integer ss
real*8 arr(3,ss)
integer i, j
real*8 temp(3)

        do i = 1, ss-1
            do j = 1, ss-i
                if (arr(1,j) > arr(1,j+1)) then
                    ! Swap elements
                    temp(:) = arr(:,j)
                    arr(:,j) = arr(:,j+1)
                    arr(:,j+1) = temp(:)
                end if
            end do
        end do
    
end subroutine bubble_sort

