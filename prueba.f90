integer i, dimi, iaux, iaux2, p

dimi=3

do i=-3,5

p = abs(i) + 10
iaux = mod(i-1+p*dimi,dimi)+1
iaux2 = abs(mod(floor(float(i-1)/float(dimi)),2))
PBCREFI = iaux+(dimi-2*iaux+1)*iaux2
print*, i, PBCREFI

enddo

end
