rm tmp*

for i in L* ; do
       	echo $i > pair_$i.dat 
	cd $i 
	for j in dim* ; do

        ll=`tail -1 $j/F_tot_gcanon.dat`
	echo $ll

        if [ "$ll" != "" ] ; then

            echo $i $j `tail -1 $j/F_tot_gcanon.dat` | sed 's/dim_//' | sed 's/L//' >> ../tmp_$i.dat 
	
        fi
	
         done 
	 cd .. 
	 sort -k2n tmp_$i.dat >> pair_$i.dat
 done

rm tmp*
