listN='1.74
1.48
1.38
1.31
1.24
1.12
1.06
1
0.88
0.77
0.72
0.67
0.56'

list='5.5
5.75
6
6.25
6.5
6.75
7
7.25
7.5
7.75
8
8.25
8.5
8.75
9
9.25
9.5
9.75
10
10.25
10.5
23.5
23.75
24'


for j in $listN ; do
      mkdir L$j
      cd L$j

for i in $list ; do
      mkdir dim_$i
      cd dim_$i
      rm *.dat
      cp ../../tosubmit.sh .
      cp ../../DEFINITIONS.txt .
      sed -i "s/_NAME/dim$i/g" tosubmit.sh
      sed -i "s/_DIM/$i/g" DEFINITIONS.txt
      sed -i "s/_RS/$j/g" DEFINITIONS.txt
      sbatch tosubmit.sh
      cd ..
done
      cd ..
done

