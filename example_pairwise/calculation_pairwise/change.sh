sed -i 's/infile 0/infile 2/g' DEFINITIONS.txt
sed -i 's/nst 4/nst 1/g' DEFINITIONS.txt
sed -i '52,54d' DEFINITIONS.txt
mv out.004.dat in.in
rm *.dat
