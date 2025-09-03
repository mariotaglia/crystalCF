import os

os.mkdir("up")
os.chdir("up")
os.system("cp ../DEFINITIONS.txt .")
os.system("~/develop/crystalCF/crystalCF")

os.system("cat F_tot_gcanon.dat | awk '{print $2}' > FF")
os.system('for i in system.*.dat ; do cat $i | grep "phisolv =" ; done > P')
os.system('for i in system.*.dat ; do cat $i | grep "musolv" ; done > M')
os.system("cat P | awk '{print $3}' > PP")
os.system("cat M | awk '{print $3}' > MM")
os.system('paste MM PP FF > MPF')

os.chdir("..")
os.mkdir("down")
os.chdir("down")
os.system("cp ../DEFINITIONS.txt .")
os.system("cp ../up/out.out in.in")


with open('DEFINITIONS.txt', 'r') as file:
        file_content = file.read()

modified_file_content = file_content.replace('nkp -20 0.01 0.0', 'nkp 0.0 0.01 -20.0')
modified_file_content = modified_file_content.replace('infile 0', 'infile 2')

with open('DEFINITIONS.txt', 'w') as file:
        file.write(modified_file_content)

os.system("~/develop/crystalCF/crystalCF")
os.system("cat F_tot_gcanon.dat | awk '{print $2}' > FF")
os.system('for i in system.*.dat ; do cat $i | grep "phisolv =" ; done > P')
os.system('for i in system.*.dat ; do cat $i | grep "musolv" ; done > M')
os.system("cat P | awk '{print $3}' > PP")
os.system("cat M | awk '{print $3}' > MM")
os.system('paste MM PP FF > MPF')

os.chdir("..")

