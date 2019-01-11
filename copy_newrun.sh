f=${1:0:-4}_$2.dat
cp $1 $f
d=$(expr $4 - $3)
sed -i "2s/.*/$2/" $f
sed -i "3s/.*/$3 $4 $d/" $f

