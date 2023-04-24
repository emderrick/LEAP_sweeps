#!/usr/bin/bash
for file in */*.gbk
do
out="${file//prokkaoutput*/_EPSPS.txt}"
name="${out//_EPSPS.txt/_EPSPS}"
nameout="${name//_EPSPS/name.txt}"
out2="${nameout//name.txt/_FINAL.txt}"
out3="${nameout//name.txt/_nowhite.txt}"
out4="${out3//nowhite.txt/EPSPS.fa}"
grep -A14 "3-phosphoshikimate 1-carboxyvinyltransferase" $file > $out
echo $name > $nameout
sed -i 's/^/>/' $nameout
cat $nameout $out > $out2
cat $out2 | tr -d "[:blank:]" > $out3
sed 's/^.*translation="//' $out3 > $out4
sed -i '/t/d' $out4
sed -i '/[/]/d' $out4
sed -i '/CDS/d' $out4
sed -i 's/"//g' $out4
sed -i '/ORIGIN/d' $out4
sed -i '/^$/d' $out4
done

rm *EPSPS.txt
rm *nowhite.txt
rm *name.txt
rm *FINAL.txt
cat *EPSPS.fa > ALL_EPSPS.fa

