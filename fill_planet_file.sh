filename=planet0.dat

dum=$(wc -l $filename)

nlin=$(echo $dum | cut -d'p' -f 1)

n=0


while read line; do
#reading each line
#echo "$n"
#echo "$line"
n=$((n+1))
if [[ $n -eq $nlin ]]
then
#echo "$line"
string="$(echo "$line")"
fi

done < $filename


dum=$((nlin-1))
dum2=$dum

for i in {1..5}
do
dum2=$((dum2+1))
echo "${string/$dum/$dum2}" 
done

