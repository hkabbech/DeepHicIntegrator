echo "chromosome: $1"
echo "file: $2"
juicetools dump observed KR $2 $1 $1 BP 25000 $1.txt
sed -i 's/^.*NaN.*\n//g;s/^.*NaN\n//g;s/^/chr20\t/g' $1.txt