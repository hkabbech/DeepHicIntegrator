# cell=GSE63525_HUVEC
# cell=GSE63525_K562
cell=GSE63525_IMR90

input_path=~/Papantonis_Integrative/data/1_binaries/hic
input_file=$input_path\/$cell\_combined_30.hic

output_path=~/Papantonis_Integrative/data/2_sparse_matrices/sparse_hic/$cell
mkdir $output_path

for i in 1 5 10 20 22
do
	output_file=$output_path\/chr$i\_$cell\.txt
	output_file_bis=$output_path\/chr$i\_$cell\_bis.txt
	java -Xms512m -Xmx2048m -jar ~/jar/juicer_tools.jar dump observed KR $input_file chr$i chr$i BP 25000 $output_file
	# sed "s/^.*NaN.*//g; /^$/d; s/^/chr$i\t/g" $output_file > $output_file_bis
	sed -i "s/^.*NaN.*//g; /^$/d; s/^/chr$i\t/g" $output_file
	echo chr$i done
done
