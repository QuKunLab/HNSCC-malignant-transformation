# Each tumour sample is processed against its paired normal sample as the reference
pairs=(
  "DZOE2023121932-b1-T1 DZOE2023121932-b1-N1"
  "DZOE2023121932-b1-T2 DZOE2023121932-b1-N2"
  "DZOE2023121932-b1-T3 DZOE2023121932-b1-N3"
  "DZOE2024011405-b1-T1 DZOE2024011405-b1-N1"
)

for pair in "${pairs[@]}"
do
set -- $pair
TUMOR=$1
NORMAL=$2
cnvkit.py batch\
 -p 16\
 /WES/rmdup/${TUMOR}.mappings.sorted.bam\
 -n /WES/rmdup/${NORMAL}.mappings.sorted.bam\
 -t /WES/reference/human_v6_hg38_clean.bed\
 -f /WES/pack/fa/hg38clean.fa\
 --annotate ./cnvkit-master/data/refFlat_hg38.txt\
 --access ./cnvkit-master/data/access-10kb.hg38.bed\
 --output-reference output/${TUMOR}.reference.cnn -d output/\
 --scatter --diagram
done

for sample in DZOE2023121932-b1-T1 DZOE2023121932-b1-T2 DZOE2023121932-b1-T3 DZOE2024011405-b1-T1
do
cnvkit.py segment\
 -o ./output/${sample}.mappings.sorted.cns\
 ./output/${sample}.mappings.sorted.cnr\
 -p 12 -m hmm-tumor
done
