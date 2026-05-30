cnvkit.py batch\
 -p 16\
 /WES/rmdup/*T*bam\
 -n /WES/rmdup/*N*bam\
 -t /WES/reference/human_v6_hg38_clean.bed\
 -f /WES/pack/fa/hg38clean.fa\
 --annotate ./cnvkit-master/data/refFlat_hg38.txt\
 --access ./cnvkit-master/data/access-10kb.hg38.bed\
 --output-reference output/reference.cnn -d output/\
 --scatter --diagram

for sample in DZOE2023121932-b1-T1 DZOE2023121932-b1-T2 DZOE2023121932-b1-T3 DZOE2024011405-b1-T1
do
cnvkit.py segment\
 -o ./output/${sample}.mappings.sorted.cns\
 ./output/${sample}.mappings.sorted.cnr\
 -p 12 -m hmm-tumor
done

