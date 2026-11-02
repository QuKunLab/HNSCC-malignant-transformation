# using picard to MarkDuplicates
for i in $(ls ./bwa/*bam)
do
bn=$(basename $i)
picard MarkDuplicates\
 -I $i\
 -M rmdup/$bn.matrics\
 -O rmdup/$bn\
 --REMOVE_DUPLICATES true &
done
