#!/usr/bin/bash


Read_pairs_num=100000  # the number of read-pairs for each sample

echo "copying the raw experimental RNAseq folder to the current directory..."
cp -r /localdisk/data/BPSM/ICA1/fastq . 

##################### Code block 1 : Integrity check ###################################
echo "Quick check for integrity is processing..."
rm -f integrity_check.txt  # in case the output file was already there
touch integrity_check.txt
# check every raw seq file in the folder fastq/, if they all have exactly 100k*4=400k lines
ls fastq/ | grep fq.gz | while read file;do
  line=$(zcat fastq/${file} | wc -l)
  # if a file does not have 400k lines, print it out to the integrity_check.txt file.
  if [ ${line} != $((${Read_pairs_num}*4)) ];then
    echo "Seq ${file} is not intact, having ${line} lines." >> integrity_check.txt
  fi
done 
defeat=$(cat integrity_check.txt | wc -l) # get the lines(number of seq files not having 400k lines).
# if all the seqs have 400k lines, here the txt should be empty (0 rows)
if [ ${defeat} == 0 ];then
  echo "quick check complete, All seq files are intact." >> integrity_check.txt
else
  echo "quick check complete, ${defeat} RNAseq files are not intact." >> integrity_check.txt
fi

##################### Code block 2 : Quality check ###################################
mkdir report # fastqc does not create the folder automatically
fastqc fastq/*fq.gz -o report --extract  # perform fastqc and output to report folder, unzipping all the files.
rm -f summary_per_item.tmp summary_per_file.tmp
cat report/*fastqc/summary.txt | cut -f2 | sort | uniq | while read item;do
fail=$(grep "${item}" report/*fastqc/summary.txt | grep -c FAIL)
warn=$(grep "${item}" report/*fastqc/summary.txt | grep -c WARN)
pass=$(grep "${item}" report/*fastqc/summary.txt | grep -c PASS)
echo -e "${item}\t${fail}\t${warn}\t${pass}" >> summary_per_item.tmp
done
echo -e "Items\tFAIL\tWARN\tPASS" > summary_per_item.txt
sort -t$'\t' -k2,2r -k3,3r summary_per_item.tmp >> summary_per_item.txt
#
ls report | grep fastqc$ | while read file;do
fail=$(grep -c FAIL report/${file}/summary.txt)
warn=$(grep -c WARN report/${file}/summary.txt)
pass=$(grep -c PASS report/${file}/summary.txt)
echo -e "${file}\t${fail}\t${warn}\t${pass}" >> summary_per_file.tmp
done
echo -e "Seq_file\tFAIL\tWARN\tPASS" > summary_per_file.txt
sort -t$'\t' -k2,2r -k3,3r summary_per_file.tmp >> summary_per_file.txt

rm -f summary_per_item.tmp summary_per_file.tmp


cp -r /localdisk/data/BPSM/ICA1/Tcongo_genome/ .
########## Code block 3 : Alignment ########
bowtie2-build Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz  Tcongo_genome/tc
mkdir Align_bam
ls fastq | grep 1.fq.gz > fq1.tmp
ls fastq | grep 2.fq.gz > fq2.tmp
paste fq1.tmp fq2.tmp | while read fq1 fq2;do
bowtie2 -p 50 -x Tcongo_genome/tc -1 fastq/${fq1} -2 fastq/${fq2} | samtools view -bS | samtools sort > Align_bam/${fq1:0:8}.bam
done
rm -f fq1.tmp fq2.tmp

########## Code block 4 : Counts data ########
cp /localdisk/data/BPSM/ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed .
mkdir counts
ls Align_bam | while read bamfile;do
bedtools coverage -a TriTrypDB-46_TcongolenseIL3000_2019.bed -b Align_bam/${bamfile} > counts/${bamfile:0:8}cov.out
done


########## Code block 5 : Mean for groups ########
mkdir groups
rm -f groups/*.tmp
tail -n +2 fastq/Tco.fqfiles | sort -k2,2 -k4,4 -k5,5 -k3,3 | while read sname stype rep time treat end1 end2;do
if [ ${rep} == 1 ];then
filename=${stype}_${time}h_${treat}
touch groups/${filename}_file.tmp
fi
echo "${end1:0:8}" >> groups/${filename}_file.tmp
done
cut -f1,4 TriTrypDB-46_TcongolenseIL3000_2019.bed > gene.tmp
cut -f5 TriTrypDB-46_TcongolenseIL3000_2019.bed > description.tmp

replicate=0
rm -f groups/rep*
ls groups | grep file.tmp | while read groupfile;do  # groupfile:Clone1_0h_Uninduced_file.tmp 
cat groups/${groupfile} | while read seq_num;do  # seq_num: Tco-5053
((replicate+=1))
cut -f6 counts/${seq_num}cov.out > groups/rep${replicate}.tmp
done
paste groups/rep* | awk 'BEGIN{sum=0}{for(i=1;i<=NF;i++){sum+=$i}{print sum/NF; sum=0}}' > groups/${groupfile%file*}mean.tmp
paste gene.tmp groups/${groupfile%file*}mean.tmp description.tmp > groups/${groupfile%file*}mean.txt
rm -f groups/rep*
done
sed -i '1i Position\tGene\tMean_exp\tDescription'  groups/*mean.txt
rm -f groups/*file.tmp

########## Code block 6 : Fold change comparisons ########
## 6.1 Specific time point. Ref: WT uninduced
mkdir fixed_time
tail -n +2 fastq/Tco.fqfiles | cut -f4 | sort | uniq | while read time;do   # 0,24,48
ref="WT_${time}h_Uninduced_mean.tmp"
ls groups | grep WT_${time}h_.*.tmp | sort -t"_" -k3,3r > fixed_time/order.tmp
ls groups | grep Clone.*_${time}h_.*.tmp | sort -t"_" -k1,1 -k3,3r >> fixed_time/order.tmp

paste gene.tmp > fixed_time/${time}h.txt
heading="Position\tGene"
while read file;do  # WT_24h_Uninduced_mean.tmp
paste groups/${ref} groups/${file} | awk 'BEGIN{FS="\t"}{if($1>0){print $2/$1}else{if("'${file}'"=="'${ref}'"){print 1}else{print "Inf"}}}' > fixed_time/${file%mean*}f.tmp
paste fixed_time/${time}h.txt fixed_time/${file%mean*}f.tmp > tmp && mv tmp fixed_time/${time}h.txt
heading="${heading}\t${file%_mean*}"
done < fixed_time/order.tmp
paste fixed_time/${time}h.txt description.tmp | sort -k4,4r > tmp && mv tmp fixed_time/${time}h.txt
sed -i '1i '${heading}'\tDescription'  fixed_time/${time}h.txt

done


## 6.2 Specific types. Ref: 0h uninduced
mkdir fixed_type
tail -n +2 fastq/Tco.fqfiles | cut -f2 | sort | uniq | while read types;do   # WT, Clone1,Clone2,Clonen
ref="${types}_0h_Uninduced_mean.tmp"
ls groups | grep ${types}.*.tmp | sort -t"_" -k2,2 -k3,3r > fixed_type/order.tmp

paste gene.tmp > fixed_type/${types}.txt
heading="Position\tGene"
while read file;do  # WT_24h_Uninduced_mean.tmp
paste groups/${ref} groups/${file} | awk 'BEGIN{FS="\t"}{if($1>0){print $2/$1}else{if("'${file}'"=="'${ref}'"){print 1}else{print "Inf"}}}' > fixed_type/${file%mean*}f.tmp
paste fixed_type/${types}.txt fixed_type/${file%mean*}f.tmp > tmp && mv tmp fixed_type/${types}.txt
heading="${heading}\t${file%_mean*}"
done < fixed_type/order.tmp
paste fixed_type/${types}.txt description.tmp | sort -k4,4r > tmp && mv tmp fixed_type/${types}.txt
sed -i '1i '${heading}'\tDescription'  fixed_type/${types}.txt

done

## 6.3 Uninduced. Ref: 0h uninduced
mkdir Uninduced
ref="WT_0h_Uninduced_mean.tmp"
ls groups | grep WT.*Uninduced.*.tmp > Uninduced/order.tmp
ls groups | grep Clone.*Uninduced.*.tmp >> Uninduced/order.tmp

paste gene.tmp > Uninduced/Uninduced.txt
heading="Position\tGene"
while read file;do  # WT_24h_Uninduced_mean.tmp
paste groups/${ref} groups/${file} | awk 'BEGIN{FS="\t"}{if($1>0){print $2/$1}else{if("'${file}'"=="'${ref}'"){print 1}else{print "Inf"}}}' > Uninduced/${file%mean*}f.tmp
paste Uninduced/Uninduced.txt Uninduced/${file%mean*}f.tmp > tmp && mv tmp Uninduced/Uninduced.txt
heading="${heading}\t${file%_mean*}"
done < Uninduced/order.tmp
paste Uninduced/Uninduced.txt description.tmp | sort -k4,4r > tmp && mv tmp Uninduced/Uninduced.txt
sed -i '1i '${heading}'\tDescription'  Uninduced/Uninduced.txt

rm -f *.tmp groups/*.tmp fixed_time/*.tmp  fixed_type/*.tmp  Uninduced/*.tmp 

mkdir user_custom