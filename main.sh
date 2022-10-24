#!/usr/bin/bash

Threads=50
Read_pairs_num=100000  # the number of read-pairs for each sample

echo "copying the raw experimental RNAseq folder to the current directory..."
cp -r /localdisk/data/BPSM/ICA1/fastq . 

##################################################################################
###################### Code block 1 : Integrity check ############################
##################################################################################
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

##################################################################################
##################### Code block 2 : Quality check ###############################
##################################################################################
mkdir -p report # fastqc does not create the folder automatically
fastqc fastq/*fq.gz -t ${Threads} -o report --extract  # perform fastqc and output to report folder, unzipping all the files.
rm -f summary_per_item.tmp summary_per_file.tmp
# get the quality summary for all the seqs once a time horizontally, read the items the fastqc has checked(e.g.Per base sequence content) in loop
cat report/*fastqc/summary.txt | cut -f2 | sort | uniq | while read item;do
  # from all the quality summary, get how many FAIL/WARN/PASS are there for each checked item, and write to a temp file.
  fail=$(grep "${item}" report/*fastqc/summary.txt | grep -c FAIL)
  warn=$(grep "${item}" report/*fastqc/summary.txt | grep -c WARN)
  pass=$(grep "${item}" report/*fastqc/summary.txt | grep -c PASS)
  echo -e "${item}\t${fail}\t${warn}\t${pass}" >> summary_per_item.tmp
done
echo -e "Items\tFAIL\tWARN\tPASS" > summary_per_item.txt
# sort in descending order according to the number of FAIL, then WARN; write to the final file.
sort -t$'\t' -k2,2r -k3,3r summary_per_item.tmp >> summary_per_item.txt

# the similar process as above, but read the seq files in loop, and list how many FAIL/WARN/PASS for each seq file.
ls report | grep fastqc$ | while read file;do
  fail=$(grep -c FAIL report/${file}/summary.txt)
  warn=$(grep -c WARN report/${file}/summary.txt)
  pass=$(grep -c PASS report/${file}/summary.txt)
  echo -e "${file}\t${fail}\t${warn}\t${pass}" >> summary_per_file.tmp
done
echo -e "Seq_file\tFAIL\tWARN\tPASS" > summary_per_file.txt
sort -t$'\t' -k2,2r -k3,3r summary_per_file.tmp >> summary_per_file.txt

rm -f summary_per_item.tmp summary_per_file.tmp

##################################################################################
##################### Code block 3 : Alignment ###################################
##################################################################################
cp -r /localdisk/data/BPSM/ICA1/Tcongo_genome/ .  # copy the whole genome of Trypanosoma congolense
# build the index for the genome, output the index files with a prefix of "tc"
bowtie2-build Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz  Tcongo_genome/tc -- threads ${Threads}
mkdir -p Align_bam
# get all the end1 and end2 seq file respectively
ls fastq | grep 1.fq.gz > fq1.tmp
ls fastq | grep 2.fq.gz > fq2.tmp
# paste them so that bowtie2 can recursively align each paired-end data with whole genome.
paste fq1.tmp fq2.tmp | while read fq1 fq2;do
  # Using 50 threads. the result sam files are then transformed to bam files with -b,(-S:improve compatibility). Bam files are sorted so as to improve the speed for following steps.  
  bowtie2 -p ${Threads} -x Tcongo_genome/tc -1 fastq/${fq1} -2 fastq/${fq2} | samtools view -bS | samtools sort > Align_bam/${fq1:0:8}.bam
  # output files are stored with the name of e.g. Tco-5053.bam
done
rm -f fq1.tmp fq2.tmp

##################################################################################
###################### Code block 4 : Counts data ################################
##################################################################################
cp /localdisk/data/BPSM/ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed . # get the gene information file
mkdir -p counts
ls Align_bam | while read bamfile;do
  # for each seq file, generate how many overlaps do each gene in bed file have. This will append colunms containing alignment information to the original bed file as output files
  bedtools coverage -a TriTrypDB-46_TcongolenseIL3000_2019.bed -b Align_bam/${bamfile} > counts/${bamfile:0:8}cov.out
  # output files are stored with the name of e.g. Tco-5053cov.out
done

##################################################################################
###################### Code block 5 : Mean for groups ############################
##################################################################################
mkdir -p groups
rm -f groups/*.tmp
######## 5.1 some preparation: generate group files and get group info and others.
# print the group info file excluding the head line, and sort. k2:type,k4:time,k5:treatment,k3:replicate. Sorted file will be at the order of groups with replicate from 1 to n
tail -n +2 fastq/Tco.fqfiles | sort -k2,2 -k4,4 -k5,5 -k3,3 | while read sname stype rep time treat end1 end2;do
  # encounter replicate 1, means this is a line starting a new group against to lines above. 
  if [ ${rep} == 1 ];then
    # paste a filename according to this group parameter and then create a file for the group.
    filename=${stype}_${time}h_${treat}  # e.g WT_24h_Uninduced
    touch groups/${filename}_file.tmp # e.g. WT_24h_Uninduced_file.tmp
  fi
  # get and slice the seq filename of this line. all the filenames from the same group will end up in the same group file.
  # e.g. WT_24h_Uninduced_file.tmp contains 3 lines: Tco-6890 Tco-6600 Tco-6037
  echo "${end1:0:8}" >> groups/${filename}_file.tmp
done
# get the gene position, uid, description column from bed file, for whenever need to paste into other files later.
cut -f1,4 TriTrypDB-46_TcongolenseIL3000_2019.bed > gene.tmp
cut -f5 TriTrypDB-46_TcongolenseIL3000_2019.bed > description.tmp

######### 5.2 calculate the mean count/exp level.
replicate=0
rm -f groups/rep*
# read the group file containing which Tco- are in that group, in loop.
ls groups | grep file.tmp | while read groupfile;do  # e.g. WT_24h_Uninduced_file.tmp 
  # for a specific group, read the seq filename in that group in loop
  cat groups/${groupfile} | while read seq_num;do  # e.g.: Tco-5053
    ((replicate+=1))
    # each file(replication) will have a unique filename for its count data, by adding a variable "replicate"
    # count data are in 6th colunm of the outputfile.
    cut -f6 counts/${seq_num}cov.out > groups/rep${replicate}.tmp
  done
  # for a specific group, paste the count data column by column, calculate the mean.
  paste groups/rep* | awk 'BEGIN{sum=0}{for(i=1;i<=NF;i++){sum+=$i}{print sum/NF; sum=0}}' > groups/${groupfile%file*}mean.tmp # outputfile e.g.: WT_24h_Uninduced_mean.tmp
  # add the gene and description columns to the final txt file. 
  paste gene.tmp groups/${groupfile%file*}mean.tmp description.tmp > groups/${groupfile%file*}mean.txt
  rm -f groups/rep*
done
# add a heading line to all the final txt file
sed -i '1i Position\tGene\tMean_exp\tDescription'  groups/*mean.txt
rm -f groups/*file.tmp

##################################################################################
###################### Code block 6 : Fold change comparisons ####################
##################################################################################
######### 6.1 Specific time points. WT Uninduced group is set as reference
dirname="Fold_change/fixed_time"
mkdir -p ${dirname}
# This is to get all the time point the experiments have done, from the info file Tco.fqfiles
tail -n +2 fastq/Tco.fqfiles | cut -f4 | sort | uniq | while read time;do  # $time=0,24,48...
  ref="WT_${time}h_Uninduced_mean.tmp" # indentify the reference group, whose data will be nomorlised to 1.
  # for each time point, generate a file that indicates the order of the groups. This is the column order shown in final output file. e.g. for 0h groups: WT_0h_Uninduced(ref),Clone1_0h_Uninduced,Clone1_0h_Induced,Clone2_0h_Uninduced,Clone2_0h_Induced,
  ls groups | grep WT_${time}h_.*.tmp | sort -t"_" -k3,3r > ${dirname}/order.tmp
  ls groups | grep Clone.*_${time}h_.*.tmp | sort -t"_" -k1,1 -k3,3r >> ${dirname}/order.tmp
  # initialise the final output file, first paste gene info columns into it. 
  paste gene.tmp > ${dirname}/${time}h.txt
  heading="Position\tGene"  # to initialise the heading line, setting a variable.
  # read the mean count files from a specific group in order from order.tmp file, in loop
  while read file;do  # e.g.: WT_0h_Uninduced_mean.tmp
    # the column 1 of the pasted text is always the mean data of reference group. if the file in loop is exactly the reference, then all the data greater than 0 will be transformed to 1. if the count is 0 in reference group, then it will be converted to a sufficient small number(0.0001) to allow the calculation of fold change.
    paste groups/${ref} groups/${file} | awk 'BEGIN{FS="\t"}{if($1==0){$1=0.0001}{print $2/$1}}' > ${dirname}/${file%mean*}f.tmp # outputfile e.g. Clone1_0h_Uninduced_f.tmp
    # append the fold change data column to the final file.
    paste ${dirname}/${time}h.txt ${dirname}/${file%mean*}f.tmp > tmp && mv tmp ${dirname}/${time}h.txt
    heading="${heading}\t${file%_mean*}" # append the group name to heading variable.
  done < ${dirname}/order.tmp
  # add the description to the final file, and sort in descending order of the 4th column(in this case is Clone1_0h_Uninduced)
  paste ${dirname}/${time}h.txt description.tmp | sort -k4,4nr > tmp && mv tmp ${dirname}/${time}h.txt
  sed -i '1i '${heading}'\tDescription'  ${dirname}/${time}h.txt # add a heading line
done  # time loop


######### 6.2 Specific types. ref group: 0h_Uninduced
# The total process is exactly the same as 6.1, but here the type is in loop.
dirname="Fold_change/fixed_type"
mkdir -p ${dirname}
tail -n +2 fastq/Tco.fqfiles | cut -f2 | sort | uniq | while read types;do  # WT, Clone1,Clone2,Clonen
  ref="${types}_0h_Uninduced_mean.tmp"
  ls groups | grep ${types}.*.tmp | sort -t"_" -k2,2 -k3,3r > ${dirname}/order.tmp
  paste gene.tmp > ${dirname}/${types}.txt
  heading="Position\tGene"
  while read file;do
    paste groups/${ref} groups/${file} | awk 'BEGIN{FS="\t"}{if($1==0){$1=0.0001}{print $2/$1}}' > ${dirname}/${file%mean*}f.tmp
    paste ${dirname}/${types}.txt ${dirname}/${file%mean*}f.tmp > tmp && mv tmp ${dirname}/${types}.txt
    heading="${heading}\t${file%_mean*}"
  done < ${dirname}/order.tmp
  paste ${dirname}/${types}.txt description.tmp | sort -k4,4nr > tmp && mv tmp ${dirname}/${types}.txt
  sed -i '1i '${heading}'\tDescription'  ${dirname}/${types}.txt
done # type loop

########## 6.3 Uninduced. ref group: WT_0h
# The whole process is same as above, but here does not include a big loop to search the multiple fixed groups like time points or types ("Uninduced" is specified).
dirname="Fold_change/Uninduced"
mkdir -p ${dirname}
ref="WT_0h_Uninduced_mean.tmp"
ls groups | grep WT.*Uninduced.*.tmp > ${dirname}/order.tmp
ls groups | grep Clone.*Uninduced.*.tmp >> ${dirname}/order.tmp
paste gene.tmp > ${dirname}/Uninduced.txt
heading="Position\tGene"
while read file;do
  paste groups/${ref} groups/${file} | awk 'BEGIN{FS="\t"}{if($1==0){$1=0.0001}{print $2/$1}}' > ${dirname}/${file%mean*}f.tmp
  paste ${dirname}/Uninduced.txt ${dirname}/${file%mean*}f.tmp > tmp && mv tmp ${dirname}/Uninduced.txt
  heading="${heading}\t${file%_mean*}"
done < ${dirname}/order.tmp
paste ${dirname}/Uninduced.txt description.tmp | sort -k4,4nr > tmp && mv tmp ${dirname}/Uninduced.txt
sed -i '1i '${heading}'\tDescription'  ${dirname}/Uninduced.txt

rm -f *.tmp groups/*.tmp Fold_change/fixed_time/*.tmp  Fold_change/fixed_type/*.tmp  Fold_change/Uninduced/*.tmp 