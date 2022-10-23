#!/usr/bin/bash

####################### User personalisation ###########################
rm -f *.tmp
ls groups | awk 'BEGIN{FS="_";OFS="_";num=0}{num=num+1; print num " " $1,$2,$3}' > all_sample.tmp
group_num=$(cat all_sample.tmp | wc -l)
while true; do
    cat all_sample.tmp
    read -r -p "Which group would you like to be the reference? " answer
    if [[ ${answer} =~ [0-9]+ ]] && [[ ${answer} -le ${group_num} ]] && [[ ${answer} -gt 0 ]];then
        Ref=$(sed -n ${answer}p all_sample.tmp)
        break
    else    
       echo "Please type a valid number."
    fi
done
while true; do 
      cat all_sample.tmp
      read -r -p "Which groups would you like to compare with? (press q to finish) " answer
      if [[ ${answer} =~ [0-9]+ ]] && [[ ${answer} -le ${group_num} ]] && [[ ${answer} -gt 0 ]];then
        sed -n ${answer}p all_sample.tmp >> choice.tmp
      else if [[ ${answer} =~ [Qq] ]];then
      break
      else
        echo "Please type a valid number."
      fi
      fi
done
chosen_num=$(cat choice.tmp | wc -l)
while true; do
      cat choice.tmp | awk 'BEGIN{num=0}{num=num+1; print num " " $2}' 
      read -r -p "Which one would you like to sort by descending?" answer
      if [[ ${answer} =~ [0-9]+ ]] && [[ ${answer} -le ${chosen_num} ]] && [[ ${answer} -gt 0 ]];then
        Sortby=${answer}
        break
      else    
       echo "Please type a valid number."
      fi
done

Ref=${Ref#* }
Filename=${Ref}${Sortby}
count=0
while read num group;do
((count+=1))
Filename="${Filename}${num}"
Com[${count}]=${group}
done < choice.tmp
Filename="${Filename}.txt"

heading="Position\tGene\t${Ref}"
tail -n +2  groups/${Ref}_mean.txt | cut -f1,2,3 | awk 'BEGIN{FS="\t";OFS="\t"}{$3=1;print $0}' > user_custom/${Filename}
tail -n +2  groups/${Ref}_mean.txt | cut -f4 > user_custom/description.tmp
tail -n +2  groups/${Ref}_mean.txt | cut -f3 > user_custom/ref.tmp

for group in ${Com[@]};do
tail -n +2  groups/${group}_mean.txt | cut -f3 > user_custom/${group}.tmp
paste user_custom/ref.tmp user_custom/${group}.tmp | awk 'BEGIN{FS="\t"}{if($1>0){print $2/$1}else{print "Inf"}}' > user_custom/${group}_f.tmp
heading="${heading}\t${group}"
paste user_custom/${Filename} user_custom/${group}_f.tmp > tmp && mv tmp user_custom/${Filename}
done
sindex=$((${Sortby}+3))
paste user_custom/${Filename} user_custom/description.tmp | sort -k${sindex},${sindex}r > tmp && mv tmp user_custom/${Filename}
sed -i '1i '${heading}'\tDescription'  user_custom/${Filename}
rm -f user_custom/*.tmp *.tmp

