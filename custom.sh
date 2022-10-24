#!/usr/bin/bash

# if the folder not exists, create a new folder.
if [ ! -d "user_custom" ];then
mkdir user_custom
fi

####################### User personalisation ###########################
rm -f *.tmp
# list all the group name from groups folder, and add a number in front of them, display to user.
ls groups | awk 'BEGIN{FS="_";OFS="_";num=0}{num=num+1; print num " " $1,$2,$3}' > all_sample.tmp
group_num=$(cat all_sample.tmp | wc -l) # get the number of all the available groups
while true; do  # ask user to choose a reference group
  cat all_sample.tmp
  read -r -p "Which group would you like to be the reference? " answer
  # if the input is a number, less than number of groups and greater than 0
  if [[ ${answer} =~ [0-9]+ ]] && [[ ${answer} -le ${group_num} ]] && [[ ${answer} -gt 0 ]];then
      Ref=$(sed -n ${answer}p all_sample.tmp)  # store the choice to a variable
      break
  else    
     echo "Please type a valid number."
  fi
done
while true; do # ask user to choose comparison groups, same as above.
  cat all_sample.tmp
  read -r -p "Which groups would you like to compare with? (press q to finish) " answer
  if [[ ${answer} =~ [0-9]+ ]] && [[ ${answer} -le ${group_num} ]] && [[ ${answer} -gt 0 ]];then
    # store the chosen group to a file.
    sed -n ${answer}p all_sample.tmp >> choice.tmp
  else if [[ ${answer} =~ [Qq] ]];then
      break
    else
      echo "Please type a valid number."
    fi
  fi
done
chosen_num=$(cat choice.tmp | wc -l)  # the number of chosen comparison groups
while true; do # ask user to choose a comparison group to sort
  cat choice.tmp | awk 'BEGIN{num=0}{num=num+1; print num " " $2}' 
  read -r -p "Which one would you like to sort by descending?" answer
  if [[ ${answer} =~ [0-9]+ ]] && [[ ${answer} -le ${chosen_num} ]] && [[ ${answer} -gt 0 ]];then
    Sortby=${answer}
    break
  else    
   echo "Please type a valid number."
  fi
done

Ref=${Ref#* } # eliminate the number before the group name
Filename=${Ref}${Sortby} # initialise the filename
count=0
unset Com
# read the chosen groups
while read num group;do
  ((count+=1))
  Filename="${Filename}_${num}"
  # store each chosen group to an array "Com", staring from index 1.
  Com[${count}]=${group}
done < choice.tmp
Filename="${Filename}.txt" # final file name generated, according to user's selection above.

heading="Position\tGene\t${Ref}" # initialise the heading, the reference group is at the 3rd column.
# initialise the final output file. get the mean count file of reference group, and set the counts as 1.
tail -n +2  groups/${Ref}_mean.txt | cut -f1,2,3 | awk 'BEGIN{FS="\t";OFS="\t"}{$3=1;print $0}' > user_custom/${Filename}
tail -n +2  groups/${Ref}_mean.txt | cut -f4 > user_custom/description.tmp # get the description column
tail -n +2  groups/${Ref}_mean.txt | cut -f3 > user_custom/ref.tmp # get the original mean count of the reference group.

# read all the chosen groups in loop, 
for group in ${Com[@]};do  # e.g.: WT_0h_Uninduced
  # get the mean count of that group to a temp file.
  tail -n +2  groups/${group}_mean.txt | cut -f3 > user_custom/${group}.tmp
  # paste to reference group mean count, and calculate fold change. 
  paste user_custom/ref.tmp user_custom/${group}.tmp | awk 'BEGIN{FS="\t"}{if($1==0){$1=0.0001}{print $2/$1}}' > user_custom/${group}_f.tmp
  heading="${heading}\t${group}" # connect the group name to the heading.
  # add the calculated data to the final output file
  paste user_custom/${Filename} user_custom/${group}_f.tmp > tmp && mv tmp user_custom/${Filename}
done
# the sort index. the Sortby is taken from user input, it is the number of a comparision group. there are 3 columns ahead of it(2 gene info columns and a reference).
sindex=$((${Sortby}+3))
# paste the description column, and sort by the chosen column(a comparison group).
paste user_custom/${Filename} user_custom/description.tmp | sort -k${sindex},${sindex}nr > tmp && mv tmp user_custom/${Filename}
sed -i '1i '${heading}'\tDescription'  user_custom/${Filename}  # add a heading
echo "Complete. The result was generated in /user_custom/${Filename}"

rm -f user_custom/*.tmp *.tmp