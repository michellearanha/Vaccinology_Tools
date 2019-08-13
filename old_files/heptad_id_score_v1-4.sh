rm *.faa *.needle *.dat *.out *.csv
mkdir RESULTS
fasta_all_filename="ntc6_emm_marcoil.fasta"
fasta_headers_file="ntc6_fasta_headers_file.txt"
heptad_register_file="ntc6_register.txt"


# take the header names fromt he fasta file
awk 'NR % 2 == 1' $fasta_all_filename > $fasta_headers_file
# remove the > symbol from the fasta header
sed 's/^.//' $fasta_headers_file >tempfile && mv tempfile $fasta_headers_file


criteria1="Similarity"
criteria2="Identity"
seq_identity_threshold=45
seq_len=50
window_size=15
gap=1 # gap between windows
n_windows=$((((seq_len - window_size)/gap) + 1))

#array contains sequence headers as in fasta file
#declare -a emmtypes=("emm_1_1-50" "emm_2_1-50" "emm_9_1-50" "emm_15_1-50" "emm_49_1-50" "emm_50_1-50" "emm_73_1-50" "emm_77_1-50" "emm_84_1-50" "emm_89_1-50" "emm_102_1-50" "emm_112_1-50" "emm_114_1-50" "emm_118_1-50" "emm_124_1-50" "emm_175_1-50" "emm_183_1-50" "emm_227_1-50" "emm_232_1-50" "emm_238_1-50" "emm_239_1-50")

declare -a emmtypes
emmtypes=(`cat "$fasta_headers_file"`)

#array contains start heptad register of each sequence, each number denotes a heptad alphabet 
#1=a, 2=b, 3=c, 4=d, 5=e, 6=f, 7=g 
#declare -a heptad_register=("1" "7" "5" "1" "7" "1" "6" "1" "1" "1" "5" "4" "6" "1" "6" "5" "2" "5" "5" "1" "1")

declare -a heptad_register
heptad_register=(`cat "$heptad_register_file"`)

#array contains heptad positions
declare -a heptad=("a" "b" "c" "d" "e" "f" "g")

printf "%s\n" "${emmtypes[@]}" > emm.csv
sed 's/_1-50//g' emm.csv > temp1.csv && mv temp1.csv emm.csv
sed 's/_//g' emm.csv > temp1.csv && mv temp1.csv emm.csv



for i in "${heptad[@]}"
do

touch ""$i"".faa

done

#Code logic- Give heptad register for each emm type and then extract sequence at each heptad for the sequence and store in the fasta file for that heptad for instance in d.faa
for z in "${emmtypes[@]}"
do

grep -A1 ""$z"" "$fasta_all_filename" >> "$z".faa
s=`sed -n '2p' "$z".faa`;


#https://unix.stackexchange.com/questions/161958/extract-every-nth-character-from-a-string
temp_heptad_array=($(
    for ((i=0; i<7; i++))
    do
        new=${s:$i:1}
        for ((j=i+7; j<${#s}; j=j+7))
        do 
            new="$new${s:$j:1}"
        done
        echo "$new"
    done
    ))

echo "${temp_heptad_array[@]}"


value=$z

for q in "${!emmtypes[@]}"; do
   if [[ "${emmtypes[$q]}" = "${value}" ]]; then
       echo "${q}";
       		if [ ${heptad_register[$q]} == 1 ]
		then
		echo ">${emmtypes[$q]}" >> a.faa
		echo "${temp_heptad_array[0]}" >> a.faa
		echo ">${emmtypes[$q]}" >> b.faa
                echo "${temp_heptad_array[1]}" >> b.faa
		echo ">${emmtypes[$q]}" >> c.faa
                echo "${temp_heptad_array[2]}" >> c.faa
		echo ">${emmtypes[$q]}" >> d.faa
                echo "${temp_heptad_array[3]}" >> d.faa
		echo ">${emmtypes[$q]}" >> e.faa
                echo "${temp_heptad_array[4]}" >> e.faa
		echo ">${emmtypes[$q]}" >> f.faa
                echo "${temp_heptad_array[5]}" >> f.faa
		echo ">${emmtypes[$q]}" >> g.faa
                echo "${temp_heptad_array[6]}" >> g.faa
		fi
   	
		if [ ${heptad_register[$q]} == 2 ]
                then
                echo ">${emmtypes[$q]}" >> b.faa
                echo "${temp_heptad_array[0]}" >> b.faa
                echo ">${emmtypes[$q]}" >> c.faa
                echo "${temp_heptad_array[1]}" >> c.faa
		echo ">${emmtypes[$q]}" >> d.faa
                echo "${temp_heptad_array[2]}" >> d.faa
		echo ">${emmtypes[$q]}" >> e.faa
                echo "${temp_heptad_array[3]}" >> e.faa
		echo ">${emmtypes[$q]}" >> f.faa
                echo "${temp_heptad_array[4]}" >> f.faa
		echo ">${emmtypes[$q]}" >> g.faa
                echo "${temp_heptad_array[5]}" >> g.faa
		echo ">${emmtypes[$q]}" >> a.faa
                echo "${temp_heptad_array[6]}" >> a.faa
		fi
		
		if [ ${heptad_register[$q]} == 3 ]
                then
                echo ">${emmtypes[$q]}" >> c.faa
                echo "${temp_heptad_array[0]}" >> c.faa
                echo ">${emmtypes[$q]}" >> d.faa
                echo "${temp_heptad_array[1]}" >> d.faa
		echo ">${emmtypes[$q]}" >> e.faa
                echo "${temp_heptad_array[2]}" >> e.faa
		echo ">${emmtypes[$q]}" >> f.faa
                echo "${temp_heptad_array[3]}" >> f.faa
		echo ">${emmtypes[$q]}" >> g.faa
                echo "${temp_heptad_array[4]}" >> g.faa
		echo ">${emmtypes[$q]}" >> a.faa
                echo "${temp_heptad_array[5]}" >> a.faa
		echo ">${emmtypes[$q]}" >> b.faa
                echo "${temp_heptad_array[6]}" >> b.faa
		fi

		if [ ${heptad_register[$q]} == 4 ]
                then
                echo ">${emmtypes[$q]}" >> d.faa
                echo "${temp_heptad_array[0]}" >> d.faa
		echo ">${emmtypes[$q]}" >> e.faa
                echo "${temp_heptad_array[1]}" >> e.faa                
		echo ">${emmtypes[$q]}" >> f.faa
                echo "${temp_heptad_array[2]}" >> f.faa
		echo ">${emmtypes[$q]}" >> g.faa
                echo "${temp_heptad_array[3]}" >> g.faa
		echo ">${emmtypes[$q]}" >> a.faa
                echo "${temp_heptad_array[4]}" >> a.faa
		echo ">${emmtypes[$q]}" >> b.faa
                echo "${temp_heptad_array[5]}" >> b.faa
		echo ">${emmtypes[$q]}" >> c.faa
                echo "${temp_heptad_array[6]}" >> c.faa
		fi		
		
		if [ ${heptad_register[$q]} == 5 ]
                then
                echo ">${emmtypes[$q]}" >> e.faa
                echo "${temp_heptad_array[0]}" >> e.faa
		echo ">${emmtypes[$q]}" >> f.faa
                echo "${temp_heptad_array[1]}" >> f.faa
                echo ">${emmtypes[$q]}" >> g.faa
                echo "${temp_heptad_array[2]}" >> g.faa
		echo ">${emmtypes[$q]}" >> a.faa
                echo "${temp_heptad_array[3]}" >> a.faa
		echo ">${emmtypes[$q]}" >> b.faa
                echo "${temp_heptad_array[4]}" >> b.faa
		echo ">${emmtypes[$q]}" >> c.faa
                echo "${temp_heptad_array[5]}" >> c.faa
		echo ">${emmtypes[$q]}" >> d.faa
                echo "${temp_heptad_array[6]}" >> d.faa
		fi

		if [ ${heptad_register[$q]} == 6 ]
                then
                echo ">${emmtypes[$q]}" >> f.faa
                echo "${temp_heptad_array[0]}" >> f.faa
		echo ">${emmtypes[$q]}" >> g.faa
                echo "${temp_heptad_array[1]}" >> g.faa                
		echo ">${emmtypes[$q]}" >> a.faa
                echo "${temp_heptad_array[2]}" >> a.faa
		echo ">${emmtypes[$q]}" >> b.faa
                echo "${temp_heptad_array[3]}" >> b.faa
		echo ">${emmtypes[$q]}" >> c.faa
                echo "${temp_heptad_array[4]}" >> c.faa
		echo ">${emmtypes[$q]}" >> d.faa
                echo "${temp_heptad_array[5]}" >> d.faa		
		echo ">${emmtypes[$q]}" >> e.faa
                echo "${temp_heptad_array[6]}" >> e.faa		
		fi
	
		if [ ${heptad_register[$q]} == 7 ]
                then
                echo ">${emmtypes[$q]}" >> g.faa
                echo "${temp_heptad_array[0]}" >> g.faa
		echo ">${emmtypes[$q]}" >> a.faa
                echo "${temp_heptad_array[1]}" >> a.faa                
		echo ">${emmtypes[$q]}" >> b.faa
                echo "${temp_heptad_array[2]}" >> b.faa
		echo ">${emmtypes[$q]}" >> c.faa
                echo "${temp_heptad_array[3]}" >> c.faa
		echo ">${emmtypes[$q]}" >> d.faa
                echo "${temp_heptad_array[4]}" >> d.faa
		echo ">${emmtypes[$q]}" >> e.faa
                echo "${temp_heptad_array[5]}" >> e.faa
		echo ">${emmtypes[$q]}" >> f.faa
                echo "${temp_heptad_array[6]}" >> f.faa
		fi

	fi
   done



unset ${temp_heptad_array[@]}
done



for i in "${emmtypes[@]}"
do
cp emm.csv temp.csv
touch HEPTAD_ID_""$i"".dat

#for 1 line of context After the match
#####################added_for_subset_makes_a_seleted_types_file_each_M_sequence##############################
/Users/ma0/software_installed/EMBOSS-6.6.0/bin/bin/needle -asequence "$i".faa -bsequence $fasta_all_filename -sprotein1 Y -sprotein2 Y  stdout -gapopen 10.0 -gapextend 0.5>> "$i".needle
grep $criteria2 "$i".needle >> "$i"_id_per_raw.dat
awk -F"[()]" '{print $2}' "$i"_id_per_raw.dat >> percent_"$criteria2"_final_"$i".out
paste <(cat emm.csv) <(cat percent_"$criteria2"_final_"$i".out)| column -s $'\t' -t > temp3.csv && mv temp3.csv percent_"$criteria2"_final_"$i".out
sed 's/%//g' percent_"$criteria2"_final_"$i".out >temp3.csv && mv temp3.csv percent_"$criteria2"_final_"$i".out
awk  -v var="${seq_identity_threshold}" '$2>var' percent_"$criteria2"_final_"$i".out  > high_id_w_"$i".out
cp high_id_w_"$i".out temp_high_id.dat
list=( $(awk '{print $1}' ./temp_high_id.dat) )
printf "%s\n" "${list[@]}" > selected_types.txt
rm temp_high_id.dat
################################################################

########################for heptad position to position comparison #####################
for j in "${heptad[@]}"
do

grep -A1 ""$i"" "$j".faa >> "$i"_"$j".faa

/Users/ma0/software_installed/EMBOSS-6.6.0/bin/bin/needle -asequence "$i"_"$j".faa -bsequence "$j".faa -sprotein1 Y -sprotein2 Y  stdout -gapopen 10.0 -gapextend 0.5>> "$i"_"$j".needle
grep "$criteria2" "$i"_"$j".needle >> "$i"_"$j"_"$criteria2"_w_all_raw.dat
awk -F"[()]" '{print $2}' "$i"_"$j"_"$criteria2"_w_all_raw.dat >> percent_"$criteria2"_final_"$i"_"$j".out
sed 's/%//g' percent_"$criteria2"_final_"$i"_"$j".out > temp22.csv && mv temp22.csv percent_"$criteria2"_final_"$i"_"$j".out

#########################################################################################
paste <(cat temp.csv) <(cat percent_"$criteria2"_final_"$i"_"$j".out) | column -s $'\t' -t > temp3.csv && mv temp3.csv temp.csv 
done
mv temp.csv HEPTAD_ID_""$i"".dat
echo -e "emm_types\ta."$i"\tb."$i"\tc."$i"\td."$i"\te."$i"\tf."$i"\tg."$i"" | cat - HEPTAD_ID_""$i"".dat > /tmp/out && mv /tmp/out HEPTAD_ID_""$i"".dat
grep -w -F -f selected_types.txt HEPTAD_ID_""$i"".dat > THRESHOLD_SELECTED_HEPTAD_ID_""$i"".dat
mv HEPTAD_ID_""$i"".dat THRESHOLD_SELECTED_HEPTAD_ID_""$i"".dat RESULTS/
done


cp emm.csv temp5.csv
echo -e "emm_types" | cat - temp5.csv > /tmp/out && mv /tmp/out temp5.csv
cp temp5.csv RESULTS/
cd RESULTS/

for i in "${emmtypes[@]}"
do
awk '{$1=""; print $0}' HEPTAD_ID_""$i"".dat > raw_data_to_cluster_""$i"".dat
paste -d " " temp5.csv raw_data_to_cluster_""$i"".dat >  temp6.csv && mv temp6.csv temp5.csv
done




#for i in "${emmtypes[@]}"
#do 
#awk '{$1=""; print $0}' HEPTAD_ID_""$i"".dat > raw_data_to_cluster_""$i"".dat 
#paste <(cat temp5.csv) <(cat raw_data_to_cluster_""$i"".dat) | column -s $'\t' -t > temp6.csv && mv temp6.csv temp5.csv
#done

mv temp5.csv HEPTAD_PERCENT_ID_FEATURES_TO_CLUSTER.dat
rm raw*.dat
