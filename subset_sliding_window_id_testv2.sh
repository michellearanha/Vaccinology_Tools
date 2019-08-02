rm *.out
set="117_emm_types"
fasta_all_filename="117_emm.fasta"
fasta_headers_file="fasta_headers_file.txt"

#test ntc6
#declare -a emmtypes=("emm_1_1-50" "emm_2_1-50" "emm_9_1-50" "emm_15_1-50" "emm_49_1-50" "emm_50_1-50" "emm_73_1-50" "emm_77_1-50" "emm_84_1-50" "emm_89_1-50" "emm_102_1-50" "emm_112_1-50" "emm_114_1-50" "emm_118_1-50" "emm_124_1-50" "emm_175_1-50" "emm_183_1-50" "emm_227_1-50" "emm_232_1-50" "emm_238_1-50" "emm_239_1-50")

# test e4
#declare -a emmtypes=("emm_3_1-50" "emm_5_1-50" "emm_11_1-50" "emm_49_1-50" "emm_114_1-50" "emm_102_1-50" "emm_77_1-50" "emm_28_1-50" "emm_232_1-50" "emm_175_1-50" "emm_124_1-50" "emm_169_1-50" "emm_109_1-50" "emm_88_1-50" "emm_84_1-50" "emm_73_1-50" "emm_2_1-50" "emm_8_1-50" "emm_22_1-50" "emm_89_1-50" "emm_112_1-50")

#test small set
#declare -a emmtypes=("emm_1_1-50" "emm_2_1-50" "emm_9_1-50"  "emm_227_1-50" "emm_238_1-50")

#this command does not work because apple will not upgrade to bash4
#mapfile -t emmtypes < $fasta_headers_file 

# take the header names fromt he fasta file 
awk 'NR % 2 == 1' $fasta_all_filename > $fasta_headers_file
# remove the > symbol from the fasta header
sed 's/^.//' $fasta_headers_file >tempfile && mv tempfile $fasta_headers_file


declare -a emmtypes
emmtypes=(`cat "$fasta_headers_file"`)

printf "%s\n" "${emmtypes[@]}" > emm.csv

sed 's/_1-50//g' emm.csv > temp1.csv && mv temp1.csv emm.csv
sed 's/_//g' emm.csv > temp1.csv && mv temp1.csv emm.csv

criteria="Identity"
seq_identity_threshold=45
seq_len=50
window_size=8
gap=3 # gap between windows
n_windows=$((((seq_len - window_size)/gap) + 1))


for x in `seq 1 $n_windows`
do
echo "$x";
done>> x_axis_raw.out
sed -e $'1i\\\nWindow#' x_axis_raw.out > x_axis.out

################### start iteration for each sequence ##########################
for i in "${emmtypes[@]}"
do
cp emm.csv temp.csv

echo "${emmtypes[@]}" 

grep -A1 ""$i"" "$fasta_all_filename" >> "$i".faa


#### added new - for list subset ######
/Users/ma0/software_installed/EMBOSS-6.6.0/bin/bin/needle -asequence "$i".faa -bsequence $fasta_all_filename -sprotein1 Y -sprotein2 Y  stdout -gapopen 10.0 -gapextend 0.5>> "$i".needle

grep $criteria "$i".needle >> "$i"_id_per_raw.dat
awk -F"[()]" '{print $2}' "$i"_id_per_raw.dat >> percent_id_final_"$i".out
paste <(cat emm.csv) <(cat percent_id_final_"$i".out)| column -s $'\t' -t > temp3.csv && mv temp3.csv percent_id_final_"$i".out
sed 's/%//g' percent_id_final_"$i".out >temp3.csv && mv temp3.csv percent_id_final_"$i".out
awk  -v var="${seq_identity_threshold}" '$2>var' percent_id_final_"$i".out  > high_id_w_"$i".out
cp high_id_w_"$i".out temp_high_id.dat
list=( $(awk '{print $1}' ./temp_high_id.dat) )
printf "%s\n" "${list[@]}" > selected_types.txt
rm *.dat
### end new #####

#####iterate over windows######
##################################
for j in `seq 1 $n_windows`
do
l=$((j+(j-1)*(gap-1)))
k=$((l+$window_size-1))
sed -n 2p "$i".faa | cut -c${l}-${k} >"$i"_"window""$j".faa


/Users/ma0/software_installed/EMBOSS-6.6.0/bin/bin/needle -asequence "$i"_"window""$j".faa -bsequence $fasta_all_filename -sprotein1 Y -sprotein2 Y  stdout -gapopen 10.0 -gapextend 0.5>> "$i"_"window""$j".needle

grep "$criteria" "$i"_"window""$j".needle >>"$i"_"window""$j"_raw_matches.dat
cut -c18-19 "$i"_"window""$j"_raw_matches.dat>> "$i"_"window""$j"_correct_matches.dat

paste <(cat temp.csv) <(cat "$i"_"window""$j"_correct_matches.dat) | column -s $'\t' -t >temp2.csv 
mv temp2.csv temp.csv

done 

mv temp.csv "$i"_all_windows.out

#####to transpose output for plots#######
awk '
{
    for (m=1; m<=NF; m++)  {
        a[NR,m] = $m
    }
}
NF>p { p = NF }
END {
    for(n=1; n<=p; n++) {
        str=a[1,n]
        for(m=2; m<=NR; m++){
            str=str" "a[m,n];
        }
        print str
    }
}' "$i"_all_windows.out > "$i"_all_windows_to_plot.out
##### end to transpose #########


paste <(cat x_axis.out) <(cat "$i"_all_windows_to_plot.out) |column -s $'\t' -t > plot_"$i"_all_windows.out
file_to_plot="plot_"$i"_all_windows.out"

#gnuplot  -persist -e "set title '${file_to_plot}' noenhanced ;plot for [t=2:'$((${#emmtypes[@]}+1))'] '${file_to_plot}' using 1:t with lines title columnhead"

#### subset plotting #####
./colnum.sh plot_"$i"_all_windows.out  selected_types.txt > new_"$i".out
subset_file_to_plot="new_"$i".out"
#paste <(cat x_axis.out)<(cat new_"$i".out) | column -s $'\t' -t > temp5.out 
#mv temp5.out new_"$i".out
#gnuplot  -persist -e "set title '${subset_file_to_plot}' noenhanced ;plot for [t=2:'$((${#list[@]}+1))'] '${subset_file_to_plot}' using 1:t with lines title columnhead"

rm *.faa
done

rm *.dat  
#mkdir Folder_needle_files
#mv *.needle Folder_needle_files
mkdir Folder_needle_files_"$set"
mv *.needle Folder_needle_files_"$set"

# testing cat percentid files to make a matrix that shows the overall sequence identity many to many after global alignment
mkdir Matrix_overal_seq_id_many_to_many
cp percent_id_final*.out Matrix_overal_seq_id_many_to_many/
cp emm.csv temp6.csv
echo -e "emm_types" | cat - temp6.csv > /tmp/out && mv /tmp/out temp6.csv
cp temp6.csv Matrix_overal_seq_id_many_to_many/
cd Matrix_overal_seq_id_many_to_many/

for i in "${emmtypes[@]}"
do
awk '{$1=""; print $0}' percent_id_final_""$i"".out > seq_id_many_to_many_""$i"".dat
paste -d " " temp6.csv seq_id_many_to_many_""$i"".dat >  temp7.csv && mv temp7.csv temp6.csv
done

cp temp6.csv OVERALL_SEQ_ID_MANY_TO_MANY.dat
cd ..
rm *.csv


