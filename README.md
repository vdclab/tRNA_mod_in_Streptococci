## Description
This project contains commands and scripts that haven been generated and used by the de Crecy lab at University of Florida for codon usage analysis in the manuscript:

“tRNA Modification Landscapes in Streptococci: Shared Losses and Clade-Specific Adaptations”

Authors: Ho-Ching Tiffany Tsui, Chi-Kong Chan, Yifeng Yuan, Roba Elias, Jingjing Sun, Virginie Marchand, Marshall Jaroch, Guangxin Sun, Irfan Manzoor, Ana Kutchuashvili, Yuri Motorin, Grazyna Leszczynska, Kinda Seaton, Kelly C. Rice, Manal A. Swairjo, Malcolm E. Winkler, Peter C. Dedon and Valérie de Crécy-Lagard

## Dependencies
MASH v2.3 https://github.com/marbl/mash

seqkit v2.8.0 https://bioinf.shenwei.me/seqkit/

python v3.12 https://www.python.org/

## Usage
Preparation
1. Retrive CDSs sequences of S. pnuemoniae strains and other Streptococcus strains in fasta format from NCBI (GCF_xxxx.ffn and GCA_xxx.ffn files) and put them in the working directory, for example $dir=work_dir, where GCF_xxx and GCA_xxx are NCBI assembly IDs.

2. Prepare two txt file of NCBI assembly IDs for the S.pnuemoniae and other Streptococcus, for example, list_Spne.txt and list_other.txt.

3. calculate the distance between S. pnuemoniae strains and other Streptococcus genomes using MASH.

```bash
# -p threads
# -l Lines in each <input> specify paths to sequence files, one per line.
mash sketch -p 4 -o Spne_sketch -l list_Spne.txt

mash sketch -p 4 -o other_sketch -l list_other.txt

mash dist -p 4 Spne_sketch.msh other_sketch.msh > mash_result.txt
```

4. Clean result by removing some genomes with queDEC, queF, TGT, queA, queG/H (genome id in remove_Qpattern_gid.txt) and removing duplicates.

```bash
sed 's#working_dir/##g' mash_result.txt | sort -k2,2n > mash_result_sort.txt
grep -v -E '^\s*$' remove_Qpattern_gid.txt | grep -F -v -f - mash_result_sort.txt > mash_result_clean.txt
```

5. Clean the MASH result. Then, calculate the distance between genomes so the sum of distance is minimal.

```
python select_minimal_value_per_gid1.py mash_result_clean.txt mash_res_minimum_by_gidSpne.txt

python select_pair_with_global_minimal_distance.py mash_result_clean.txt mash_res_global_minimum.txt

cut -d$'\t' -f2 mash_result_global_minimum.txt > other_gid.txt

cut -d$'\t' -f2 mash_res_minimum_by_gidSpne.txt >> other_gid.txt
```

6. Count codons in ffn files
```
dir_w=working_dir

mkdir ${dir_w}/ffn_prep || true

# remove CDSs less than 10 amino acids.
# -l, print sequences in lower case.
# -m, print sequences >= 10 aa (33 nts).
for gid in $(cat Spne_gid.txt other_gid.txt) ; do
  seqkit seq -g -l -m 33 ${dir_seq}/${gid%.fna}.ffn > ${dir_w}/ffn_prep/${gid%.fna}.ffn
done

# re-format ffn files.
for gid in $(cat other_gid.txt) ; do
  ffn=${dir_w}/ffn_prep/${gid}.ffn
  # prepare sequences.
  sed -i '/^>/ s/ .*$/#/g' $ffn      # replace space after fig number with #.
  sed -i ':a;N;$!ba;s/\n//g' $ffn    # delet all line breaks.
  sed -i 's/>/\n>/g' $ffn            # make each > a new line.
  sed -i 's/#/\n/g' $ffn             # make # a new line.
  sed -i '/^>/! s/^...//' $ffn       # delete start codon.
  sed -i '/^>/! s/...$//' $ffn       # delete stop codon.
  sed -i '/^>/! s/.\{3\}/& /g' $ffn  # insert space every triplet.
  sed -i '/^$/d' $ffn                # remove empty lines.
done

# count codons.
for gid in $(cat other_gid.txt) ; do
  ffn=${dir_w}/ffn_prep/${gid}.ffn
  python 1_pycodon_count.py ${ffn} ${dir_w}
done

for gid in $(cat other_gid.txt) ; do
  file=${dir_w}/py_count/${gid}.ffn_count.txt
  python3 2_count2CDScodon.py ${file}
done
```

7. Write Q codons' profile to the output file.
```
output=Qcodon_stats.txt
dir_w=working_dir

echo -e "gid\ttotal_aau\ttotal_aac\ttotal_cau\ttotal_cac\ttotal_gau\ttotal_gac\ttotal_uau\ttotal_uac\ttotal_nau\ttotal_nac\ttotal_nnn\tratio_nau/nac" > "${output}"

for file in $(ls ${dir_w}/*_CDScodon.tsv); do
  gid="$(basename "$file" | sed 's/.ffn_CDScodon.tsv//')"
  total_aau=$(awk -F '\t' 'FNR==2{print $5}' "${file}" | paste -sd+ | bc )
  total_aac=$(awk -F '\t' 'FNR==2{print $3}' "${file}" | paste -sd+ | bc )
  total_cau=$(awk -F '\t' 'FNR==2{print $21}' "${file}" | paste -sd+ | bc )
  total_cac=$(awk -F '\t' 'FNR==2{print $19}' "${file}" | paste -sd+ | bc )
  total_gau=$(awk -F '\t' 'FNR==2{print $37}' "${file}" | paste -sd+ | bc )
  total_gac=$(awk -F '\t' 'FNR==2{print $35}' "${file}" | paste -sd+ | bc )
  total_uau=$(awk -F '\t' 'FNR==2{print $53}' "${file}" | paste -sd+ | bc )
  total_uac=$(awk -F '\t' 'FNR==2{print $51}' "${file}" | paste -sd+ | bc )
  total_nau=$(awk -F '\t' 'FNR==2{print $71}' "${file}" )
  total_nac=$(awk -F '\t' 'FNR==2{print $72}' "${file}" )
  total_nnn=$(awk -F '\t' 'FNR==2{print $73}' "${file}" )
  ratio_nau_c=$(echo "scale=6; ${total_nau} / ${total_nac}" | bc)
  echo -e "${gid}\t${total_aau}\t${total_aac}\t${total_cau}\t${total_cac}\t${total_gau}\t${total_gac}\t${total_uau}\t${total_uac}\t${total_nau}\t${total_nac}\t${total_nnn}\t${ratio_nau_c}" >> "${output}"
done
```


## Help and Issues
Please contact Yifeng Yuan at yuanyifeng@ufl.edu

## Authors
Yifeng Yuan, Ph.D. Valerie de Crecy-Lagard (Principal Investigator)

Version History
v1.0 --2024/02/26 --creation of the project.
