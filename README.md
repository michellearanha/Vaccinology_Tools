# Vaccinology Package
Determine cross-reactivity between peptide sequences based on short sliding windows and/or helical wheel homology of coiled-coil proteins

# Dependencies:
1. Linux or macOS based operating system
2. EMBOSS program http://emboss.open-bio.org/html/adm/ch01s01.html

# Installation Instructions

**1. Download EMBOSS program and install.**

(A helpful link with EMBOSS installation instructions: https://www.shengweihou.com/blog/install_emboss)

**2. Download script files from repository**

**3. edit ~/.bash_profile **

`vi 	~/.bash_profile
export PATH=$PATH:/directory_containing_script_files
  
source ~/.bash_profile`

**4. (A) To run sliding window**

**slide -f {fasta file with sequences} -w {window size} -g {gap size} -l {length of sequences} -o {output folder}**

**(B) To run heptad homology program**

`heptad_id [-f <file:fasta file with sequences>][-r <file:heptad register file>][-c <String:Criteria-Identity/Similarity>][-t<threshold: Percentage>][-E <Directory_of_EMBOSS_program>] [-A <alignment_algorithm>]`


**Options for heptad_id**
`
-f [<.fasta>] (seq.fasta) (Input)

File with sequences listed in fasta format

-r [<.txt>/<.dat>/...] (register.txt) (Input)

File with a list of heptad registers of all sequences (predicted or actual heptad site of the first residue of each sequence

-c [<type_enum>] (Identity) (Input)
  
Criteria : Identity, Similarity

-t [<float_percentage>] (45) (Optional)
  
threshold percentage:  sequences that share "criteria" greater than threshold will be saved in a separate folder.

-E [<dir_path>] (Input)

Directory of EMBOSS binaries

-A [<type_enum>] (water) (Input)
  
Global alignment - Needleman-Wunsch algorithm : needle

Local alignment - Waterman algorithm : water`


# Example of command to run heptad_id
`heptad_id -f seq.fasta -r register.txt -c Identity -t 25 -E /Users/XXX/EMBOSS-6.6.0/bin -A needle`

`heptad_id -f seq.fasta -r register.txt -c Similarity -t 60 -E /Users/XXX/EMBOSS-6.6.0/bin -A water`


# Tutorial with examples

  

What the program does:

*****Sliding window*******

Scanning for matches between sliding k-mer sequence of vaccine type and entire sequence of non-vaccine types (Sliding window approach)
Consider a set of potentially cross reactive protein sequences that are related within species (e.g. antigenic M proteins of different strains of S.pyogenes ) or between species (e.g. antigenic proteins of different viruses belonging to the Flaviviridae virus family such as Dengue, Zika, West Nile that are known induce cross reactive antibodies). Each sequence in the given set is iteratively made the reference sequence and the tool uses a sliding window approach to calculate a pairwise alignment and the number of residues that are identical between successive fragments of the reference sequence and the sequence of the remaining proteins in the set. Matches in terms of sequence similarity as a criterion can also be chosen instead of sequence identity. The subsequence size/ window length/ fragment size (k) for the reference sequence can be adjusted, and the number of amino acids by which the window moves down the reference sequence can be adjusted by adjusting the gap length (m). This gap also determines the degree of overlap between successive refence sequence fragments. The sliding window direction is from the N-terminal to the C-terminal. The number of identical residues between vaccine and non-vaccine type is plotted as a function of window number. If desired, an overall sequence identity or similarity cutoff can be chosen to limit the number of sequences compared. 

*********Helical wheel homology***********
A standard or canonical coiled coil structure is built by two or more helices twisting around each other forming bundles with their side chains interlocking in a ‘knobs’ into ‘holes’ packing. The regular meshing of knobs into holes packing requires recurrent positions of the side-chains every seven residues along the helix interface. This seven-residue sequence repeat is called a heptad repeat and the positions in the heptad repeat are labelled a-g. The core forming positions (a and d) are usually occupied by hydrophobic residues whereas the remaining, solvent exposed positions (b, c, e, f and g) are dominated by hydrophilic residues. For proteins with coiled-coil secondary structure, machine learned models implemented in webservers exist that can assign the residues to the heptad pattern seen in coiled-coil proteins. We compared the assignment of individual residues in a coiled coil sequence to the heptad by the program MARCOIL (website) and found that it was identical to the actual heptad assignment of three M protein coiled coil crystal structures (pdb ids: 2oto, 5hyt, 5hzp) that were of interest to us. Similar assignments can also be done using NCOILS, PCOILS. 
The tool requires three inputs 1) a text file with the sequences of all vaccine and non-vaccine types in the fasta format, 2) a text file with identifiers/headers of the fasta sequences without ‘>’ character 3) a file with the heptad registers of all sequences (assigned heptad position of the first residue of all sequences) which can be obtained using MARCOIL. The identity at corresponding heptad positions of each peptide/protein with every other peptide/protein in the given set is output in a tabular format. Additionally, an output is also provided as a matrix in tabular format that contains heptad identity of every sequence with every other sequence which can be clustered to obtain clusters of cross-reactive peptides. An R-based script for clustering is also provided. 

