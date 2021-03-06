# Vaccinology Tools

Determine cross-reactivity between peptide sequences based on short sliding windows and/or helical wheel homology of coiled-coil proteins

## Dependencies:

1. Linux or macOS based operating system
2. [EMBOSS program](http://emboss.open-bio.org/html/adm/ch01s01.html)
 - **_NOTE:_** The download instructions above include an FTP site. If you are unable to download from websites via FTP, try [this fedora repo](https://src.fedoraproject.org/lookaside/extras/EMBOSS/) instead.

## Installation Instructions for MacOS

These instructions are assuming you are compiling and using this software in your home directory, and that you have wget and git installed. If you are installing somewhere else, you will need to adjust the commands accordingly. If you do not have git and wget installed, please check out the documentation to install [homebrew](https://brew.sh/) for macos, as it is the simplest way to get them. 

### 1. Clone this git repo

- `git clone https://github.com/michellearanha/Vaccinology_Tools.git`

### 2. Download and Install EMBOSS.

- Download EMBOSS via [ftp](http://emboss.open-bio.org/html/adm/ch01s01.html) or using the [fedora repository](https://src.fedoraproject.org/lookaside/extras/EMBOSS/)
- Untar the tarball: `tar -xvf EMBOSS-6.6.0.tar.gz`
- In the new EMBOSS directory, run `./configure --prefix=/Users/$('whoami')/EMBOSS-6.6.0/`
- Run `make`
- Run `make install`
- Make some coffee, the last two commands may take a bit to install. 
 
### 3. Edit your ~/.bash_profile

- Need to include this line: `PATH=$PATH:</directory/containing/Vaccinology_Tools/Heptad_homology>` 
	- Assuming this is in your home directory, the easiest way to do this is to run:

	`echo "PATH=$PATH:~/Vaccinology_Tools/Heptad_homology" >> ~/.bash_profile && source ~/.bash_profile`

### 4. Ensure the heptad_id file is executable

- Check the permissions on the file to make sure you are able to execute it.
	- Easiest way to check is to `cd` to the `~/Vaccinology_Tools/Heptad_homology` directory and run `ls -la`
	- If you need to add the executable flag, just run `chmod +x heptad_id`

## Installation for Linux

#### 1. Clone this git repo

- `git clone https://github.com/michellearanha/Vaccinology_Tools.git`

#### 2. Download and Install EMBOSS

- Download EMBOSS via [ftp](http://emboss.open-bio.org/html/adm/ch01s01.html) or using the [fedora repository](https://src.fedoraproject.org/lookaside/extras/EMBOSS/)
- Install the x11 development libraries: 
	- Debian and Ubuntu:`sudo apt install libx11-dev`
	- Rhel and Centos: `sudo yum install libX11-devel`
	- **_NOTE:_** The next command will complain if you dont install these. 
- Untar the tarball: `tar -xvf EMBOSS-6.6.0.tar.gz`
- In the new EMBOSS directory, run `./configure --prefix=/home/$('whoami')/EMBOSS-6.6.0/`
- Run `make`
- Run `make install`
- Make some coffee, the last two commands may take a bit to install. 

#### 3. Edit your ~/.bashrc

- Need to include this line: `PATH=$PATH:</directory/containing/Vaccinology_Tools/Heptad_homology>` 
	- Assuming this is in your home directory, the easiest way to do this is to run:

		 `echo "PATH=$PATH:~/Vaccinology_Tools/Heptad_homology" >> ~/.bashrc && source ~/.bashrc`

#### 4. Ensure the heptad_id file is executable

- Check the permissions on the file to make sure you are able to execute it.
	- Easiest way to check is to `cd` to the `~/Vaccinology_Tools/Heptad_homology` directory and run `ls -la`
	- If you need to add the executable flag, just run `chmod +x heptad_id`

## Options for heptad_id

`-f [<.fasta>] (seq.fasta) (Input)`

 - File with sequences listed in fasta format

`-r [<.txt>/<.dat>/...] (register.txt) (Input)`

 - File with a list of heptad registers of all sequences (predicted or actual heptad site of the first residue of each sequence

`-c [<type_enum>] (Identity) (Input)`
  
 - Criteria : Identity, Similarity

`-t [<float_percentage>] (45) (Optional)`
  
 - threshold percentage:  sequences that share "criteria" greater than threshold will be saved in a separate folder.

`-E [<dir_path>] (Input)`

 - Directory of EMBOSS binaries

`-A [<type_enum>] (water) (Input)`
  
 - Global alignment - Needleman-Wunsch algorithm : needle

 - Local alignment - Waterman algorithm : water


## Examples 

#### Running a sliding window

`sliding_debug1 [-f <file:fasta file with sequences>] [-w <window size>] [-g <gap size>]  [-t <threshold: percentage>] [-c <String:Criteria-Identity/Similarity>]  [-E <Directory_of_EMBOSS_program>] [-A <alignment_algorithm>] [-l <length of the longest sequence>] [-D <output folder name>] `

#### Running a heptad homology program

`heptad_id [-f <file:fasta file with sequences>][-r <file:heptad register file>][-c <String:Criteria-Identity/Similarity>][-t<threshold: percentage>][-E <Directory_of_EMBOSS_program>] [-A <alignment_algorithm>]`

#### Specific examples:

`sliding_debug1 -f ~/Vaccinology_Tools/Heptad_homology/seq.fasta -w 10 -g 4 -t 25 -c Identity -E ~/EMBOSS-6.6.0/bin/bin -A needle -l 40 -D ~/Vaccinology_Tools/Heptad_homology/NTC6`

`sliding_debug1 -f ~/Vaccinology_Tools/Heptad_homology/seq.fasta -w 15 -g 4 -t 50 -c Similarity -E ~/EMBOSS-6.6.0/bin/bin -A needle -l 40 -D ~/Vaccinology_Tools/Heptad_homology/NTC6`

`heptad_id -f ~/Vaccinology_Tools/Heptad_homology/seq.fasta -r ~/Vaccinology_Tools/Heptad_homology/register.txt -c Identity -t 25 -E ~/EMBOSS-6.6.0/bin -A needle`

`heptad_id -f ~/Vaccinology_Tools/Heptad_homology/seq.fasta -r ~/Vaccinology_Tools/Heptad_homology/register.txt -c Similarity -t 60 -E ~/EMBOSS-6.6.0/bin -A water`


## Definitions and Long Explanations

What the program does:

- Sliding window:

Scanning for matches between sliding k-mer sequence of vaccine type and entire sequence of non-vaccine types (Sliding window approach).

Consider a set of potentially cross reactive protein sequences that are related within species (e.g. antigenic M proteins of different strains of S.pyogenes ) or between species (e.g. antigenic proteins of different viruses belonging to the Flaviviridae virus family such as Dengue, Zika, West Nile that are known induce cross reactive antibodies). Each sequence in the given set is iteratively made the reference sequence and the tool uses a sliding window approach to calculate a pairwise alignment and the number of residues that are identical between successive fragments of the reference sequence and the sequence of the remaining proteins in the set. Matches in terms of sequence similarity as a criterion can also be chosen instead of sequence identity. The subsequence size/ window length/ fragment size (k) for the reference sequence can be adjusted, and the number of amino acids by which the window moves down the reference sequence can be adjusted by adjusting the gap length (m). This gap also determines the degree of overlap between successive refence sequence fragments. The sliding window direction is from the N-terminal to the C-terminal. The number of identical residues between vaccine and non-vaccine type is plotted as a function of window number. If desired, an overall sequence identity or similarity cutoff can be chosen to limit the number of sequences compared. 

- Helical wheel homology:

A standard or canonical coiled coil structure is built by two or more helices twisting around each other forming bundles with their side chains interlocking in a ‘knobs’ into ‘holes’ packing. The regular meshing of knobs into holes packing requires recurrent positions of the side-chains every seven residues along the helix interface. This seven-residue sequence repeat is called a heptad repeat and the positions in the heptad repeat are labelled a-g. The core forming positions (a and d) are usually occupied by hydrophobic residues whereas the remaining, solvent exposed positions (b, c, e, f and g) are dominated by hydrophilic residues. For proteins with coiled-coil secondary structure, machine learned models implemented in webservers exist that can assign the residues to the heptad pattern seen in coiled-coil proteins. We compared the assignment of individual residues in a coiled coil sequence to the heptad by the program MARCOIL (website) and found that it was identical to the actual heptad assignment of three M protein coiled coil crystal structures (pdb ids: 2oto, 5hyt, 5hzp) that were of interest to us. Similar assignments can also be done using NCOILS, PCOILS. 

The tool requires three inputs 1) a text file with the sequences of all vaccine and non-vaccine types in the fasta format, 2) a text file with identifiers/headers of the fasta sequences without ‘>’ character 3) a file with the heptad registers of all sequences (assigned heptad position of the first residue of all sequences) which can be obtained using MARCOIL. The identity at corresponding heptad positions of each peptide/protein with every other peptide/protein in the given set is output in a tabular format. Additionally, an output is also provided as a matrix in tabular format that contains heptad identity of every sequence with every other sequence which can be clustered to obtain clusters of cross-reactive peptides. An R-based script for clustering is also provided. 

