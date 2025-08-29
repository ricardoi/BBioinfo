# BBioinfo
Begining installations for mac users to search for viruses in HTS data

## Setting up the terminal
Open the terminal (open the spotlight and search `terminal`) and Install xcode by pasting this command 
```bash
xcode-select --install
```
Check your installation with `xcode-select -p`, if you get this line `/Applications/Xcode.app/Contents/Developer` your installation was succesfull.

Now, we need to install `Homebrew`, a package manager to install software in Mac, for linux users, this is similar to `apt get`.  
```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
Test if brew was installed succesfully with this command `brew install wget`. 

If everything is working, now you can start installing software related to bioinformatics quite easily.

## Begining Software Installation for Bionformatics
Begin installing [Spades](http://ablab.github.io/spades/), a sequece assembler great for viruses, prokaryots and small eukaryots. 
Install using the following command.
```bash
brew install spades
```
If installation was succesful, try typing `spades.py` in the terminal.
You will see the following information 
```bash
SPAdes genome assembler v4.2.0

Usage: spades.py [options] -o <output_dir>

Basic options:
  -o <output_dir>             directory to store all the resulting files (required)
  --isolate                   this flag is highly recommended for high-coverage isolate and multi-cell data
  --sc                        this flag is required for MDA (single-cell) data
  --meta                      this flag is required for metagenomic data
 ...
  -h, --help                  prints this usage message
  -v, --version               prints version
```

Now you can try the next command with your fastq files. Plug in your computer, it will drain your battery otherwise.
```bash
spades.py --rnaviral -1 file1_R1.fq.gz -2 file1_R2.fq.gz -t 6  -o results_out
```
The flag `rnaviral` will enables virus assembly module from RNA-Seq data. The `-1` and `-2` indicates the input data for paired ends reads, the `-t` indicates the number of processors, and `-o` indicates the output directory to save the results. In this folder, you will be have the results for the assembly.

--
Copy the scaffold.fasta files into a new `spades_results` folder. Subset the `fasta` file using `subsetfastas.bash` targeting the expected genome size range infecting your host, e.g. 0.6 to 16 kb for mycoviruses.
You need to clone this repository as follows
```bash
git clone https://github.com/ricardoi/BBioinfo/
cd BBioinfo/scripts/
chmod +x subsetfastas.bash
```
Now you can execute `subsefastas.bash` with `./PATHtoBBioinfo/scripts/subsetfasta.bash 600 16000 scaffolds.fasta` and it will create an output `scaffolds_subset.fasta`.
Type the following:
```bash
./subsetfastas.bash
Usage: subsetfastas.bash MIN MAX file1.fasta [file2.fasta ...]
Example: subsetfastas.bash 800 20000 sample1.fasta sample2.fasta
.fasta, .fa, .fas, .faa, .fna, are accepted extensions
```

## Constructing a mycovirus database using `esearch` from NCBI
Create a database with all queries assigned to a family within the mycovirus group.
```bash
esearch -db nuccore -query \
'(Mitoviridae[Organism] OR Narnaviridae[Organism] OR Partitiviridae[Organism] OR Chrysoviridae[Organism] OR Totiviridae[Organism] OR Polymycoviridae[Organism] OR Fusariviridae[Organism] OR Hypoviridae[Organism] OR Endornaviridae[Organism] OR Quadriviridae[Organism] OR Yadokariviridae[Organism] OR Botourmiaviridae[Organism] OR Mymonaviridae[ORGANISM] OR Megabirnaviridae[Organism] OR Alternaviridae[Organism] AND (complete genome[Title]))' \
| efetch -format fasta > mycoviruses_by_family.fasta
```
Create a database with all complete genomes assigned to a family within the mycovirus group
```bash
esearch -db nuccore -query \
'(Mitoviridae[Organism] OR Narnaviridae[Organism] OR Partitiviridae[Organism] OR Chrysoviridae[Organism] OR Totiviridae[Organism] OR Polymycoviridae[Organism] OR Fusariviridae[Organism] OR Hypoviridae[Organism] OR Endornaviridae[Organism] OR Quadriviridae[Organism] OR Yadokariviridae[Organism] OR Botourmiaviridae[Organism] OR Mymonaviridae[ORGANISM] OR Megabirnaviridae[Organism] OR Alternaviridae[Organism] AND (complete genome[Title]))' \
| efetch -format fasta > mycoviruses_by_family_completegenome.fasta
```
Create a database with all protein queries assigned to a family within the mycovirus group.
```bash
esearch -db protein -query \
'(Mitoviridae[Organism] OR Narnaviridae[Organism] OR Partitiviridae[Organism] OR Chrysoviridae[Organism] OR Totiviridae[Organism] OR Polymycoviridae[Organism] OR Fusariviridae[Organism] OR Hypoviridae[Organism] OR Endornaviridae[Organism] OR Quadriviridae[Organism] OR Yadokariviridae[Organism] OR Botourmiaviridae[Organism] OR Mymonaviridae[ORGANISM] OR Megabirnaviridae[Organism] OR Alternaviridae[Organism])' \
| efetch -format fasta > mycoviruses_by_family_protein.fasta
```
Create a database with the protein queries from the complete genomes assigned to a family within the mycovirus group.
```bash
esearch -db nuccore -query \
'((Mitoviridae[Organism] OR Narnaviridae[Organism] OR Partitiviridae[Organism] OR Chrysoviridae[Organism] OR Totiviridae[Organism] OR Polymycoviridae[Organism] OR Fusariviridae[Organism] OR Hypoviridae[Organism] OR Endornaviridae[Organism] OR Quadriviridae[Organism] OR Yadokariviridae[Organism] OR Botourmiaviridae[Organism] OR Mymonaviridae[Organism] OR Megabirnaviridae[Organism] OR Alternaviridae[Organism])) 
 AND "complete genome"[Title]' \
| elink -target protein \
| efetch -format fasta > mycovirus_proteins_from_complete_genomes.faa
```


Create the databases for complete genome, proteins from complete genomes, all nucleotide sequences and proteins
```bash
makeblastdb -in mycoviruses_by_family.fasta -dbtype nucl -out virusDBnc -parse_seqids -hash_index -title "Mycoviruses (nucl)"

makeblastdb -in mycoviruses_by_family_protein.fasta -dbtype prot -out virusDBaa -parse_seqids -hash_index -title "Mycoviruses (aa)"

makeblastdb -in mycoviruses_by_family-completegenome.fasta -dbtype nucl -out virusDBgeno -parse_seqids -hash_index -title "Mycoviruses (geno)"

makeblastdb -in mycovirus_proteins_from_complete_genomes.faa -dbtype prot -out virusDBprot -parse_seqids -hash_index -title "Mycoviruses (prot)"
```

BLAST
```bash
 blastn -query file_scaffolds.fasta -db virusDBnc -task dc-megablast -evalue 1e-5 -max_target_seqs 5 -num_threads 6 -outfmt '6 qseqid sseqid pident length qcovs evalue bitscore staxids sscinames scomnames sskingdoms stitle' -out file_blastn_res.tsv"
```










