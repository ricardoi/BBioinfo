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

