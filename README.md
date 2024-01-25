# CATree - Candida Auris Tree
## Phylogenetic tree building pipeline for _Candida auris_

Rotation project in Dr. Russell Corbett-Detig's lab for 2023-2024 academic year

Pipeline to generate a Phylogenetic tree for Candida auris.

Uses Snakemake, python scripts, and UShER or MAPLE

# Workflow
1. Snakefile
   	-
	- Downloads sequence using [fasterq](https://github.com/ncbi/sra-tools)
	- Aligns C. auris sequence to the reference genome using [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)
		- [Reference seq](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_003013715.1/)
	- Calls variants using GATK's [HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller)
 	- Only accepts SRA files with 2 fasta files (01/24/2024) 

2. Python scripts
	- 	
	- Because UShER doesn't take input with multiple chromosomes, the chromosome and position information need to be merged as one big chromosome
 	- '**merge_contigs_bed.py**' re-writes the **bed** files to have one big chromosome named 'NC_072814.1' with a length combining that of all seven chromosomes
  		- use '**run_mergebed.py**' to run multiple samples at once in the command line
    - '**vcf_to_diff_script.py**' converts **vcf** to **diff**. The program was written by [Lily Karim](https://github.com/lilymaryam/parsevcf) and was modified to make changes to **vcf** before the conversion
    	- This step is required for MAPLE and recommended for UShER to save storage space
    	- '**merge_contigs_vcf.py**' re-writes the **vcf** files to have one big chromosome named 'NC_072814.1' with a length combining that of all seven chromosomes
     	- use '**run_vcftodiff.py**' to run multiple samples at once in the command line

3. Tree Building
   -
   	- Concatenate **diff** files from multiple sample into a single **diff** file
   	- Use a blank tree as an initial tree
   		- Add `(ref);` and save as a Newick tree like `tree.nwk` 
	- '**UShER**' uses maximum parsimony to place samples [UShER wiki](https://usher-wiki.readthedocs.io/en/latest/index.html)
 		- Outputs: .pb and .nh 
   		- `usher-sampled -t PATH/initial_tree.nwk --diff PATH/combined_diff.diff --ref PATH/reference.fasta -o PATH/output_tree.pb -d PATH`
   	 	- Once tree is made, convert to .jsonl format using `usher_to_taxonium -i PATH/tree.pb -o PATH/tree.jsonl` 
   	- '**MAPLE**' is a Likelihood-based phylogenetic analysis tool that places the sample at the node with the highest score [MAPLE github](https://github.com/NicolaDM/MAPLE)
   		- Output: .nh
	   	- Need pypy3
   	 	- `pypy3 scripts/MAPLEv0.3.6.py --input  PATH/combined_diff.diff --reference PATH/reference.fasta --output OUTPUT/PATH`

4. Tree visualization
   -
  	- Upload the tree file (.pb, .nh, or .jsonl) to Taxonium to visualize the tree [taxonium](https://taxonium.org/)
   	- Submit meta data to use filters
   	- Run treenome viewer by using `usher_to_taxonium` command with annotation file (Genome GBFF)
   		- `usher_to_taxonium -i PATH/tree.pb -o PATH/tree_chromosome_1.jsonl.gz --genbank PATH/genome_chromosom_1.gbff`
   	 	- if multiple chromosome present, need to break the gbff into the different chromosome and run command for each chromosome 
