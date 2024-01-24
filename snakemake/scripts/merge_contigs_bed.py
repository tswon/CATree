# merge_contigs_bed.py
"""
A program designed to modify the contig positions in a bed file 
to start from the end of the previous contig.
"""

import argparse

#this script requires bed files
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, type=str,help='input bed file to convert')
parser.add_argument('-o', '--output', required=True, type=str, help='output bed file with merged contigs')


args = parser.parse_args()
input_file = args.input
output_file = args.output


def merge_contigs(input, output):
    new_start_pos_2 = 3148135
    new_start_pos_3 = new_start_pos_2 + 2554418
    new_start_pos_4 = new_start_pos_3 + 2336890
    new_start_pos_5 = new_start_pos_4 + 1318327
    new_start_pos_6 = new_start_pos_5 + 1007026
    new_start_pos_7 = new_start_pos_6 + 1004684
        
    with open(input, 'r') as infile, open(output, 'w') as outfile:
        for line in infile:
            columns = line.strip().split('\t')

            contig = columns[0]
            start_pos = int(columns[1])
            end_pos = int(columns[2])
            coverage = columns[3]

            if contig == 'NC_072812.1':
                contig = 'NC_07281X.1'
        
            elif contig == 'NC_072813.1':
                contig = 'NC_07281X.1'
                start_pos += new_start_pos_2
                end_pos += new_start_pos_2
                
            elif contig == 'NC_072814.1':
                contig = 'NC_07281X.1'
                start_pos += new_start_pos_3
                end_pos += new_start_pos_3
            
            elif contig == 'NC_072815.1':
                contig = 'NC_07281X.1'
                start_pos += new_start_pos_4
                end_pos += new_start_pos_4
                
            elif contig == 'NC_072816.1':
                contig = 'NC_07281X.1'
                start_pos += new_start_pos_5
                end_pos += new_start_pos_5
                
            elif contig == 'NC_072817.1':
                contig = 'NC_07281X.1'
                start_pos += new_start_pos_6
                end_pos += new_start_pos_6
            
            elif contig == 'NC_072818.1':
                contig = 'NC_07281X.1'
                start_pos += new_start_pos_7
                end_pos += new_start_pos_7
                
            outfile.write(f"{contig}\t{start_pos}\t{end_pos}\t{coverage}\n")


if __name__ == "__main__":
              
    merge_contigs(input_file, output_file)