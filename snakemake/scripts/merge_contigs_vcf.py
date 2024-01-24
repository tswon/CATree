# merge_contigs_vcf.py
"""
A program designed to modify the contig positions in a VCF file 
to start from the end of the previous contig.
"""

from cyvcf2 import VCF, Writer
import vcfpy
import os
import argparse

#this script requires individual VCFs
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, type=str,help='input vcf file to convert')
parser.add_argument('-o', '--output', required=False, type=str, default=None, help='output vcf file with merged contigs')


args = parser.parse_args()
input_file = args.input

if args.output is None:
    output_file = input_file.replace('.vcf', '_merged.vcf')
else:
    output_file = args.output


def merge_contigs(input, output):
    vcf = VCF(input)
    w = Writer(output, vcf)
      
    # # Modify the header and add new contig info
    new_contig_line = '##contig=<ID=NC_07281x.1,length=12249773>'
    w.add_to_header(new_contig_line)
              
    
    # new start positions for different cotings
    # set_pos set the POS to the given 1-based position    
    new_start_pos_2 = 3148135 - 1
    new_start_pos_3 = new_start_pos_2 + 2554418
    new_start_pos_4 = new_start_pos_3 + 2336890
    new_start_pos_5 = new_start_pos_4 + 1318327
    new_start_pos_6 = new_start_pos_5 + 1007026
    new_start_pos_7 = new_start_pos_6 + 1004684
    
    
    # make modifications to the start positions for each contig
    for variant in vcf:
        if variant.CHROM.startswith('NC_072812.1'):
            variant.CHROM = variant.CHROM.replace('NC_072812.1', 'NC_07281X.1')
            w.write_record(variant)
            
        elif variant.CHROM.startswith('NC_072813.1'):
            new_pos_2 = variant.POS + new_start_pos_2
            variant.set_pos(new_pos_2) 
            variant.CHROM = variant.CHROM.replace('NC_072813.1', 'NC_07281X.1')
            w.write_record(variant)
            
        elif variant.CHROM.startswith('NC_072814.1'):
            new_pos_3 = variant.POS + new_start_pos_3
            variant.set_pos(new_pos_3) 
            variant.CHROM = variant.CHROM.replace('NC_072814.1', 'NC_07281X.1')
            w.write_record(variant)
            
        elif variant.CHROM.startswith('NC_072815.1'):
            new_pos_4 = variant.POS + new_start_pos_4
            variant.set_pos(new_pos_4)
            variant.CHROM = variant.CHROM.replace('NC_072815.1', 'NC_07281X.1') 
            w.write_record(variant)
            
        elif variant.CHROM.startswith('NC_072816.1'):
            new_pos_5 = variant.POS + new_start_pos_5
            variant.set_pos(new_pos_5)
            variant.CHROM = variant.CHROM.replace('NC_072816.1', 'NC_07281X.1')  
            w.write_record(variant)
            
        elif variant.CHROM.startswith('NC_072817.1'):
            new_pos_6 = variant.POS + new_start_pos_6
            variant.set_pos(new_pos_6)
            variant.CHROM = variant.CHROM.replace('NC_072817.1', 'NC_07281X.1')  
            w.write_record(variant)
            
        elif variant.CHROM.startswith('NC_072818.1'):
            new_pos_7 = variant.POS + new_start_pos_7
            variant.set_pos(new_pos_7)
            variant.CHROM = variant.CHROM.replace('NC_072818.1', 'NC_07281X.1')  
            w.write_record(variant)

      
    # close the vcf files    
    w.close()    
    vcf.close()
        
        
def remove_contig_lines(input, output):
    # Create a reader for the input VCF
    reader = vcfpy.Reader.from_path(input)

    # Create a new header without '##contig' lines
    new_header_lines = [line for line in reader.header.lines if not line.key.startswith('contig')]
    new_header = vcfpy.Header(new_header_lines, reader.header.samples)

    # Create a writer for the output VCF with the modified header
    writer = vcfpy.Writer.from_path(output, header=new_header)    

    # Iterate through the records and write them to the output VCF
    for record in reader:
        writer.write_record(record)

    # Close the reader and writer
    reader.close()
    writer.close()


def add_contig(input, output):
    vcf = VCF(input)
    w = Writer(output, vcf)
      
    # # Modify the header and add new contig info
    new_contig_line = '##contig=<ID=NC_07281x.1,length=12249773>'
    w.add_to_header(new_contig_line)   
    
    for variant in vcf:
        w.write_record(variant)
        
    # close the vcf files    
    w.close()    
    vcf.close() 

if __name__ == "__main__":
    temp_file1 = input_file + 'temp1.vcf'
    temp_file2 = input_file + 'temp2.vcf'
                
    merge_contigs(input_file, temp_file1)    
    remove_contig_lines(temp_file1, temp_file2)
    add_contig(temp_file2, output_file)
    
    
    try:
        # Remove temporary files (if any)
        os.remove(temp_file1)
        os.remove(temp_file2)
    except OSError as e:
        print(f"Error removing temporary files: {e}") 


