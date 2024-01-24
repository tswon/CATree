import os
import argparse
import subprocess

#this script requires individual VCFs
parser = argparse.ArgumentParser()
parser.add_argument('-bd', '--BED_directory', required=True, type=str,help='path to directory of bed files')
parser.add_argument('-wd', '--working_directory', required=True, type=str,help='path to directory of merged bed files')
parser.add_argument('-sl', '--SRA_list_file', required=True, type=str,help='file with the list of SRA want to be processed')
args = parser.parse_args()
bd = args.BED_directory
wd = args.working_directory
sl = args.SRA_list_file

# make output directory if doesn't exist
os.makedirs(wd, exist_ok=True)

with open(sl, 'r') as SRA_list:
    for sra in SRA_list:
        sra = sra.strip()
        bed_filename = f"aligned_{sra}.bed"
        bed_path = os.path.join(bd, bed_filename)
        
        # bed_files = [f for f in os.listdir(bd) if f.endswith(".bed")]
        
        if os.path.exists(bed_path):
            # for bed in bed_filename:
            # input_prefix = bed_filename.split("aligned_")[-1]
            bm = f"{sra}_merged.bed"
            
            # current_bed_path = os.path.join(bd, bed)
            output_path = os.path.join(wd, bm)
            
            # check if the file already exists
            if os.path.exists(output_path):
                print(f"Skipping {sra}: Output file {output_path} already exists.")
                continue        
            
            else:
                arg = f'python scripts/merge_contigs_bed.py -i {bed_path} -o {wd}/{bm}'
                print('arg', arg)
                
                try:
                    subprocess.run(arg, shell=True, check=True)
                except subprocess.CalledProcessError as e:
                    # If subprocess run fails, print the error, delete the specific output file, and move on
                    print(f"Error: {e}")
                    print(f"Deleting {output_path}")
                    os.remove(output_path)
        else:
            print(f"Bed file for SRA {sra} not found in {bd}")
    


# for bed in bed_files:
#     input_prefix = bed.split("aligned_")[-1]
#     bm = input_prefix.replace('.bed', '') + "_merged.bed"
    
#     bed_path = os.path.join(bd, bed)
#     output_path = os.path.join(wd, bm)
    
#     # check if the file already exists
#     if os.path.exists(output_path):
#         print(f"Skipping {bed}: Output file {output_path} already exists.")
#         continue        
    
#     else:
#         arg = f'python scripts/merge_contigs_bed.py -i {bed_path} -o {wd}/{bm}'
#         print('arg', arg)
        
#         try:
#             subprocess.run(arg, shell=True, check=True)
#         except subprocess.CalledProcessError as e:
#             # If subprocess run fails, print the error, delete the specific output file, and move on
#             print(f"Error: {e}")
#             print(f"Deleting {output_path}")
#             os.remove(output_path)