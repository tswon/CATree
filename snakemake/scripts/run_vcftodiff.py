import os
import argparse
import gzip
import logging
import subprocess

#this script requires individual VCFs
parser = argparse.ArgumentParser()
parser.add_argument('-vd', '--VCF_directory', required=True, type=str,help='path to directory of single-sample VCFs')
parser.add_argument('-wd', '--working_directory', required=True, type=str, help='directory for all outputs (make sure this directory will have enough space!!!!)')
parser.add_argument('-bd', '--bedgraph_directory', required=True, type=str, help="path to directory of bed coverage file (bedgraph) for vcf")
parser.add_argument('-sl', '--SRA_list_file', required=True, type=str,help='file with the list of SRA want to be processed')

args = parser.parse_args()
vd = args.VCF_directory
wd = args.working_directory
# wd = "diff_files"
bd = args.bedgraph_directory
sl = args.SRA_list_file

with open(sl, 'r') as SRA_list:
    for sra in SRA_list:
        sra = sra.strip()
        
        vcfs = [f"{sra}.vcf.gz"]
        beds = [f"{sra}_merged.bed"]
        diffs = [f"{sra}.diff"]
        
        for vcf, bed, diff in zip(vcfs, beds, diffs):
            diff_path = os.path.join(wd, diff)
            bed_path = os.path.join(bd, bed)
            vcf_path = os.path.join(vd, vcf)
            
            # Debugging statements
            print(f"Current working directory: {os.getcwd()}")
            print(f"Checking for existence of {diff_path}")
            print(vcf, bed)
        
            # arg = f'python scripts/vcf_to_diff_script.py -v {os.path.join(vd, vcf)} -d {wd} -bed {os.path.join(bd, bed)}'
            # print('arg', arg)    
            
            # check if the .diff file already exists
            if os.path.exists(diff_path) and diff_path.endswith(".diff"):
                print(f"Skipping {vcf}: Output file {diff_path} already exists.")
                continue
            
            elif not os.path.exists(bed_path):
                print(f"Skipping {vcf}: Bedgraph file {bed_path} does not exists.")
                continue
            
            elif not os.path.exists(vcf_path):
                print(f"Skipping {vcf}: vcf file {vcf_path} does not exists.")
                continue                

            else:
                arg = f'python scripts/vcf_to_diff_script.py -v {os.path.join(vd, vcf)} -d {wd} -bed {os.path.join(bd, bed)}'
                subprocess.run(arg, shell=True, check=True)
                print('arg', arg)
                print("Finished")
        

# # get a list of files in vd with .vcf.gz extension
# vcfs = [f for f in os.listdir(vd) if f.endswith(".vcf.gz")]

# # create a list of corresponding BED files by changing the extension
# beds = [os.path.splitext(os.path.splitext(f)[0])[0] + "_merged.bed" for f in vcfs]

# # create a list of corresponding DIFF files
# diffs = [os.path.splitext(os.path.splitext(f)[0])[0] + ".vc.diff" for f in vcfs]


# for vcf, bed, diff in zip(vcfs, beds, diffs):
#     diff_path = os.path.join(wd, diff)
#     bed_path = os.path.join(bd, bed)
    
#     # Debugging statements
#     print(f"Current working directory: {os.getcwd()}")
#     print(f"Checking for existence of {diff_path}")
#     print(vcf, bed)
 
#     # arg = f'python scripts/vcf_to_diff_script.py -v {os.path.join(vd, vcf)} -d {wd} -bed {os.path.join(bd, bed)}'
#     # print('arg', arg)    
    
#     # check if the .diff file already exists
#     if os.path.exists(diff_path) and diff_path.endswith(".vc.diff"):
#         print(f"Skipping {vcf}: Output file {diff_path} already exists.")
#         continue
    
#     elif not os.path.exists(bed_path):
#         print(f"Skipping {vcf}: Bedgraph file {bed_path} does not exists.")
#         continue

#     else:
#         arg = f'python scripts/vcf_to_diff_script.py -v {os.path.join(vd, vcf)} -d {wd} -bed {os.path.join(bd, bed)}'
#         subprocess.run(arg, shell=True, check=True)
#         print('arg', arg)
#         print("Finished")

    
    # try:
    #     subprocess.run(arg, shell=True, check=True)
    #     print('arg', arg)
    #     print("Finished")
    # except subprocess.CalledProcessError as e:
    #     # If subprocess run fails, print the error, delete the specific output file, and move on
    #     print(f"Error: {e}")
    #     print(f"Deleting {diff_path}")
    #     os.remove(diff_path)
        