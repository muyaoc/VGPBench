import os
import subprocess
import time
import numpy as np

my_dir = os.path.dirname(__file__)
input_dir = os.path.join(my_dir,"converted-hpo_files_2021_clean")
input_hpo_dir = os.path.join(input_dir,"converted_txts")#//os.path.join(input_dir,"HPO")
input_gene_dir = os.path.join(input_dir,"gene")
res_dir = os.path.join(my_dir,"result-hpo_files_2021_clean")

# python3 phen2gene.py -f example/HPO_sample.txt -v -out out/prioritizedgenelist
# python3 phen2gene.py -f example/case_1-gene_ABCC9.txt -out out/ABCC9 -l example/TRY.txt
# case_1-gene_ABCC9.txt:HPO  ; TRY.txt:gene

def find_gene_txt(gene):
      for path, dirs, files in os.walk(input_gene_dir):
        # print(path, dirs, files)
        for file in files:
            filename = file
            find_gene = filename.split("-")[1][5:]
            if gene == find_gene:
                gene_filepath = os.path.join(path,file)
                return gene_filepath

           

def sort_hpo_gene(HPO_filename):
     HPO_filepath = os.path.join(input_hpo_dir,HPO_filename)
    #  gene = HPO_filename.split("-")[1][5:-4]
     gene = HPO_filename.split("_")[1][5:-4]

     gene_filepath = find_gene_txt(gene)
     dir_name = HPO_filename[:-4]
     out_dir = os.path.join(res_dir,"HPO+GENE",dir_name)
    #  out_dir = os.path.join(res_dir,gene)
     if not os.path.exists(out_dir):
        os.makedirs(out_dir)
     cmd_line = ['python3', 'phen2gene.py', '-f', HPO_filepath, '-out', out_dir, '-l', gene_filepath]
     result = subprocess.run(cmd_line, stdout=subprocess.PIPE, text=True)
     print(result.stdout)

def sort_HPO(HPO_filename):
    HPO_filepath = os.path.join(input_hpo_dir,HPO_filename)
    # gene = HPO_filename.split("-")[1][5:-4]
    gene = HPO_filename.split("_")[1][5:-4]

    # case = HPO_filename.split("-")[0][5:]
    dir_name = HPO_filename[:-4]
    out_dir = os.path.join(res_dir,"HPO",dir_name)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    cmd = 'python3 phen2gene.py -f example/HPO_sample.txt -v -out out/prioritizedgenelist'
    cmd_line = ['python3', 'phen2gene.py', '-f', HPO_filepath, '-v', '-out', out_dir]
    result = subprocess.run(cmd_line, stdout=subprocess.PIPE, text=True)
    print(result.stdout)

def walklog(input_HPO_dir):
    for path, dirs, files in os.walk(input_HPO_dir):
        # print(path, dirs, files)
        for file in files:
            filename = file
            # sort_hpo_gene(filename)
            sort_HPO(filename)
            print("-------------" + filename + " is done-------------")

def main():
    test_HPOname = "case_1-gene_ABCC9.txt"
    test_HPOpath = os.path.join(input_hpo_dir,test_HPOname)
    # sort(test_HPOpath)
    # sort_HPO(test_HPOpath)
    time_start = time.time()
    walklog(input_hpo_dir)
    time_end = time.time() 
    time_sum = time_end - time_start
    print("筛选总共用时：", time_sum)


main()