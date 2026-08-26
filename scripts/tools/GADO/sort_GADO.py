import json
import requests
import os
import csv
import time
import re
import glob
import pandas as pd
import logging
import sys

# GET https://www.genenetwork.nl/api/v1/prioritization/HP:0001874,HP:0001419,HP:0002718,HP:0004313,HP:0000951

import ssl


def setup_logger(log_file_path):
    """
    设置logger，同时输出到终端和文件
    """
    # 创建logger
    logger = logging.getLogger('GADO_Logger')
    logger.setLevel(logging.INFO)

    # 清除之前的handlers（避免重复）
    if logger.handlers:
        logger.handlers.clear()

    # 创建文件handler
    file_handler = logging.FileHandler(log_file_path, mode='w', encoding='utf-8')
    file_handler.setLevel(logging.INFO)

    # 创建控制台handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)

    # 设置格式
    formatter = logging.Formatter('%(message)s')
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)

    # 添加handlers
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger


def test_example():
    url = 'https://www.genenetwork.nl/api/v1/prioritization/'
    hpo = 'HP:0001874,HP:0001419,HP:0002718,HP:0004313,HP:0000951'
    URL_HPO = url + hpo
    requests.DEFAULT_RETRIES = 15
    response = requests.get(URL_HPO, verify=False)
    TEXT = response.text
    new_text = eval(TEXT)
    print(new_text)


def call_CADO(HPO):
    url = 'https://www.genenetwork.nl/api/v1/prioritization/'
    # hpo = 'HP:0001874,HP:0001419,HP:0002718,HP:0004313,HP:0000951'
    URL_HPO = url + HPO
    requests.DEFAULT_RETRIES = 25
    response = requests.get(URL_HPO, verify=False)
    TEXT = response.text
    answer = eval(TEXT)
    result = answer['results']
    return result


def deal_CADO(answer, target_gene, select_num=10, threshold=1.0):
    gene_output = ''
    flag = 0
    gene_score_dict = {}
    for dict in answer:
        gene_line = dict['gene']
        gene_name = gene_line['name'].strip()
        gene_score = gene_line['genePredScore']
        gene_score_dict[gene_name] = float(gene_score)
        # if float(gene_score) > 1.0:
        #     gene_output += gene_name + "(" + gene_score + ");"
        #     if gene_name == target_gene:
        #         flag = 1
    return gene_score_dict


def select(gene_score_dict, target_gene, select_num=115, threshold=1.0):
    flag = 0
    # new_dict = gene_score_dict
    new_dict = sorted(gene_score_dict.items(), key=lambda item: item[1], reverse=True)
    # sorted(gene_score_dict.items(), key=lambda item: item[1])
    #
    # select_list = new_dict[0:select_num]
    select_list = new_dict[0:]
    gene_output = ''
    for line in select_list:
        i = line
        gene = line[0]
        score = line[1]
        gene_line = gene + "(" + str(score) + ");"
        gene_output += gene_line
        if gene == target_gene:
            flag = 1

    return gene_output[:-1], flag


def GADO_sort(HPOcsv_filepath):
    title = ["case_ID", "target_gene", "input_HPO", "output_gene([paper,score])", "flag"]
    res_list = [title]
    with open(HPOcsv_filepath, "r", errors='ignore') as csv_file:
        reader = csv.reader(csv_file)
        reader.__next__()
        for line in reader:
            case = int(line[0].strip())
            TARGET_gene = line[1].strip()
            phenotype = line[5].strip()
            phenotype = phenotype.replace(';', ',')
            try:
                answer = call_CADO(phenotype)
                gene_score_dict = deal_CADO(answer, TARGET_gene)
                gene_output, flag = select(gene_score_dict, TARGET_gene, select_num=115)
                write_line = [case, TARGET_gene, phenotype, gene_output, flag]
                res_list.append(write_line)
                print(write_line)
            except:
                gene_output = "Internal Server Error"
                flag = 500
                write_line = [TARGET_gene, phenotype, gene_output, flag]
                res_list.append(write_line)
                print(TARGET_gene + " " + gene_output + "!!!")
    return res_list


def GADO_sort(txt_filepath):
    title = ["case_ID", "target_gene", "input_HPO", "output_gene([paper,score])", "flag"]
    res_list = [title]

    # 1. 从文件名提取 case_ID 和 target_gene
    basename = os.path.basename(txt_filepath)  # 例如 "case-74_gene-IL17F.txt"
    match_case = re.search(r'case-(\d+)', basename)
    match_gene = re.search(r'_gene-([^.]+)', basename)
    if match_case and match_gene:
        case_id = int(match_case.group(1))
        target_gene = match_gene.group(1).strip()
    else:
        raise ValueError(f"文件名 {basename} 必须包含 'case-数字' 和 '_gene-基因名' 格式")

    # 2. 读取文件内容，提取所有 HP: 开头的表型
    with open(txt_filepath, "r", errors='ignore') as f:
        content = f.read().strip()
    tokens = content.split()
    hpo_list = [token for token in tokens if token.startswith('HP:')]
    if not hpo_list:
        raise ValueError(f"文件 {txt_filepath} 中未找到任何 HP: 开头的表型")
    # 原代码将分号替换为逗号，此处直接用逗号连接
    phenotype = ','.join(hpo_list)

    # 3. 调用 CADO 相关函数（与原逻辑完全相同）
    try:
        answer = call_CADO(phenotype)
        gene_score_dict = deal_CADO(answer, target_gene)
        gene_output, flag = select(gene_score_dict, target_gene, select_num=115)
        write_line = [case_id, target_gene, phenotype, gene_output, flag]
        res_list.append(write_line)
        print(write_line)
    except Exception as e:
        # 异常时也输出完整的 5 列，避免列数不一致
        gene_output = "Internal Server Error"
        flag = 500
        write_line = [case_id, target_gene, phenotype, gene_output, flag]
        res_list.append(write_line)
        print(f"{target_gene} {gene_output} !!!")

    return res_list


def GADO_sort_new(txt_filepath, logger=None):
    """
    修改后的GADO_sort_new，支持传入logger记录输出
    """
    title = ["case_ID", "target_gene", "input_HPO", "output_gene([paper,score])", "flag"]
    res_list = [title]

    # 从文件名提取 case_ID 和 target_gene
    basename = os.path.basename(txt_filepath)
    match_case = re.search(r'case-(\d+)', basename)
    match_gene = re.search(r'_gene-([^.]+)', basename)
    if match_case and match_gene:
        case_id = int(match_case.group(1))
        target_gene = match_gene.group(1).strip()
    else:
        error_msg = f"文件名 {basename} 必须包含 'case-数字' 和 '_gene-基因名'"
        if logger:
            logger.error(error_msg)
        raise ValueError(error_msg)

    # 从文件内容提取所有 HP: 开头的表型，并用逗号连接
    with open(txt_filepath, "r", errors='ignore') as f:
        content = f.read().strip()
    tokens = content.split()
    hpo_list = [token for token in tokens if token.startswith('HP:')]
    if not hpo_list:
        error_msg = f"文件 {txt_filepath} 中未找到任何 HP: 开头的内容"
        if logger:
            logger.error(error_msg)
        raise ValueError(error_msg)
    phenotype = ','.join(hpo_list)

    try:
        answer = call_CADO(phenotype)
        gene_score_dict = deal_CADO(answer, target_gene)
        gene_output, flag = select(gene_score_dict, target_gene, select_num=115)
        write_line = [case_id, target_gene, phenotype, gene_output, flag]
        res_list.append(write_line)

        # 使用logger输出
        output_msg = str(write_line)
        if logger:
            logger.info(output_msg)
        else:
            print(output_msg)

    except Exception as e:
        # 异常时也输出 5 列，其中 case_id 和 target_gene 仍可获取
        gene_output = "Internal Server Error"
        flag = 500
        write_line = [case_id, target_gene, phenotype, gene_output, flag]
        res_list.append(write_line)

        error_msg = f"{target_gene} {gene_output} !!!"
        if logger:
            logger.error(error_msg)
        else:
            print(error_msg)

    return res_list


def process_gado(folder_path, output_csv=None, log_file=None):
    """
    遍历文件夹中所有 case-*_gene-*.txt 文件，调用 GADO_sort() 处理，
    合并所有结果并返回列表，可选保存为 CSV。
    新增：记录所有终端输出到log_file（如果提供）
    """
    # 如果提供了log_file路径，设置logger
    logger = None
    if log_file:
        logger = setup_logger(log_file)
        logger.info(f"开始处理文件夹: {folder_path}")
        if output_csv:
            logger.info(f"输出CSV文件: {output_csv}")
        logger.info(f"日志文件: {log_file}")
        logger.info("-" * 50)

    pattern = os.path.join(folder_path, "case-*_gene-*.txt")
    txt_files = glob.glob(pattern)

    all_rows = []
    header = None
    time_sum = 0

    # 记录处理的文件数量
    total_files = len(txt_files)
    if logger:
        logger.info(f"找到 {total_files} 个文件待处理")
        logger.info("-" * 50)

    for idx, file_path in enumerate(txt_files, 1):
        progress_msg = f"Processing: {file_path} ({idx}/{total_files})"
        if logger:
            logger.info(progress_msg)
        else:
            print(progress_msg)

        time_start = time.time()
        res = GADO_sort_new(file_path, logger=logger)
        time_end = time.time()
        time_sim = time_end - time_start

        time_msg = f"Processed: {file_path}; time: {time_sim:.2f}s"
        if logger:
            logger.info(time_msg)
        else:
            print(time_msg)

        time_sum = time_sum + time_sim

        if header is None:
            header = res[0]
            all_rows.append(header)
        all_rows.extend(res[1:])

    final_msg = f"All test file is processed; total time: {time_sum:.2f}s"
    if logger:
        logger.info("-" * 50)
        logger.info(final_msg)
    else:
        print(final_msg)

    if output_csv and all_rows:
        df = pd.DataFrame(all_rows[1:], columns=all_rows[0])
        df.to_csv(output_csv, index=False)
        save_msg = f"Results saved to {output_csv}"
        if logger:
            logger.info(save_msg)
        else:
            print(save_msg)

    return all_rows


def save_csv(data_list, save_dir, save_name="sortBy_GADO-115-oridata.csv"):
    save_path = os.path.join(save_dir, save_name)
    with open(save_path, 'w', newline="") as f:
        writer = csv.writer(f)
        for line in data_list:
            writer.writerow(line)
    print("saved")


def generate_GADO_input(datapath, save_dir):
    with open(datapath, "r", errors='ignore') as csv_file:
        reader = csv.reader(csv_file)
        reader.__next__()
        for line in reader:
            case = line[0].strip()
            TARGET_gene = line[1].strip()
            phenotype = line[5].strip()
            phenotype_new = phenotype.replace(';', '\t')
            sample_name = "case-" + str(case) + "_gene-" + TARGET_gene
            new_line = sample_name + '\t' + phenotype_new

            txtpath = os.path.join(save_dir, sample_name + '.txt')
            with open(txtpath, 'w') as txt_file:
                txt_file.write(new_line)
            print("ok")


def try_again(result_filepath):
    title = ["target_gene", "input_HPO", "output_gene([paper,score])", "flag"]
    new_res_list = [title]
    with open(result_filepath, 'r', errors='ignore') as f:
        reader = csv.reader(f)
        reader.__next__()
        for line in reader:
            TARGET_gene = line[0]
            HPO = line[1]
            ori_flag = line[3]
            if ori_flag == '500':
                try:
                    answer = call_CADO(HPO)
                    gene_score_dict = deal_CADO(answer, TARGET_gene)
                    gene_output, flag = select(gene_score_dict, TARGET_gene, select_num=115)
                    write_line = [TARGET_gene, HPO, gene_output, flag]
                    new_res_list.append(write_line)
                    print(write_line)
                except:
                    gene_output = "Internal Server Error"
                    flag = 500
                    write_line = [TARGET_gene, HPO, gene_output, flag]
                    new_res_list.append(write_line)
                    print(TARGET_gene + " " + gene_output + "!!!")
    return new_res_list


def deal_cmd_res(res_dir):
    title = ["target_gene", "case", "output_gene(rank)", "flag"]
    new_res_list = [title]
    for path, dirs, files in os.walk(res_dir):
        # print(path, dirs, files)
        for file in files:
            if file.endswith('txt'):
                filename = file
                target_gene = filename.split('-')[-1]
                target_gene = target_gene[:-4].strip()
                case = filename.split('_')[0][5:]
                if "sample" in file:
                    continue
                filepath = os.path.join(path, file)
                with open(filepath, 'r', errors='ignore') as file:
                    next(file)
                    gene_line = ''
                    flag = 0
                    for line in file:
                        line_list = line.split('\t')
                        gene = line_list[1].strip()
                        if gene == '':
                            gene = line_list[0]
                        rank = line_list[2]
                        score = line_list[3]
                        # if float(rank) < 11:
                        if float(score) >= 5:
                            if target_gene == gene:
                                flag = 1
                            gene_rank = gene + "(" + rank + ")"
                            gene_line += (gene_rank + ';')
                    write_line = [target_gene, case, gene_line, flag]
                    new_res_list.append(write_line)
    return new_res_list


if __name__ == '__main__':
    # test()
    test_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "test_data", "hpo_files")
    input_dirname = "converted-GRCh37"
    input_dir = os.path.join(test_dir, input_dirname)
    output_filename = input_dirname + "-result-2.csv"
    save_dir = os.path.join(os.path.dirname(__file__), "result")
    save_input_dir = os.path.join(test_dir, "GADO_input43")
    ori_filename = "已发表的诊断案例-HPO.csv"
    # ori_filename = "HPO-from-data_phenotype2.csv"
    # ori_filename = "GADO-ORIGINDATA.csv"
    ori_filepath = os.path.join(test_dir, ori_filename)
    number_file_name = "used_samples_combine.txt"
    number_file_path = os.path.join(test_dir, number_file_name)

    # generate_GADO_input(ori_filepath,save_input_dir)
    # cmd_res_dir = os.path.join(test_dir,"GADO_result_paper")
    # res_list = deal_cmd_res(cmd_res_dir)
    # save_csv(res_list, save_dir, save_name="sortBy_CMD-GADO-10-paper.csv")

    # 修改：生成与output_filename同名的txt日志文件路径
    output_txt_filename = input_dirname + "-result-2.txt"
    output_log_path = os.path.join(save_dir, output_txt_filename)

    # 确保保存目录存在
    os.makedirs(save_dir, exist_ok=True)

    # 调用process_gado，传入日志文件路径
    output_csv_path = os.path.join(save_dir, output_filename)
    process_gado(input_dir, output_csv=output_csv_path, log_file=output_log_path)

    # res_list = GADO_sort(ori_filepath)

    # save_csv(res_list,save_dir,save_name= "sortBy_GADO-115-oridata2.csv")
    print("ok")

    # target_gene = "EIF2AK3"
    # HPO = "HP:0100651,HP:0001953"
    # answer = call_CADO(HPO)
    # GENE,flag = deal_CADO(answer,target_gene)
    # save_csvfilepath = os.path.join(save_dir,'GADO-phennotype2-2.csv')
    # new_res = try_again(save_csvfilepath)
    # save_csv(new_res,save_dir,'GADO-phennotype2-3.csv')