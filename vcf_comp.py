import argparse
import logging
import os
import time
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description='Compare two VCF files. Find genotype changes.')


parser.add_argument('-1', '--vcf_1', type=str,
                    help='vcf_1.', required=True)
parser.add_argument('-2', '--vcf_2', type=str,
                    help='vcf_1.', required=True)
parser.add_argument('--only_merge', action='store_true',
                    help='only_merge flag.')
parser.add_argument('-w', '--work_dir', type=str,
                    help='Working dir absolute path.', required=True)

def change_mkdir(out_dir_path):
    if not os.path.exists(out_dir_path):
        os.mkdir(out_dir_path)
    os.chdir(out_dir_path)

def read_vcf(vcf_file_name):
    # #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  smpl.bam
    with open(vcf_file_name, "r") as vcf_file:
        for line in vcf_file:
            if line.startswith("#CHROM"):
                vcf_names = [x for x in line.strip().split('\t')]
                break
    vcf_names[-1] = 'GT:PL'
    vcf_df = pd.read_csv(vcf_file_name, comment='#', delim_whitespace=True, header=None, names=vcf_names)
    return vcf_df

def read_vcf_txt(vcf_file_name):
    vcf_df = pd.read_csv(vcf_file_name, sep='\t')
    return vcf_df

def merge_two_vcf(vcf_1, vcf_2):
    columns_on = ['#CHROM', 'POS']
    
    vcf_1.sort_values(by=columns_on, inplace=True)
    vcf_2.sort_values(by=columns_on, inplace=True)
    merged_vcf = pd.merge(vcf_1, vcf_2, how='inner', on=columns_on)
    return merged_vcf

def proc_merged_vcf(merged_vcf, min_homozig_cov = 10, min_MQ = 40):
    diff_gt_merged_vcf = merged_vcf[merged_vcf['ALT_x'] != merged_vcf['ALT_y']].copy(deep = True)
    passed_variants = []
    VAF_x_col = []
    VAF_y_col = []
    for info_x, info_y, gt_x, gt_y in zip(diff_gt_merged_vcf['INFO_x'], diff_gt_merged_vcf['INFO_y'],
                                            diff_gt_merged_vcf['GT:PL_x'], diff_gt_merged_vcf['GT:PL_y']):
        info_x = info_x.strip().split(';')
        info_y = info_y.strip().split(';')
        indel_flag = False
        if 'INDEL' in info_x or 'INDEL' in info_y:
            passed_variants.append(True)
            info_x = { info.split('=')[0] : info.split('=')[1] for info in info_x[1:] }
            info_y = { info.split('=')[0] : info.split('=')[1] for info in info_y[1:] }
            indel_flag = True
        else:
            info_x = { info.split('=')[0] : info.split('=')[1] for info in info_x }
            info_y = { info.split('=')[0] : info.split('=')[1] for info in info_y }

        # Create new cols with VAF
        info_x_DP4 = list(map(float, info_x['DP4'].split(',')))
        info_x_DP4_sum = sum(info_x_DP4)
        VAF_x = 0.0
        if info_x_DP4_sum > 0:
            VAF_x = (info_x_DP4[2] + info_x_DP4[3])/info_x_DP4_sum
        VAF_x_col.append(VAF_x)
        info_y_DP4 = list(map(float, info_y['DP4'].split(',')))
        info_y_DP4_sum = sum(info_y_DP4)
        VAF_y = 0.0
        if info_y_DP4_sum > 0:
            VAF_y = (info_y_DP4[2] + info_y_DP4[3])/info_y_DP4_sum
        VAF_y_col.append(VAF_y)

        if indel_flag:
            continue

        gt_x_0_1 = gt_x.split('/')
        gt_y_0_1 = gt_y.split('/')
        gt_x_0_1[1] = gt_x_0_1[1].split(':')[0]        
        gt_y_0_1[1] = gt_y_0_1[1].split(':')[0]

        info_x['DP'] = float(info_x['DP'])
        info_y['DP'] = float(info_y['DP'])
        info_x['MQ'] = float(info_x['MQ'])
        info_y['MQ'] = float(info_y['MQ'])
        if (gt_x_0_1[0] == gt_x_0_1[1] and info_x['DP'] < min_homozig_cov) or \
            (gt_y_0_1[0] == gt_y_0_1[1] and info_y['DP'] < min_homozig_cov):
            passed_variants.append(False)
        elif info_x['MQ'] < min_MQ or info_y['MQ'] < min_MQ:
            passed_variants.append(False)
        else:
            passed_variants.append(True)

    diff_gt_merged_vcf['VAF_x'] = VAF_x_col
    diff_gt_merged_vcf['VAF_y'] = VAF_y_col
    diff_gt_merged_vcf = diff_gt_merged_vcf[passed_variants]
    return diff_gt_merged_vcf
        

def main():
    args = parser.parse_args()  
    change_mkdir(args.work_dir)
    logging.basicConfig(filename=args.work_dir + '/vcf_comp.log', level=logging.INFO)
    name_vcf_file_1 = os.path.splitext(os.path.basename(args.vcf_1))[0]
    name_vcf_file_2 = os.path.splitext(os.path.basename(args.vcf_2))[0]
    vcf_df_1 = None
    vcf_df_2 = None
    if args.only_merge:
        vcf_df_1 = read_vcf_txt(args.vcf_1)
        vcf_df_2 = read_vcf_txt(args.vcf_2)
    else:
        vcf_df_1 = read_vcf(args.vcf_1)
        vcf_df_2 = read_vcf(args.vcf_2)
    merged_vcf = merge_two_vcf(vcf_df_1, vcf_df_2)
    if args.only_merge:
        merged_vcf.to_csv(args.work_dir + '/' + '_'.join([name_vcf_file_1, name_vcf_file_2]) + '.txt', sep='\t', index = False)
    else:
        diff_gt_merged_vcf = proc_merged_vcf(merged_vcf)
        diff_gt_merged_vcf.to_csv(args.work_dir + '/' + '_'.join([name_vcf_file_1, name_vcf_file_2]) + '.txt', sep='\t', index = False)


if __name__ == "__main__":
    main()