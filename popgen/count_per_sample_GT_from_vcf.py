import gzip
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="统计VCF文件的杂合度信息 (原生Python高效版)")
    parser.add_argument('-i', '--input', required=True, help="输入VCF文件路径 (支持 .vcf 和 .vcf.gz)")
    parser.add_argument('-o', '--output', required=True, help="输出统计结果文件路径")
    return parser.parse_args()

def process_vcf(input_file, output_file):
    # 根据后缀判断是否为压缩文件
    if input_file.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt' # text mode
    else:
        open_func = open
        mode = 'r'

    samples = []
    # stats_matrix: 每个元素是一个列表 [Het, HomRef, HomAlt, Missing]
    stats_matrix = []
    
    print(f"正在处理文件: {input_file} ...")

    with open_func(input_file, mode) as f:
        for line in f:
            # 跳过元数据行
            if line[:2] == '##':
                continue
            
            # 解析表头行获取样本名
            if line[:2] == '#C':
                parts = line.strip().split('\t')
                samples = parts[9:]
                print(f"检测到 {len(samples)} 个样本.")
                print(samples)
                # 初始化统计矩阵: [Het, HomRef, HomAlt, Missing]
                stats_matrix = [[0, 0, 0, 0] for _ in range(len(samples))]
                print("初始化完成，开始统计...")
                continue

            parts = line.strip().split('\t')
            
            # 假定GT在FORMAT字段的第一个位置
            # 遍历样本列 (从索引9开始)
            # 使用 enumerate 获取 stats_matrix 对应的索引
            for idx, sample_data in enumerate(parts[9:]):
                # 快速提取 GT 字段 (冒号前的部分)
                gt_str = sample_data.split(':', 1)[0]
                # 假定都是二倍体
                a1 = gt_str[0]
                a2 = gt_str[-1]
                if a1 == '.' or a2 == '.':
                    stats_matrix[idx][3] += 1 # Missing
                    continue
                if a1 == a2:
                    if a1 == '0':
                        stats_matrix[idx][1] += 1 # HomRef
                    else:
                        stats_matrix[idx][2] += 1 # HomAlt
                else:
                    stats_matrix[idx][0] += 1 # Het


    # 写入结果
    with open(output_file, 'w') as out:
        # 写入表头
        header = ["Sample", "Het_Count", "HomRef_Count", "HomAlt_Count", "Missing_Count", "Heterozygosity_Ratio"]
        out.write('\t'.join(header) + '\n')

        for i, sample in enumerate(samples):
            het = stats_matrix[i][0]
            hom_ref = stats_matrix[i][1]
            hom_alt = stats_matrix[i][2]
            missing = stats_matrix[i][3]
            
            # 计算杂合度 (Het / (HomRef + HomAlt))
            total_hom = hom_ref + hom_alt
            if total_hom > 0:
                ratio = het / total_hom
            else:
                ratio = 0.0

            line_out = [
                sample,
                str(het),
                str(hom_ref),
                str(hom_alt),
                str(missing),
                f"{ratio:.8f}"
            ]
            out.write('\t'.join(line_out) + '\n')
    print(f"统计结果已写入: {output_file}")


if __name__ == "__main__":
    args = parse_args()
    process_vcf(args.input, args.output)
