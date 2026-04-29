import phlash
import pickle
import os
import argparse  # 用于解析命令行参数

def load_chromosome_lengths(genome_file):
    """加载染色体长度，返回 {chr1: 长度, chr2: 长度, ...} 格式"""
    chrom2len = {}
    with open(genome_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                chrom_name = parts[0]
                chrom_length = parts[1]
                chrom2len[chrom_name] = chrom_length
    return chrom2len

def process_sample(vcffile, sample, chrom2len, params):
    """
    处理单个指定样本
    params: 包含所有运行参数的命名空间对象
    """
    print(f"=====================================")
    print(f"开始处理样本: {sample}")
    print(f"使用VCF文件: {vcffile}")
    print(f"=====================================")
    
    # 收集该样本的所有染色体数据
    chroms_data = []
    chromosomes = range(params.chrom_start, params.chrom_end + 1)
    for chr_num in chromosomes:
        chrom_name = f"chr{chr_num}"  # 染色体名称：chr1, chr2...
        try:
            # 检查染色体长度是否存在
            if chrom_name not in chrom2len:
                print(f"  警告：{chrom_name} 的长度在基因组文件中未找到，跳过")
                continue
            
            # 获取染色体长度
            chrom_length = chrom2len[chrom_name]
            # 提取单个样本单个染色体的数据（指定完整区域）
            contig_data = phlash.contig(
                vcffile,
                samples=sample.strip().split(','),
                region=f"{chrom_name}:1-{chrom_length}"
            )
            chroms_data.append(contig_data)
            print(f"  成功加载 {chrom_name}:1-{chrom_length}")
        except Exception as e:
            print(f"  警告：加载 {chrom_name} 时出错 - {str(e)}")
            continue  # 跳过出错的染色体
    
    # 检查是否有有效数据
    if not chroms_data:
        print(f"  ❌ 样本 {sample} 没有有效染色体数据，跳过拟合")
        return
    
    # 拟合模型
    print(f"\n  开始拟合 {sample} 的数据...")
    test_data = chroms_data[0]
    train_data = chroms_data[1:]
    results = phlash.fit(
        data=train_data,
        test_data=test_data,
        mutation_rate=params.mutation_rate
    )
    
    # 创建输出目录（如果不存在）
    os.makedirs(params.outdir, exist_ok=True)
    
    # 保存结果到pickle文件
    output_path = os.path.join(params.outdir, f"phlash_results_{sample}.pkl")
    with open(output_path, 'wb') as f:
        pickle.dump(results, f)
    
    print(f"  ✅ 样本 {sample} 处理完成！")
    print(f"  结果已保存到: {output_path}")
    print(f"=====================================\n")
    return output_path

def main():
    # 解析命令行参数（所有默认值直接在这里定义）
    parser = argparse.ArgumentParser(
        description='Phlash单样本分析工具 - 支持命令行指定VCF文件和样本ID',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # 必选参数
    parser.add_argument('--vcf', required=True, help='输入的VCF文件路径（.gz，需有索引）')
    parser.add_argument('--samples', required=True, help='样本ID,多样本用逗号分割（sample1,sample2,sample3）')
    parser.add_argument('--genome', required=True, help='基因组长度文件路径（格式：每行 "chromosome_name length"）')
    
    # 可选参数（默认值直接在参数中定义）
    parser.add_argument(
        '--outdir', 
        default="./phlash_results", 
        help='结果输出目录（默认：./phlash_results）'
    )
    parser.add_argument(
        '--mutation-rate',
        type=float,
        default=1.17e-8,
        help='突变率（默认：1.17e-08）'
    )
    parser.add_argument(
        '--chrom-start', 
        type=int, 
        default=1, 
        help='起始染色体编号（默认：1）'
    )
    parser.add_argument(
        '--chrom-end', 
        type=int, 
        default=29, 
        help='结束染色体编号（默认：29）'
    )
    
    args = parser.parse_args()
    
    # 打印本次运行的所有参数信息
    print("=====================================")
    print("Phlash 单样本分析 - 运行参数汇总")
    print("=====================================")
    print(f"输入VCF文件:    {args.vcf}")
    print(f"分析样本ID:     {args.samples}")
    print(f"基因组长度文件: {args.genome}")
    print(f"结果输出目录:   {args.outdir}")
    print(f"突变率:         {args.mutation_rate}")
    print(f"染色体范围:     chr{args.chrom_start} ~ chr{args.chrom_end}")
    print("=====================================\n")
    
    # 加载染色体长度
    print(f"正在加载染色体长度信息...")
    chrom2len = load_chromosome_lengths(args.genome)
    print(f"成功加载 {len(chrom2len)} 条染色体的长度信息\n")
    
    # 处理指定的单个样本（将所有参数通过args传递）
    try:
        process_sample(args.vcf, args.samples, chrom2len, args)
    except Exception as e:
        print(f"❌ 样本 {args.samples} 处理失败: {str(e)}")
        exit(1)
    
    print("🎉 所有分析完成！")

if __name__ == '__main__':
    main()