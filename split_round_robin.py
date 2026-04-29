import argparse
import sys
import os
from itertools import cycle

def main():
    parser = argparse.ArgumentParser(description="文件轮询拆分工具 (Round-Robin Splitter)")
    parser.add_argument('-i', '--input', required=True, help="输入文件路径")
    parser.add_argument('-p', '--prefix', required=True, help="输出文件的前缀 (例如: split_out)")
    parser.add_argument('-n', '--number', type=int, required=True, help="要拆分成多少个文件")
    
    args = parser.parse_args()

    input_path = args.input
    num_splits = args.number
    output_prefix = args.prefix

    if num_splits <= 0:
        print("错误：拆分数量必须大于 0")
        sys.exit(1)

    # 1. 自动推断扩展名
    # 如果输入是 'data.vcf'，ext 就是 '.vcf'
    # 如果输入没有扩展名，ext 就是空字符串
    _, ext = os.path.splitext(input_path)

    print(f"准备将文件拆分为 {num_splits} 份...")

    # 2. 打开所有输出文件的句柄
    file_handles = []
    output_filenames = []
    
    try:
        for i in range(num_splits):
            # 自动构建文件名: 前缀 + _ + 编号 + 原扩展名
            # 例如: prefix_1.txt
            fname = f"{output_prefix}_{i+1}{ext}"
            
            # 使用 'w' 模式打开
            f = open(fname, 'w', encoding='utf-8')
            
            file_handles.append(f)
            output_filenames.append(fname)
            
        print(f"已创建输出文件: {output_filenames[0]} ... {output_filenames[-1]}")

    except OSError as e:
        print(f"创建文件失败 (可能是打开文件过多): {e}")
        # 清理已打开的文件
        for f in file_handles:
            if f: f.close()
        sys.exit(1)

    # 3. 核心循环：轮询分发
    # 创建一个无限循环的迭代器：File1 -> File2 -> File3 -> File1 ...
    file_cycler = cycle(file_handles)
    
    count = 0
    try:
        with open(input_path, 'r', encoding='utf-8', errors='ignore') as f_in:
            for line in f_in:
                # 获取下一个文件句柄并写入
                target_file = next(file_cycler)
                target_file.write(line)
                count += 1
                
                # 简单的进度提示
                if count % 1000000 == 0:
                    print(f"已处理 {count} 行...", end='\r')

    except Exception as e:
        print(f"\n处理过程中出错: {e}")
    finally:
        # 4. 关闭所有文件
        for f in file_handles:
            f.close()

    print(f"\n完成！总共处理了 {count} 行。")
    print("生成的文件列表:")
    for name in output_filenames:
        print(f" - {name}")

if __name__ == "__main__":
    main()