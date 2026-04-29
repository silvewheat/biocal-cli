import argparse
import sys

def count_lines(filepath):
    """
    第一遍扫描：计算文件总行数。
    使用 buffer 读取以提高大文件读取速度。
    """
    count = 0
    buffer_size = 1024 * 1024 * 1024  # 1GB buffer
    try:
        with open(filepath, 'rb') as f:
            while True:
                buffer = f.read(buffer_size)
                if not buffer:
                    break
                count += buffer.count(b'\n')
        # 处理最后一行可能没有换行符的情况（虽然 SNP 文件通常都有）
        # 但这种简单的二进制计数通常足够精确用于行抽样
        return count
    except Exception as e:
        print(f"Error checking file size: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="大文件均匀抽样工具 (Streaming Mode)")
    parser.add_argument('-i', '--input', required=True, help="输入文件路径")
    parser.add_argument('-o', '--output', required=True, help="输出文件路径")
    parser.add_argument('-n', '--number', type=int, required=True, help="需要抽取的行数")
    
    args = parser.parse_args()

    input_path = args.input
    output_path = args.output
    target_n = args.number

    print(f"[Phase 1] 正在扫描文件以统计总行数: {input_path} ...")
    total_lines = count_lines(input_path)
    
    if total_lines == 0:
        print("错误：文件为空。")
        sys.exit(0)

    print(f" -> 检测到总行数: {total_lines}")

    # 边界情况处理：如果请求数量 >= 总行数，直接复制
    if target_n >= total_lines:
        print(f"警告: 目标数量 ({target_n}) >= 总行数，将复制整个文件。")
        with open(input_path, 'r', encoding='utf-8', errors='ignore') as f_in, \
             open(output_path, 'w', encoding='utf-8') as f_out:
            for line in f_in:
                f_out.write(line)
        return

    print(f"[Phase 2] 正在均匀抽取 {target_n} 行数据...")
    
    with open(input_path, 'r', encoding='utf-8', errors='ignore') as f_in, \
         open(output_path, 'w', encoding='utf-8') as f_out:
        
        # 核心算法：如何决定当前行是否抓取？
        # 我们需要抓取 N 个点，分布在 0 到 Total-1 之间。
        # 第 k 个采样点的理想索引位置 = floor(k * Total / Target)
        
        samples_taken = 0
        current_line_idx = 0
        
        # 计算第一个需要抓取的目标行索引
        # 公式：next_target_idx = int(samples_taken * total_lines / target_n)
        # 初始 samples_taken = 0, 所以第一个目标是第 0 行
        next_target_idx = 0 

        for line in f_in:
            # 如果当前行号匹配目标行号
            if current_line_idx == next_target_idx:
                f_out.write(line)
                samples_taken += 1
                
                # 如果已经抽够了，停止读取（节省时间）
                if samples_taken >= target_n:
                    break
                
                # 计算下一个目标行号
                # 使用乘法而非累加步长，是为了避免浮点数累加产生的精度误差
                next_target_idx = int(samples_taken * total_lines / target_n)
            
            current_line_idx += 1

    print(f"完成！已成功抽取 {samples_taken} 行并保存至: {output_path}")

if __name__ == "__main__":
    main()