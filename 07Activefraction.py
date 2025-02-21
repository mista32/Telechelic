import os

# 文件夹路径
folder_path = './SCluster'

# 初始化计数器
total_files = 0
total_ratio_sum = 0.0

# 遍历文件夹中的所有文件
for filename in os.listdir(folder_path):
    if filename.startswith('F'):  # 确保文件名以F开头且后续为数字    and filename.isdigit()
        file_path = os.path.join(folder_path, filename)
        total_files += 1
        
        # 初始化计数器
        count_zero_sixth_column = 0
        
        # 读取文件并统计第6个数据为0的行数
        with open(file_path, 'r') as file:
            for line in file:
                data = line.strip().split()  # 假设数据以空格分隔
                if len(data) >= 6 and data[5] == '0':  # 确保数据至少有6列且第6列是0
                    count_zero_sixth_column += 1
        
        # 计算每个文件的比例 m/2000
        ratio = count_zero_sixth_column / 2000.0
        print(ratio)
        # 累加比例
        total_ratio_sum += ratio

# 计算比例的平均值
average_ratio = total_ratio_sum / total_files

print(f"5001个文件中，第六个数据为0的行比例的平均值为: {average_ratio}")