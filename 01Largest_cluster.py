import os

# 存储所有第二列最大值的列表
max_values = []

# 遍历5001个文件
for i in range(5001):
    filename = f"F{i:04d}"
    filepath = os.path.join("/public2/home/yucao/research/Telechelic/SH050_T1500/03Melt/DC028/001/Cluster/SCluster", filename)  # 替换为包含这些文件的文件夹路径

    # 读取文件并提取第二列数据
    with open(filepath, "r") as file:
        lines = file.readlines()
        column2_values = [float(line.split()[1]) for line in lines]

    # 获取第二列数据的最大值并添加到列表中
    max_value = max(column2_values)
    max_values.append(max_value)

# 计算所有最大值的平均值
average_max_value = sum(max_values) / len(max_values)

# 打印平均值
print("5001个文件的第二列最大值的平均值为:", average_max_value)