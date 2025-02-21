import os

# 存储每个Fxxxx文件遍历后的平均值
averages = []

# 遍历SCluster文件夹下的文件
folder_path = "SSCluster"
for filename in os.listdir(folder_path):
    if filename.startswith("F"):
        sixth_col_to_second_col = {}
        with open(os.path.join(folder_path, filename), 'r') as file:
            prev_sixth_col = None
            total_sum = 0
            count = 0
            for line in file:
                columns = line.strip().split()
                if len(columns) == 6:
                    sixth_col = int(columns[5])
                    second_col = int(columns[1])
                    if second_col > 3:  # 仅当第二列数据大于3时进行后续操作
                        if sixth_col != 0:
                            if prev_sixth_col is None:
                                sixth_col_to_second_col[sixth_col] = second_col
                            elif sixth_col != prev_sixth_col:
                                sixth_col_to_second_col[sixth_col] = sixth_col_to_second_col.get(sixth_col, 0) + second_col
                            prev_sixth_col = sixth_col
                            total_sum += second_col
                            count += 1

            if count != 0:
                total_sum = sum(sixth_col_to_second_col.values())
                num_unique_values = len(sixth_col_to_second_col)
                if num_unique_values != 0:
                    average = total_sum / num_unique_values
                else:
                    average = 0
                averages.append(average)

# 计算所有Fxxxx文件遍历后的平均值的平均值
if len(averages) != 0:
    total_average = sum(averages) / len(averages)
else:
    total_average = 0

print("每个Fxxxx文件遍历后的平均值:", averages)
print("所有Fxxxx文件遍历后的平均值的平均值:", total_average)
