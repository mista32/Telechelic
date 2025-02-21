import os

# 存储每个段中第一列数据大于等于1001和小于等于1000的比例
segment_ratios = []

# 遍历SCluster文件夹下的文件
folder_path = "SCluster"
for filename in os.listdir(folder_path):
    if filename.startswith("F"):
        segments = {}
        with open(os.path.join(folder_path, filename), 'r') as file:
            current_segment = None
            for line in file:
                columns = line.strip().split()
                if len(columns) == 6:
                    sixth_col = int(columns[5])
                    first_col = int(columns[0])
                    if sixth_col != 0:
                        if sixth_col not in segments:
                            segments[sixth_col] = []
                        segments[sixth_col].append(first_col)

            for segment_data in segments.values():
                count_greater_1000 = sum(1 for data in segment_data if data >= 1001)
                count_less_equal_1000 = sum(1 for data in segment_data if data <= 1000)
                if len(segment_data) != 0:
                    segment_ratio = count_greater_1000 / len(segment_data)
                    segment_ratios.append(segment_ratio)

# 计算所有段的比例的平均值
if len(segment_ratios) != 0:
    average_segment_ratio = sum(segment_ratios) / len(segment_ratios)
else:
    average_segment_ratio = 0

# 计算所有文件的平均值的平均值
if len(segment_ratios) != 0:
    total_average_segment_ratio = sum(segment_ratios) / len(segment_ratios)
else:
    total_average_segment_ratio = 0

print("每个段中第一列数据大于等于1001和小于等于1000的比例:", segment_ratios)
print("所有段的比例的平均值:", average_segment_ratio)
print("所有文件的平均值的平均值:", total_average_segment_ratio)