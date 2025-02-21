import os

# 定义需要统计的第一列数据的值
targets = [100, 300, 500, 700, 900, 1100, 1300, 1500, 1700, 1900]

result_array = [[0] * len(targets) for _ in range(5001)]

# 遍历文件夹中的文件
folder_path = './SCluster'
for file_index, filename in enumerate(sorted(os.listdir(folder_path))):
    if filename.startswith('F'):  # 简单的文件名验证       and filename.isdigit()[1:]
        file_path = os.path.join(folder_path, filename)
        
        with open(file_path, 'r') as file:
            for line_index, line in enumerate(file):
                if line.strip():
                    data = list(map(float, line.split()))

                    if data[0] in targets:
                        target_index = targets.index(data[0])
                        
                        if data[5] == 0:
                            result_array[file_index][target_index] = 0
                        else:
                            result_array[file_index][target_index] = 1

# 打印结果数组（或者你可以保存到文件中）
for row in result_array:
    print(row)

with open('Sticky.txt', 'w') as output_file:
    for row in result_array:
        output_file.write(' '.join(map(str, row)) + '\n')