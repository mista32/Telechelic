import os

# 存储每个文件的dangle变量和loop变量
realfree_values = []
dangle_values = []
loop_values = []

# 遍历SCluster文件夹下的文件
folder_path = "SCluster"
for filename in os.listdir(folder_path):
    if filename.startswith("F"):
        segments = {}
        segments0 = []
        realfree = 0
        dangle = 0
        loop = 0
        with open(os.path.join(folder_path, filename), 'r') as file:
            current_segment = None
            for line in file:
                columns = line.strip().split()
                if len(columns) == 6:
                    first_col = int(columns[0])
                    sixth_col = int(columns[5])

                    # 第六列数据为0且第一列小于等于1000的行
                    if sixth_col == 0 and first_col <= 1000:
                        dangle += 1                             # 统计单离子的总个数
                        segments0.append(first_col)             # 将单个离子的簇挑选出来并统计离子ID
                        continue

                    if sixth_col not in segments:              
                        segments[sixth_col] = []
                    segments[sixth_col].append(first_col)

            for segment_data in segments.values():
                consecutive_odd_even = False                    #在这些离子簇中  如果 连续的奇数+偶数   则为loop链 
                used_pairs = set()
                for data in segment_data:
                    if data <= 1000:
                        if data % 2 == 1:
                            if data + 1 in segment_data and (data, data + 1) not in used_pairs:
                                loop += 2
                                used_pairs.add((data, data + 1))
            #print(segments0)
            for segment0_data in segments0:
                used_pairs0 = set()                             # 在这些单离子中  如果有连续的奇数+偶数  则为free链  否则为dangling链
                if segment0_data % 2 == 1:
                    if segment0_data + 1 in segments0 and (data, data + 1) not in used_pairs0:
                        realfree += 2
                        used_pairs0.add((segment0_data, segment0_data + 1))

        realfree_values.append(realfree)
        dangle_values.append(dangle)
        loop_values.append(loop)

# 计算所有文件的dangle变量和loop变量的平均值
# print(len(dangle_values), len(loop_values))

average_realfree = sum(realfree_values) / len(realfree_values)

average_dangling = (sum(dangle_values) - sum(realfree_values)) / len(realfree_values)   

average_loop = sum(loop_values) / len(loop_values)

# print("所有文件变量的平均值:", average_free, average_loop, average_realdangle, 1000-(average_realdangle + average_free + average_loop))
print("所有文件变量的平均值: {:>8.7f} {:>8.7f} {:>8.7f} {:>8.7f}".format(average_realfree/1000, average_loop/1000, average_dangling/1000, (1000-(average_realfree + average_dangling + average_loop))/1000))
