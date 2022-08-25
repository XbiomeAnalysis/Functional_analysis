import argparse
import pandas as pd 

parser = argparse.ArgumentParser()
parser.add_argument("--input","-i",dest="input",help="input file genefamilies")
parser.add_argument("--knumber","-k",dest="knumber",help="knumber file, eg: samplename_knumber.txt")
parser.add_argument("--minpath","-m",dest="minpath",help="minpath file, eg: samplename.input")
args = parser.parse_args()

df = pd.read_table(args.input,sep="\t")

# 提取最高层级
lst_index = []
for i in df["# Gene Family"].values:
	if len(i.split("|")) > 1:
		lst_index.append(False)
	else:
		lst_index.append(True)
selected_data = df.loc[lst_index,:]
selected_data = selected_data.set_index("# Gene Family")

# 相对丰度化
## 暂时不对KNumbers的RPK值进行相对丰度化，这一步在后续分析中会进行
#selected_data_relative = selected_data.apply(lambda x: x/sum(x),axis=0)

# 提取knumber条目
k_lst = []
for i in selected_data.index.values:
	if i.startswith("K"):
		k_lst.append(True)
	elif i.startswith("UNMAPPED"):
		k_lst.append(True)
	else:
		k_lst.append(False)
selected_k = selected_data.loc[k_lst,:]

selected_k.to_csv(args.knumber,header=True,sep="\t")

with open(args.minpath,"w") as f:
	for i in range(selected_k.shape[0]):
		f.write(str(i) + "\t" + selected_k.index.tolist()[i] + "\n")
	f.close
