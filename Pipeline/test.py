import pandas as pd


del_c = pd.read_csv('/home/xinzhou/PROJECT/Edwin_Brca2/del_count.txt', sep='\t')
nor_c = pd.read_csv('/home/xinzhou/PROJECT/Edwin_Brca2/normal_count.txt', sep='\t')
all = pd.read_csv('/home/xinzhou/PROJECT/Edwin_Brca2/data_expression_median.txt', sep='\t')

s = del_c.columns.tolist() + nor_c.columns.tolist()
print(len(del_c.columns.tolist()))
r = list(set(s).intersection(set(all.columns.tolist())))

print(r)
df = all[['Entrez_Gene_Id', 'Description']+r]
df.to_csv('/home/xinzhou/PROJECT/Edwin_Brca2/data_exp.txt', sep='\t', index=False)
d = {}
for i in del_c.columns.tolist():
    d[i] = 'CNA'
for i in nor_c.columns.tolist():
    d[i] = 'WT'
fh = open('/home/xinzhou/PROJECT/Edwin_Brca2/phe.cls', 'w')
fh.write(f'{len(r)} {2} {1}\n')
fh.write(f'# CNA WT\n')
l = ['CNA']*75 + ['WT']*58
sl = ' '.join(l)
fh.write(f'{sl}')
fh.close()