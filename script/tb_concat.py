import pandas as pd

def AddCellID(TB1, HEAD):
    TB1['cell_id'] = [f"{HEAD}_{i.split('_')[0]}" for i in TB1['seq_id']]
    return TB1

TB = pd.read_csv('result/GEX_anno.csv', index_col = 0)
TB1 = pd.read_csv('result/duck_G1_VDJ.tsv', sep ='\t', index_col = 0)
TB2 = pd.read_csv('result/duck_Group3_VDJ.tsv', sep ='\t', index_col = 0)
TB3 = pd.read_csv('result/miniHA3_G1_VDJ.tsv', sep ='\t', index_col = 0)
TB4 = pd.read_csv('result/miniHA3_G3_VDJ.tsv', sep ='\t', index_col = 0)

TB['cell_id'] = TB.index
TB1 = AddCellID(TB1, 'duck_G1')
TB2 = AddCellID(TB2, 'duck_Group3')
TB3 = AddCellID(TB3, 'miniHA3_G1')
TB4 = AddCellID(TB4, 'miniHA3_G3')

TB_all = pd.concat([TB1, TB2, TB3, TB4])
TB_final = pd.merge(TB_all, TB[["cell_id", 'seurat_clusters']], how = 'left')
TB_final = TB_final.sort_values(by=["cell_id", 'seq_id'])

L_lst = TB_final.cell_id[TB_final.chain_type =='L'].to_list()
H_lst = TB_final.cell_id[TB_final.chain_type =='H'].to_list()
C_lst = [i for i in H_lst if i in L_lst]

TB_final = TB_final[TB_final.cell_id.isin(C_lst)]
TB_final.to_csv('result/Duck_numbring_class.tsv', sep = '\t')
