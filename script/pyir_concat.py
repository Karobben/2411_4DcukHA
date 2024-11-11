import pandas as pd


# read the tables
TB_VH = pd.read_csv('VH.tsv', sep = '\t', index_col = 0)
TB_VL = pd.read_csv('VL.tsv', sep = '\t', index_col = 0)
#
# Remove on annotated results
TB_VL = TB_VL[TB_VL.productive== "T"]
TB_VH = TB_VH[TB_VH.productive== "T"]

TB_VL = TB_VL.reset_index()
TB_VH = TB_VH.reset_index()

# find the ID list
TB_VL = pd.concat([TB_VL, pd.DataFrame([i.split("_",1) for i in TB_VL.sequence_id], columns = ['cell', 'ConID'])], axis = 1)
TB_VH = pd.concat([TB_VH,pd.DataFrame([i.split("_",1) for i in TB_VH.sequence_id], columns = ['cell', 'ConID'])], axis = 1)

ID1 = TB_VL.cell.value_counts().index[TB_VL.cell.value_counts() >1].to_list()
ID2 = TB_VH.cell.value_counts().index[TB_VH.cell.value_counts() >1].to_list()
ID_lst = set(ID1+ID2)


# save the sequences which has more than 1 result
TB_VL[TB_VL.cell.isin(ID_lst)].to_csv('VL_more.csv')
TB_VH[TB_VH.cell.isin(ID_lst)].to_csv('VH_more.csv')

# pair the chain 
TB_VL = TB_VL[~TB_VL.cell.isin(ID_lst)]
TB_VH = TB_VH[~TB_VH.cell.isin(ID_lst)]

pd.merge(TB_VH, TB_VL, on= 'cell', how = 'outer', suffixes = ["_H", "_L"]).to_csv('pyir_contact_result.csv')
