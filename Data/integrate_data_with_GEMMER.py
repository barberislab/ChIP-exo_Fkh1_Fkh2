import pandas as pd
import sqlite3
import ast
import math
import numpy as np
import json

def set_pathview_score(row, case, cutoff):
    
    if row["maxPeak_AD_12 "+case] >= cutoff:
        Mondeel = True
    else:
        Mondeel = False
    Venters = row['Venters 2011 ' + case[:4]]
    Ostrow = row['Ostrow 2014 ' + case[:4]]
    MacIsaac = row['MacIsaac 2006 ' + case[:4]]
        
    score = None 
    if sum([Mondeel, Venters, Ostrow, MacIsaac]) > 0:
        if sum([Mondeel, Venters, Ostrow]) == 0: # found in MacIsaac
            score = -1
        if sum([Mondeel, Venters, Ostrow]) == 1: # found 1x with ChIP
            score = 0
        if sum([Mondeel, Venters, Ostrow]) >= 2: # verified at least 2x with ChIP
            score = 1

    return score

############################
# SETTINGS
############################
### Choose a promoter setting below for TSS/ORF
promoter = 'ORF/'

### Determine percentiles per TF (LOG+STAT) or separately for each experiment
normalization_setting = 'condition_normalized/'


############################
# LOAD DATASETS
############################
# GOAL: INTEGRATE ALL DIFFERENT DATASETS INTO ONE DATAFRAME MERGED WITH GEMMER DATA

# start from GEMMER
df = pd.read_pickle('DB_GEMMER_20180731')

# All files share the same structure: column 4 (3 with zero indexing) is the gene name and 4 the quantification

# for the EBICB manuscript we use ORF and include loessPeak and OL7, OL10, AD7 and AD10
if promoter == 'ORF/':
    filenames = [ 
    'maxPeak_AD_12_Fkh1_log.bed', 'maxPeak_AD_12_Fkh1_stat.bed', 'maxPeak_AD_12_Fkh2_log.bed', 'maxPeak_AD_12_Fkh2_stat.bed',
    'MACE_Fkh1_log.bed', 'MACE_Fkh1_stat.bed', 'MACE_Fkh2_log.bed', 'MACE_Fkh2_stat.bed',
    'GEM_Fkh1_log.bed', 'GEM_Fkh1_stat.bed', 'GEM_Fkh2_log.bed', 'GEM_Fkh2_stat.bed',]


### Integrate data from all files
column_name_list = []
for file in filenames: 
    file = file[:-4] # remove .bed

    method = '_'.join(file.split('_')[:-2])
    condition = ' '.join(file.split('_')[-2:])
    column_name = method + ' ' + condition
    column_name_list.append(column_name)
    
    if 'GEM' in file or 'MACE' in file:
        path = './'+promoter+file+'.bed'
    else:
        path = './'+promoter+normalization_setting+file+'.bed'

    if 'GEM' in file:
        df_new = pd.read_csv(path, sep="\t")
        df_new = df_new.set_index('GeneSys')
        df_new = df_new.rename(columns={condition: column_name, 'Motif': 'GEM motif ' + condition}) # rename to identify method + condition + TF
        df_new = df_new[[column_name,'GEM motif ' + condition]] # drop other columns

    else:
        df_new = pd.read_csv(path, sep="\t", header=None)
        df_new = df_new.set_index(3)
        df_new = df_new.rename(columns={4: column_name}) # rename to identify method + condition + TF
        df_new = df_new[[column_name]] # drop other columns


    # merge on indices but keep only the ones in the existing df
    df = df.merge(df_new, how='outer', left_index=True, right_index=True) # outer = union of indices. 'left' would leave out genes not in GEMMER

    # if a gene occurs multiple times -> multiple peaks -> pick largest one
    # sort by the new column
    if 'MACE' in file:
        df = df.sort_values(column_name, ascending=True) # NOTE THAT FOR MACE WE NEED TO SORT ASCENDING. LOWER P-VALUE IS GOOD!
    else:
        df = df.sort_values(column_name, ascending=False) # here higher SNR is good so descending.
    
    # index.duplicates returns a Boolean ndarray with True for duplicates except the 'first' entry
    df = df[~df.index.duplicated(keep='first')]
    

############################
# FIX NAN ERRORS
############################ 
# NON-PROTEIN-CODING GENES (NOT FROM GEMMER) ARE NOT ENZYMES
# NO NAN IN STANDARD NAME
# MACISAAC, VENTERS AND OSTROW DEFAULT TO FALSE
df = df.fillna({'Primary GO term':'None','Secondary GO term':'None','is enzyme':False,'yeast7':False,'Standard name':'', 'MacIsaac 2006 Fkh1':False, 'MacIsaac 2006 Fkh2':False, 'Venters 2011 Fkh1':False, 'Venters 2011 Fkh2':False, 'Ostrow 2014 Fkh1':False, 'Ostrow 2014 Fkh2':False})
    
df['Standard name'] = [sysname if df.loc[sysname]['Standard name'] == '' else  df.loc[sysname]['Standard name'] for sysname in df.index.tolist()]
    
############################
# set discrete score for KEGG pathview consensus view
############################
# This is used to make the consensus pathview images: 0x,1x,2x shown by ChIP evidence
df['KEGG pathview Fkh1 log'] = df.apply(lambda x: set_pathview_score(x, 'Fkh1 log', 1), axis=1)
df['KEGG pathview Fkh2 log'] = df.apply(lambda x: set_pathview_score(x, 'Fkh2 log', 1), axis=1)


############################
# Order and store the dataset in pickled format
############################
df = df[['Standard name','Primary GO term','Secondary GO term','GFP abundance','GFP localization','Expression peak phase',
        'Expression peak time','is enzyme','yeast7',
        'MacIsaac 2006 Fkh1', 'Venters 2011 Fkh1', 'Ostrow 2014 Fkh1',
        'MacIsaac 2006 Fkh2', 'Venters 2011 Fkh2', 'Ostrow 2014 Fkh2',
        'KEGG name', 'KEGG KO', 'KEGG description', 'KEGG pathway','KEGG pathview Fkh1 log','KEGG pathview Fkh2 log','All GO terms','Name description',
        'Description'] + column_name_list + ['GEM motif ' + condition for condition in ['Fkh1 log', 'Fkh1 stat', 'Fkh2 log', 'Fkh2 stat']]]
df = df.sort_values(by='Standard name', ascending=True)

df.to_pickle("./df_"+promoter[:-1]+"_"+normalization_setting[:-1])


############################
# Generate the Excel file for the dataset
############################
writer = pd.ExcelWriter("./df_"+promoter[:-1]+"_"+normalization_setting[:-1]+".xlsx", engine='xlsxwriter')
workbook = writer.book

format_null = workbook.add_format({'text_wrap': True,'align':'left','font_size':10})
format_lightred = workbook.add_format({'text_wrap': True,'bg_color': '#ffcccc',
                                'font_color': '#000000','border':1,'align':'left','font_size':10})


#############################
# 0
#############################
df.to_excel(writer,sheet_name='All', index=True)
worksheet = writer.sheets['All']
worksheet.set_column('A:B',20) # Systematic name - Standard name
worksheet.set_column('C:R',17) # Primary GO - KEGG KO
worksheet.set_column('S:T',60) # KEGG description + KEGG pathways
worksheet.set_column('U:V',20) # 2x Pathview
worksheet.set_column('W:X',60) # go terms count + name description
worksheet.set_column('Y:Y',150) # description
worksheet.set_column('Z:AZ',25) # Quantification 

# freeze first row and column 1+2
worksheet.freeze_panes(1, 2)

# format rows matching a condition
# row is the row number, row_data is a tuple where the second element is the data series
for row in range(len(df.index)):
    row_data = df.iloc[row]
    worksheet.set_row(row + 1, None, format_null)


# SAVE
writer.save()