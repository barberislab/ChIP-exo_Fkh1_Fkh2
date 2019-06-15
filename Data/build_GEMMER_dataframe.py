import pandas as pd
import sqlite3
import ast
import math
import numpy as np
import json

# util functions
def create_connection(db_file):
    """ create a database connection to the SQLite database
        specified by db_file
    :param db_file: database file
    :return: Connection object or None
    """
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Exception as e:
        print(e.message, e.args)
 
    return None

def parse_html(row):
    df_list = pd.read_html(row['All GO terms']) # a list of dataframes

    df = df_list[0]
    df.rename(columns={'Unnamed: 0':'Standard name'}, inplace=True)
    df = df.set_index('Standard name')

    s = df.iloc[0] # series
    d = s.to_dict() # dictionary

    return d


############################
# get data from GEMMER
############################
database = "DB_GEMMER_20180731.db"
conn = create_connection(database)

cursor = conn.execute('SELECT standard_name, systematic_name, go_term_1, go_term_2, go_terms, is_enzyme, yeast7,\
                    GFP_localization, GFP_abundance, expression_peak_phase, expression_peak_time,\
                    name_desc, desc,\
                    KEGG_name, KEGG_KO, KEGG_description, KEGG_pathway FROM genes')
data = [list(x) for x in cursor]

# dataframe of all genes in GEMMER
df_gemmer = pd.DataFrame(data, columns=['Standard name', 'Systematic name', 'Primary GO term', 'Secondary GO term', 
                                  'All GO terms' , 'is enzyme', 'yeast7',
                                  'GFP localization', 'GFP abundance', 'Expression peak phase', 'Expression peak time',
                                  'Name description', 'Description',
                                  'KEGG name', 'KEGG KO', 'KEGG description', 'KEGG pathway'])

# # parse HTML content in 'All GO terms'
print "Parsing HTML table containing GO term counts to a dictionary"
df_gemmer['All GO terms'] = df_gemmer.apply(lambda x: parse_html(x), axis=1)

# force is enzyme to be boolean
df_gemmer['is enzyme'] = df_gemmer.apply(lambda x: bool(x['is enzyme']), axis=1)
df_gemmer['yeast7'] = df_gemmer.apply(lambda x: bool(x['yeast7']), axis=1)


df_gemmer = df_gemmer.set_index('Systematic name')

# shorter name
df = df_gemmer

############################
# Add columns for previous studies for Fkh1,2 and identify which ones they found according to SGD
############################
query = 'SELECT source, target, evidence FROM interactions WHERE SOURCE in (?,?) AND type is ?'
cursor = conn.execute(query, ['FKH1','FKH2']+['regulation'])
data = [list(x) for x in cursor]

# How many regulatory interactions do we have in SGD originating in Fkh1,2?
df_prev = pd.DataFrame(data,columns=['source', 'target', 'evidence']) # previously shown interactions SGD
df_prev_fkh1 = df_prev[df_prev['source'] == 'FKH1']
df_prev_fkh2 = df_prev[df_prev['source'] == 'FKH2']

print 'SGD lists:', len(df_prev), 'existing regulatory interactions for Fkh1,2.'
print 'For Fkh1:', len(df_prev_fkh1)
print 'For Fkh2:', len(df_prev_fkh2)

# get list of targets shown by each study separated for the two Fkh
print 'Per publication we find the following numbers of targets:'

# venters = df_prev[df_prev['evidence'].str.contains('21329885')]
# venters_fkh1 = list(set([x for y in venters[['source','target']][venters['source']=='FKH1'].values for x in y]))
# venters_fkh2 = list(set([x for y in venters[['source','target']][venters['source']=='FKH2'].values for x in y]))

# Get venters from raw 25C UTmax data
# FDR cutoffs: 0.8 (Fkh1), 1.13 (Fkh2)
df_venters = pd.read_excel('Venters_et_al_2011/25C_UTmax.xls',header=0,skiprows=[0,1,2,4,5,6,7,8,9,10,11,12,13], usecols=[0,8,9])
df_venters = df_venters.set_index('Factor')  

venters_fkh1 = df_venters[df_venters['Fkh1'] > 0.8].index.tolist()
venters_fkh2 = df_venters[df_venters['Fkh2'] > 1.13].index.tolist()
print 'Venters \t Fkh1:',len(venters_fkh1),'\t Fkh2:',len(venters_fkh2) 

print 'Missing gene labels from Venters most seem to be retrotransposons:'
print [g for g in venters_fkh1+venters_fkh2 if g not in df.index.tolist()]
print

# remove these from the list
venters_fkh1 = [g for g in venters_fkh1 if g in df.index.tolist()]
venters_fkh2 = [g for g in venters_fkh2 if g in df.index.tolist()]

macIsaac = df_prev[df_prev['evidence'].str.contains('16522208')]
# macIsaac_fkh1 = list(set([x for y in macIsaac[['source','target']][macIsaac['source']=='FKH1'].values for x in y]))
# macIsaac_fkh2 = list(set([x for y in macIsaac[['source','target']][macIsaac['source']=='FKH2'].values for x in y]))
macIsaac_fkh1 = list(set(macIsaac[macIsaac['source']=='FKH1']['target'].values))
macIsaac_fkh2 = list(set(macIsaac[macIsaac['source']=='FKH2']['target'].values))

# map standard names to systematic names
macIsaac_fkh1 = [df_gemmer[df_gemmer['Standard name'] == stdname].index[0] for stdname in macIsaac_fkh1 ]
macIsaac_fkh2 = [df_gemmer[df_gemmer['Standard name'] == stdname].index[0] for stdname in macIsaac_fkh2 ]

print 'MacIsaac \t Fkh1:',len(macIsaac_fkh1),'\t Fkh2:',len(macIsaac_fkh2) 

# add new columns for these 3 studies
df['Venters 2011 Fkh1'] = False
df['Venters 2011 Fkh2'] = False
df['Venters 2011 Fkh1'].loc[venters_fkh1] = True
df['Venters 2011 Fkh2'].loc[venters_fkh2] = True

df['MacIsaac 2006 Fkh1'] = False
df['MacIsaac 2006 Fkh2'] = False
# df.loc[df['Standard name'].isin(macIsaac_fkh1),'MacIsaac 2006 Fkh1'] = True
# df.loc[df['Standard name'].isin(macIsaac_fkh2),'MacIsaac 2006 Fkh2'] = True
# df['MacIsaac 2006 Fkh1'].loc[macIsaac_fkh1] = True
# df['MacIsaac 2006 Fkh2'].loc[macIsaac_fkh2] = True
df.at[macIsaac_fkh1,'MacIsaac 2006 Fkh1'] = True
df.at[macIsaac_fkh2,'MacIsaac 2006 Fkh2'] = True

# due to GEMMER only having protein-coding genes we are missing tK(UUU)P which was shown by MacIsaac et al.
df.loc['tK(UUU)P'] = None
df.loc['tK(UUU)P'] = None
df.at[['tK(UUU)P'],'MacIsaac 2006 Fkh1'] = True
df.at[['tK(UUU)P'],'MacIsaac 2006 Fkh2'] = True
df.at[['tK(UUU)P'],'is enzyme'] = False
df.at[['tK(UUU)P'],'yeast7'] = False

print(df.loc['tK(UUU)P'])

# Read in the table of targets from Ostrow et al.
df_ostrow = pd.read_excel("./Ostrow_et_al_2014/Table_S4.xlsx",header=None)
df_ostrow.drop([0,1,2,4],axis=1,inplace=True)

# separate Fkh1,2 results
df_ostrow_fkh1 = df_ostrow[df_ostrow[5] == 'Fkh1'][3].values
df_ostrow_fkh2 = df_ostrow[df_ostrow[5] == 'Fkh2'][3].values

print 'Ostrow et al. found:', len(df_ostrow_fkh1), 'targets for Fkh1 and',len(df_ostrow_fkh2),'targets for Fkh2'

# update dataframe with ostrow et al. results
df['Ostrow 2014 Fkh1'] = False
df['Ostrow 2014 Fkh2'] = False
df.loc[df.index.isin(df_ostrow_fkh1),'Ostrow 2014 Fkh1'] = True
df.loc[df.index.isin(df_ostrow_fkh2),'Ostrow 2014 Fkh2'] = True

# force is enzyme to be boolean
df = df.fillna({'is enzyme':False,'yeast7':False,'Standard name':'', 'MacIsaac 2006 Fkh1':False, 'MacIsaac 2006 Fkh2':False, 'Venters 2011 Fkh1':False, 'Venters 2011 Fkh2':False, 'Ostrow 2014 Fkh1':False, 'Ostrow 2014 Fkh2':False})
df['is enzyme'] = df.apply(lambda x: bool(x['is enzyme']), axis=1)
df['yeast7'] = df.apply(lambda x: bool(x['yeast7']), axis=1)
# df = df.apply(lambda x: bool(x[['is enzyme','MacIsaac 2006 Fkh1']]), axis=1)
bool_cols = ['is enzyme','yeast7','MacIsaac 2006 Fkh1','MacIsaac 2006 Fkh2','Venters 2011 Fkh1','Venters 2011 Fkh2','Ostrow 2014 Fkh1','Ostrow 2014 Fkh2']
df[bool_cols] = df[bool_cols].astype(bool)



############################
# Order and store the dataset in pickled format
############################
df = df[['Standard name','Primary GO term','Secondary GO term', 'All GO terms',
        'Name description','Description','GFP abundance','GFP localization','Expression peak phase',
        'Expression peak time','is enzyme','yeast7',
        'MacIsaac 2006 Fkh1', 'Venters 2011 Fkh1', 'Ostrow 2014 Fkh1',
        'MacIsaac 2006 Fkh2', 'Venters 2011 Fkh2', 'Ostrow 2014 Fkh2',
        'KEGG name', 'KEGG KO', 'KEGG description', 'KEGG pathway']]
df = df.sort_values(by='Standard name')

print(df.loc['tK(UUU)P'])


df.to_pickle("DB_GEMMER_20180731")