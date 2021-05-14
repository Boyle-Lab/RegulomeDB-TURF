import json
import pandas as pd
import pyBigWig
import numpy as np
import pickle
from collections import defaultdict
import re
import sys

input_json = sys.argv[1]
input_bed = sys.argv[2]
organ_list = sys.argv[3]
outfile = sys.argv[4]
path = 'Precalculated_features/' #directory for precaculated features, see README for the link to download
####################################################################
###Function to output RegDB generic features from json object###
def get_generic_features(var_json):
    '''Function to output RegDB generic features from a json object
    Args:
        var_json(json): a json object contains RegDB features
    Returns:
        out(list of str): RegDB features ('chrom','end','CHIP','DNASE','PWM','FOOTPRINT','QTL','PWM_matched','FOOTPRINT_matched')
    '''
    out = [var_json['chrom'],var_json['end']]
    feature_names = ['ChIP','DNase','PWM','Footprint','QTL','PWM_matched','Footprint_matched']
    out_features = [int(var_json['features'][feature]) for feature in feature_names] #convert to binary values
    out.extend(out_features)
    return out

###Function to get RegulomeDB organSp features DNase & FOOTPRINT from json object###
def get_organ_sp_RegDB_features(var_json,organs_query,biosample_organ_mapping = '/data/data_repo/shengchd/Data_resource/ENCODE/biosample_organ.txt'):
    '''Function to output RegDB organ specific features from a json object of a variant
    Args:
        var_json(json): a json object contains RegDB features & peaks information
        organs_query (list of str): biosample organ name(s) (e.g. ['bodily fluid','blood','brain'])
        biosample_organ_mapping (str): file containing mapping between biosample name to organ name
    Returns:
        feature_values (dict of list): key: organ names in input query organ list; values: 'chrom','end','DNASE','FOOTPRINT'
    '''
    from collections import defaultdict
    import re
    biosample_organ_dict = defaultdict(str)
    #get mapping from biosample to organ
    with open(biosample_organ_mapping,'r') as f:
        for line in f:
            if len(line.strip().split('\t')) < 2: #biosample does not have corresponding organ name
                next
            else:
                biosample,organ_name = line.strip().split('\t')
                biosample_organ_dict[biosample] = organ_name #organ_name could contain multiple organs, split by ','
    ####manually add 'H1-hESC' to biosample dict (named as 'H1' in mapping file, but 'H1-hESC' in json object from .out)###
    biosample_organ_dict['H1-hESC'] = 'embryo'
    
    feature_values = defaultdict(list) #output dict
    method_names = ['DNase-seq','Footprints']
    for organ_query in organs_query:
        results = defaultdict(int) #keys: method_names; default value is 0
        for exp in var_json['peaks']:
            if exp['method'] in method_names:
                hit_organs = biosample_organ_dict[exp['biosample_term_name']].split(',') #if no mapping to organ, then return ['']; any non-empty string not in [''];
                if organ_query in hit_organs:
                    results[exp['method']] = 1
        out = [var_json['chrom'],var_json['end']]
        out.append(results['DNase-seq'])
        out.append(results['Footprints'])
        feature_values[organ_query] = out
    return feature_values

###Function to get records from SQL database on Histone marks
def query_RegDB_histone(chrom, pos, DB_DIR):
    '''Function to query a single base on genome over RegDB SQL tables
    Args:
        chrom (int or str): query chromosome (1..22, X, Y)
        pos (int): query position (0-based)
        DB_DIR (str): directory of SQL database on Histone marks
    Returns:
        hits (dict): query results including all annotations
    '''
    
    import sqlite3
    from collections import defaultdict
    
    #remove 'chr' if exits in input chrom:
    import re
    chrom = re.sub('^chr','',chrom)
    #column name & column number in SQL tables
    mapping = {
        'Histone_Mark': {
            'Method': 0,
            'Location': '',
            'Histone Mark': 1,
            'Cell Type': 2,
            'Additional Info': 3,
            'Reference': 4
        }
    }

    dbfile = 'RegDB.'+ str(chrom) + '.db'
    query = (pos,pos)
    query_state = "SELECT * FROM (SELECT data_index.minX, data_index.maxX, data.label_id FROM data, data_index WHERE data.id=data_index.id AND minX <= ? AND maxX > ?) AS hits, labels WHERE labels.label_id = hits.label_id"
    
    conn = sqlite3.connect(DB_DIR + dbfile)
    c = conn.cursor()
    hits = defaultdict(lambda: defaultdict(list)) # create a two dimensional dict with list as values to store all query results
    for record in c.execute(query_state,query):
        min_pos,max_pos,id1,id2,display_table= record[:5]
        fields = record[5:] #annotaions start from 5th column
        if display_table in mapping.keys(): 
            for colname, column in mapping[display_table].items():
                if colname != 'Location': #for columns in fields (location not included)
                    hits[display_table][colname].append(fields[int(column)])
    conn.close()
    return hits

###Function to get organ-specific histone mark features
def calculate_organ_sp_feature_score(chrom,pos,organs_query,histone_list,DB_DIR,biosample_organ_mapping = 'biosample_organ.txt'):
    '''
    Args:
        chrom (int or str): query chromosome (1..22, X, Y)
        pos (int): query position (0-based)
        organs_query (list of str): biosample organ name(s) (e.g. ['bodily fluid','blood','brain'])
        histone_list (list of str): histone mark names
        DB_DIR (str): directory of SQL database on Histone marks
    Returns:
        feature_values (dict of list): key: organ names in input organ list; values: 'chrom','end' & binary values of each histone mark in histone list
    '''
    from collections import defaultdict
    hits = query_RegDB_histone(chrom,pos,DB_DIR)
    feature_values = defaultdict(list)
    #get mapping from biosample to organ
    biosample_organ_dict = defaultdict(str)
    with open(biosample_organ_mapping,'r') as f:
        for line in f:
            if len(line.strip().split('\t')) < 2: #biosample does not have corresponding organ name
                next
            else:
                biosample,organ_name = line.strip().split('\t')
                biosample_organ_dict[biosample] = organ_name #organ_name could contain multiple organs, split by ','
    
    for organ_query in organs_query: #doing matching for each query organ one by one -> could have a more efficient way; save to a 2-layer dict...
        histone_mark_bin = defaultdict(int) #default value is 0
        results = [chrom,pos+1] #output chrom,end
        if 'Histone_Mark' in hits.keys():
            for histone_mark, cell in zip(hits['Histone_Mark']['Histone Mark'],hits['Histone_Mark']['Cell Type']):
                hit_organs = biosample_organ_dict[cell].split(',')
                if organ_query in hit_organs:
                    histone_mark_bin[histone_mark] = 1
        results.extend([histone_mark_bin[histone_mark] for histone_mark in histone_list]) #return in order of the histone marks in list file
        feature_values[organ_query] = results
    return feature_values
####################################################################

#Read position and genotypes from input bed file
df_input = pd.read_csv(input_bed,sep='\t',header = None)
df_input.columns = ['chrom','start','end','ref','alt']

#Retrieve RegDB generic features from input json file
RegDB_generic_features = []
with open(input_json) as f:
    for line in f:
        var_json = json.loads(line.strip())
        RegDB_generic_features.append(get_generic_features(var_json))
        
#Add RegDB generic features to dataframe
df_input = pd.merge(df_input,pd.DataFrame(RegDB_generic_features,columns = ['chrom','end','CHIP','DNASE','PWM','FOOTPRINT','EQTL_2','PWM_matched','FOOTPRINT_matched']))

#Get other precalculated generic features from bigwig files
IC_change_bw = {}
IC_matched_change_bw = {}
funsig_bw = {}
for base in ['A','C','G','T']:
    IC_change_bw[base] = pyBigWig.open(path + 'IC_change_max_' + base + '.bw')
    IC_matched_change_bw[base] = pyBigWig.open(path + 'IC_matched_change_max_' + base +
 '.bw')
    funsig_bw[base] = pyBigWig.open(path + 'funsig_' + base +
 '.bw')    
df_input['IC_change'] = [np.nan_to_num(np.round(IC_change_bw[alt].values(chrom,int(end)-1,int(end))[0],4)) for chrom,end,alt in zip(df_input['chrom'],df_input['end'],df_input['alt'])]
df_input['IC_matched_change'] = [np.nan_to_num(np.round(IC_matched_change_bw[alt].values(chrom,int(end)-1,int(end))[0],4)) for chrom,end,alt in zip(df_input['chrom'],df_input['end'],df_input['alt'])]
df_input['funsig'] = [np.nan_to_num(np.round(funsig_bw[alt].values(chrom,int(end)-1,int(end))[0],6),nan=0.127229) for chrom,end,alt in zip(df_input['chrom'],df_input['end'],df_input['alt'])] #for variants without precalculated funsig scores, use mean of dbSNP variants outside DNase peaks
for base in ['A','C','G','T']:
    IC_change_bw[base].close()
    IC_matched_change_bw[base].close()
    funsig_bw[base].close()
for ChIP_sig_feature in ['ChIP_var','ChIP_max','ChIP_quantile1','ChIP_quantile2','ChIP_quantile3']:
    bw = pyBigWig.open(path + ChIP_sig_feature + '.bw')
    df_input[ChIP_sig_feature] = [np.round(np.array(bw.values(chrom,end-1,end))[0],4) for chrom,end in zip(df_input['chrom'],df_input['end'])] 
    bw.close()
    
#predict TURF generic scores
generic_features = ['CHIP','DNASE','PWM','FOOTPRINT','EQTL_2','PWM_matched','FOOTPRINT_matched','IC_change','IC_matched_change','funsig','ChIP_var','ChIP_max','ChIP_quantile1','ChIP_quantile2','ChIP_quantile3']
clf_generic = pickle.load(open('Models/TURF_generic.sav','rb'))
df_input['TURF_generic'] = clf_generic.predict_proba(df_input[generic_features])[:,1] 

#Read query organs for TURF organ-specific scores and mapping from biosample to organ names
organs_query = []
with open(organ_list) as f0:
    for line in f0:
        organs_query.append(line.strip())

#retrieve organ-specific features for each query organ
#DNASE_organSp & FOOTPRINT_organSp from 'peaks' in json objects
RegDB_organSp_features_all = defaultdict(list)
with open(input_json) as f:
    for line in f:
        var_json = json.loads(line.strip())
        temp = get_organ_sp_RegDB_features(var_json,organs_query)
        for organ in organs_query:
            RegDB_organSp_features_all[organ].append(temp[organ])
#Organ-specific features from histone marks
histone_list = ['H3K27ac','H3K36me3','H3K4me1','H3K4me3','H3K27me3']
histone_organSp_features_all = defaultdict(list)
for _,row in df_input.iterrows():
    temp = calculate_organ_sp_feature_score(row['chrom'],row['start'],organs_query,histone_list,path + 'Database_Histone_Mark_2019/')
    for organ in organs_query:
        histone_organSp_features_all[organ].append(temp[organ])
#predict TURF organ-specific scores
df_out = df_input[['chrom','start','end','ref','alt','TURF_generic']].copy()
clf_tissueSp = pickle.load(open('Models/TURF_tissue_specific.pkl','rb'))
organSp_features = ['H3K27ac_organSp',
 'H3K36me3_organSp',
 'H3K4me1_organSp',
 'H3K4me3_organSp',
 'H3K27me3_organSp',
 'DNASE_organSp',
 'FOOTPRINT_organSp']
for organ in organs_query:
    df_1 = pd.DataFrame(RegDB_organSp_features_all[organ],columns = ['chrom','end','DNASE_organSp','FOOTPRINT_organSp'])
    df_2 = pd.DataFrame(histone_organSp_features_all[organ],columns = ['chrom','end'] + [histone + '_organSp' for histone in histone_list])
    df_organSp_features = pd.merge(df_1,df_2)
    df_organSp_features = pd.merge(df_organSp_features,df_input[['chrom','end','ref','alt','TURF_generic']])
    df_organSp_features['organSp_only'] = np.mean(np.array([list(clf_tissueSp[c].predict_proba(df_organSp_features[organSp_features])[:,1]) for c in clf_tissueSp.keys()]),axis=0)
    df_organSp_features['TURF_organSp_'+ re.sub(' ','_',organ)] = np.round(df_organSp_features['organSp_only'] * df_organSp_features['TURF_generic'],4)
    df_out = pd.merge(df_out,df_organSp_features[['chrom','end','ref','alt','TURF_generic','TURF_organSp_'+ re.sub(' ','_',organ)]])
#output all TURF predictions
df_out.to_csv(outfile,sep='\t',index=False)
