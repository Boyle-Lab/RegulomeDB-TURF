import numpy as np
import pandas as pd
import pickle

#load random forest models
clf_generic = pickle.load(open('TURF_generic.sav','rb'))
clf_tissueSp = pickle.load(open('TURF_tissue_specific.pkl','rb'))
#feature lists
generic_features = ['CHIP','DNASE','PWM','FOOTPRINT','EQTL_2','PWM_matched','FOOTPRINT_matched','IC_change','IC_matched_change','funsig','ChIP_var','ChIP_max','ChIP_quantile1','ChIP_quantile2','ChIP_quantile3']
tissueSp_features = ['H3K27ac_tissueSp',
 'H3K36me3_tissueSp',
 'H3K4me1_tissueSp',
 'H3K4me3_tissueSp',
 'H3K27me3_tissueSp',
 'DNASE_tissueSp',
 'FOOTPRINT_tissueSp']
#read input variants with features assigned
df_input = pd.read_csv('example_input_features.txt',sep='\t')
#predict with TURF generic and tissue-specific scores
x_generic = df_input[generic_features]
x_tissueSp = df_input[tissueSp_features]
y_generic = clf_generic.predict_proba(x_generic)[:,1] 
y_tissueSp_only = np.mean(np.array([list(clf_tissueSp[c].predict_proba(x_tissueSp)[:,1]) for c in clf_tissueSp.keys()]),axis=0)
y_tissueSp_score = y_generic * y_tissueSp_only
df_input['TURF_generic'] = y_generic
df_input['TURF_tissueSp'] = np.round(y_tissueSp_score,4)
#output TURF scores
df_out = df_input[['chrom','end','ref','alt','TURF_generic','TURF_tissueSp']]
df_out.to_csv('example_output_scores.txt',index=False,sep='\t')
