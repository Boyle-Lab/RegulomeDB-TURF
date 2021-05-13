import pickle
import numpy as np
from sklearn.ensemble import RandomForestClassifier
import pandas as pd

#Train TURF generic rf model
train_ASB = pd.read_csv('../Training/ASB_SNVs/generic/train_ASB_FINAL.txt',sep='\t') #Load generic feature dataframe for all ASB SNVs in 6 cell lines
features = ['CHIP','DNASE','PWM','FOOTPRINT','EQTL_2','PWM_matched','FOOTPRINT_matched','IC_change','IC_matched_change','funsig','ChIP_var','ChIP_max','ChIP_quantile1','ChIP_quantile2','ChIP_quantile3']
x=np.array(train_ASB[features])
y=np.array(train_ASB['label'])
clf_generic=RandomForestClassifier(n_estimators=100,random_state=100,oob_score=True)
clf_generic.fit(x,y)
pickle.dump(clf_generic,open('TURF_generic_2.sav','wb')) 

#Train tissue-specific rf models for each ASB cell line with the 7 features in the final TURF ensemble model
features_names_cell_sp = ['H3K27ac_cellSp',
 'H3K36me3_cellSp',
 'H3K4me1_cellSp',
 'H3K4me3_cellSp',
 'H3K27me3_cellSp',
 'DNASE_cellSp',
 'FOOTPRINT_cellSp']
cells = ['GM12878','A549','H1hESC','HepG2','K562','MCF7']
df_all_dict = pickle.load(open('../Training/ASB_SNVs/tissue_specific/df_cell_sp_all_FINAL.pkl','rb')) #Load tissue-specific feature dataframes for each ASB cell line
clf = {}
for cell in cells:
    df_all = df_all_dict[cell]
    x = np.array(df_all[features_names_cell_sp])
    y = df_all['cellSp_label']
    clf[cell]=RandomForestClassifier(class_weight='balanced',n_estimators=100,random_state=100,oob_score=True)
    clf[cell].fit(x,y)
pickle.dump(clf,open('TURF_tissue_specific.pkl','wb')) 
