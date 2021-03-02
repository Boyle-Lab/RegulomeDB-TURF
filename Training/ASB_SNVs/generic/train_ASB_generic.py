import numpy as np
import pandas as pd
import pickle
from sklearn.ensemble import RandomForestClassifier

from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import average_precision_score, roc_auc_score

#train generic rf model
train_ASB = pd.read_csv('train_ASB_FINAL.txt',sep='\t')
features = ['CHIP','DNASE','PWM','FOOTPRINT','EQTL_2','PWM_matched','FOOTPRINT_matched','IC_change','IC_matched_change','funsig','ChIP_var','ChIP_max','ChIP_quantile1','ChIP_quantile2','ChIP_quantile3']
x=np.array(train_ASB[features])
y=np.array(train_ASB['label'])
clf=RandomForestClassifier(n_estimators=100,random_state=100,oob_score=True)
clf.fit(x,y)

#test on MPRA variants
test_MPRA_eqtl = pd.read_csv('test_MPRA_eqtl_FINAL.txt',sep='\t')
y_scores = clf.predict_proba(test_MPRA_eqtl[features])[:,1]
y_true = test_MPRA_eqtl['label']

print (np.round(average_precision_score(y_true, y_scores),3))
print (np.round(roc_auc_score(y_true, y_scores),3))
print (np.round(np.corrcoef(y_true,y_scores)[0,1],3))