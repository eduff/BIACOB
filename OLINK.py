#!/usr/bin/env python
# coding: utf-8

# # OLINK


import proteomic_funcs
import importlib
importlib.reload(proteomic_funcs)
from proteomic_funcs import *
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

import statsmodels.api as sm
from sklearn.base import BaseEstimator, RegressorMixin
import sklearn.model_selection
from  sklearn.svm import SVC
from sklearn.svm import SVR
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.ensemble import HistGradientBoostingClassifier 


#warnings.filterwarnings("ignore", category=DeprecationWarning).simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

class SMWrapper(BaseEstimator, RegressorMixin):
    """ A universal sklearn-style wrapper for statsmodels regressors """
    def __init__(self, model_class, fit_intercept=True):
        self.model_class = model_class
        self.fit_intercept = fit_intercept
    def fit(self, X, y):
        if self.fit_intercept:
            X = sm.add_constant(X)
        self.model_ = self.model_class(y, X)
        self.results_ = self.model_.fit()
        return self
    def predict(self, X):
        if self.fit_intercept:
            X = sm.add_constant(X)
        return self.results_.predict(X)


os.chdir('/home/eduff/biobank/Proteomics/OLINK')


# OLINK_data=pd.read_csv('olink_data2.txt',sep='\t')

# eids_OLINK=OLINK_data[OLINK_data['ins_index']==2]['eid'].unique()
# eids_SIMOA=data.index.unique()
# common_eids = np.union1d(eids_OLINK,data.index)
# np.savetxt('subjects_OLINK.txt',common_eids),fmt='%10.0f')

# OLINK=pd.DataFrame(index=np.unique(OLINK_data['eid']))
# Plists={}
# Plists[0]=[]
# Plists[2]=[]
# Plists[3]=[]

# Plists['all'] = OLINK_data['protein_id'].unique()

# for a in Plists['all']:
#     for b in [0,2,3]:
#         Plists[b].append('P_' + str(a) + '-' + str(b))
#         els=(OLINK_data['protein_id']==a)&(OLINK_data['ins_index']==b)
#         outsuff='0'
#         if b==2:
#             outsuff='pre'
#         elif b==3:
#             outsuff='post'
#         OLINK.loc[OLINK_data.loc[els,'eid'],'P_' + str(a) + '_' + outsuff] = OLINK_data.loc[els,'result'].values
        
OLINK_meta=pd.read_csv('ukb672393.csv')
OLINK_meta=OLINK_meta.set_index('eid')
OLINK_meta=OLINK_meta.loc[data.index,:]


# Plists['all'] = list(OLINK_data['protein_id'].unique())
# #  participants with too many missing vals
# parts_t0=(OLINK.loc[:,Plists[0]].notnull().sum(axis=1)>500)
# parts_t23=(OLINK.loc[:,Plists[2]].notnull().sum(axis=1)>500) & (OLINK.loc[:,Plists[3]].notnull().sum(axis=1)>500)

# common_eids=np.array([ a in data.index for a in OLINK.index])
# remain_eids = (~common_eids & parts_t0)
# common_eids = (common_eids & parts_t23 & parts_t0)
# common_eids_0 = (common_eids & parts_t23 )


# coding = pd.read_csv('coding143.tsv',sep='\t',index_col=0)
# coding.loc[:,'Protein']=np.nan
# coding.loc[:,'Protein']=coding['meaning'].str.extract('(.*);').values

# OLINK.columns=[a.replace('-post','_post') for a in OLINK.columns]
# OLINK.columns=[a.replace('-pre','_pre') for a in OLINK.columns]

# repl={}
# for a in coding.index:
#     cols=OLINK.columns[OLINK.columns.str.contains('_'+str(a)+'-')]
    
#     for col in cols:
#         repl[col]=col.replace('_'+str(a)+'-','_'+coding.loc[a,'Protein']+'_')
# OLINK=OLINK.rename(columns=repl)
# OLINK.columns=[a.replace('_2','_pre') for a in OLINK.columns]
# OLINK.columns=[a.replace('_3','_post') for a in OLINK.columns]
# OLINK.to_csv('OLINK_data.csv')



# OLINK=pd.read_csv('/home/eduff/biobank/Proteomics/OLINK/OLINK_data.csv',index_col=0)
# OLINK.columns=[a.replace('_3','_post') for a in OLINK.columns]
# OLINK.columns=[a.replace('_2','_pre') for a in OLINK.columns]
# cols=(OLINK.loc[:,OLINK.columns.str.contains('_pre')].notna().sum()>20)
# cols=list(cols[cols].index.values)

# all_OLINK=(OLINK[cols].notna().sum(axis=1))>0
# OLINK_matched=(data.loc[:,'matched']==True) & all_OLINK)
# eids_OLINK=eids_OLINK[eids_OLINK].index.values

# common_eids = np.union1d(eids_OLINK,data.index)

# #eids_OLINK=OLINK_data[OLINK_data['ins_index']==2]['eid'].unique()


OLINK=pd.read_csv('/home/eduff/biobank/Proteomics/OLINK/OLINK_data_COVID.csv',index_col=0)

eids_OLINK=OLINK.index.values # (OLINK[cols].notna().sum(axis=1))>0


# estimate number of available samples

data.loc[:,'OLINK']=False
data.loc[eids_OLINK,'OLINK']=True
((data.loc[eids_OLINK,'matched'] & (data.loc[eids_OLINK,'matched_eid'].notna())).sum())


((data.loc[eids_OLINK,'matched'] & (data.loc[eids_OLINK,'matched_eid'].notna())&all_case&all_matched).sum())
data.loc[eids_OLINK,OLINK.columns]=OLINK.loc[eids_OLINK,:]


SIMOA_assays = ['Ab40_regPl', 'Ab42_regPl', 'Ab42/Ab40_regPl', 'pTau-181_regPl', 'NfL_regPl', 'GFAP_regPl']
OLINK_assays_red=['P_NEFL','P_GFAP' ,'P_TREM2','P_AXL','P_TYRO3','P_MIF','P_C1QA', 'P_C1QTNF1','P_C1QL2','P_C4BPB','P_IL6','P_IL18']
OLINK_assays=['P_TREM2','P_AXL','P_TYRO3','P_MIF','P_IL6','P_IL18','P_VWF','P_TNF','P_IL1A','P_IL1B','P_IFNG','P_PTX3','P_CALCA','P_GZMB','P_GDNF','P_NGF','P_GDF15', 'P_TNFSF10','P_TNFRSF10A','P_TNFRSF10B', 'P_BCAN', 'P_GFAP', 'P_NEFL','P_PIGR','P_CHI3L1','P_L1CAM','P_NCAM1','Alzheimers_dementia','Vascular_dementia']
OLINK_assays_final=['P_IL1A','P_IL1B','P_IL6','P_IL18','P_IFNG','P_TNF','P_TNFSF10','P_C1QA', 'P_C4BPB','P_CHI3L1','P_MIF','P_PTX3','P_CALCA',]# 'P_TINAGL1','P_GDF15','P_ICA1','P_CDON','P_BCAN','P_ADGRG1','P_INPPL1','P_EDA2R','P_KDR','P_ADAMTS8','P_P4HB','P_ODAM','P_ATF2']
'P_C1QA',("'P_C4BPB'")
#extra_vars=['P_TREM2','P_AXL','P_TYRO3','P_MIF','P_C1QA', 'P_C4BPB','P_IL6','P_IL18','P_VWF','P_FGA','P_NRGN','P_C5','P_C8','P_TNF','P_IL12','P_IL1A','P_IL1B','P_IFNG','P_CSF2','P_OCLN','P_PTX3','P_IFNL2','P_CALCA','P_GZMB','P_C1QTNF9','P_CFH','P_APOE','P_BDNF','P_GDNF','P_NGF']

OLINK_assays_all=[a[:-4] for a in find_vars(data,'^P_.*pre$')]



OLINK_meta.loc[:,find_vars(OLINK_meta,'30901')]-8.9e+11


plt.hist(OLINK_meta.loc[:,find_vars(OLINK_meta,'30901')[0]]-8.9e+11)


# 

prots_avail=(data[(find_vars(data,'^P_.*pre$'))].notna().sum(axis=0)>500)


# apply neuroimaging cleaning to OLINK proteomics 
CONF=['Age-3.0_f','assessment_sep','assessment_sep^2','Ethnicity(White)','22001-0.0']
for a in OLINK_assays_all:
    out = remove_conf(data,a,[],flatten=True,suf_pre=True,remove_out=True,suf=False,int=True)
    data.loc[:,a+'_pre_cl'] = out[a+'_pre_cl']
    data.loc[:,a+'_post_cl'] = out[a+'_post_cl']
    data.loc[:,a+'_0_cl'] = remove_conf(data,a+'_0',[],flatten=False,suf_pre=True,remove_out=True,suf=False,int=True)

    data.loc[:,a+'_diff_cl'] = data[a+'_post_cl']-data[a+'_pre_cl']
    data.loc[:,a+'_pc_cl'] = 100*data[a+'_diff_cl']/data[a+'_pre_cl']
    


disease_weights= pd.read_csv('/home/eduff/biobank/Proteomics/OLINK/supp_13_disease_weights.csv',header=3)
disease_weights.loc[:,'Protein']=disease_weights['Protein'].str.extract('([^\.]*)\.').values
disease_weights.loc[177,'Protein']='HLA-E'

for a in disease_weights['ProteinScore'].unique():
    astr=a.replace(' ','_')
    astr=astr.replace("'","")
    weights=disease_weights.loc[disease_weights['ProteinScore'].str.contains(a),['Protein','Coefficient']]
    for b in ['0','pre_cl','post_cl']:
        score=(data.loc[eids_OLINK,'P_'+weights['Protein'].values+'_'+b]*weights.loc[:,'Coefficient'].values).sum(axis=1)
        data.loc[score.index,astr+'_'+b]=score.values
        #tmp=score.values
    data.loc[score.index,astr+'_'+'diff_cl']=data.loc[score.index,astr+'_post_cl']-data.loc[score.index,astr+'_pre_cl']
diseases=[a.replace(' ','_') for a in disease_weights['ProteinScore'].unique()]        
diseases_orig=diseases
diseases[9]="Alzheimers_dementia"
diseases[12]="Parkinsons_disease"

data.loc[eids_OLINK,'OLINK']=True

#for a in [
#data=data.rename(columns={"Alzheimer's_dementia":"Alzheimers_dementia","Parkinson's_disease":'Parkinsons_disease'})


a='Alz'
b='pre_cl'
weights=disease_weights.loc[disease_weights['ProteinScore'].str.contains(a),['Protein','Coefficient']]
(data.loc[eids_OLINK,'P_'+weights['Protein'].values+'_'+b]*weights.loc[:,'Coefficient'].values).var()
'P_GDF15', 'P_TNFSF10', 'P_BCAN', 'P_GFAP', 'P_NEFL'


[a+'_pre_cl' for a in OLINK_assays_all]+['Ab40_pre']


find_vars(data,'P_EIF4G1')


compare_prots[-5:]


# check
compare_prots=['P_MMP9_pre','P_MMP8_pre','NfL_post','P_NEFL_post','P_PTX3_post','pTau-181_pre','pTau-181_post','Ab40_pre','Ab40_post','Ab42_pre','Ab42_post']
# 'GFAP_pre','P_GFAP_pre','GFAP_post','P_GFAP_post','NfL_pre',
compare_prots=[a+'_pre_cl' for a in OLINK_assays_all]+['Ab40_regPl_pre_cl','Ab42_regPl_pre_cl']
fig, ax = plt.subplots(figsize=(10,10))
corr_map_OLINK = data.loc[:,compare_prots].corr()
#mask_m = np.tril(corr_map_OLINK,k=-1)
#ax=sns.heatmap(corr_map_OLINK,cmap="coolwarm",vmin=-.8,vmax=.8, annot=True,ax=ax)
#ax.invert_yaxis()


corr_map_OLINK.loc['Ab40_regPl_pre_cl',:].sort_values()


# corr of OLINK with Ab42/Ab40
Abcorrs=data[['Ab42/Ab40_regPl_diff_cl']+[a+'_diff_cl' for a in OLINK_assays_all]].corr()['Ab42/Ab40_regPl_diff_cl']
Abcorrs.sort_values().iloc[-20:]


for a in ['NTRK3', 'NTRK2', 'BLMH', 'CBLN4', 'PTPRN2', 'PTPRS']:
    print(len(data.loc[:,'P_'+a+'_pre_cl'].dropna()))


aa=data.loc[:,'P_MCFD2_post'].dropna().index.intersection(data.loc[:,'Ab42_post'].dropna().index)
bb=data.loc[:,'P_MCFD2_pre'].dropna().index.intersection(data.loc[:,'Ab42_pre'].dropna().index)
aa.intersection(bb)


prot_els_regPl_pre_cl
Abcorrs=data[prot_els_regPl_pre_cl +[a+'_pre_cl' for a in diseases]].corr()[prot_els_regPl_pre_cl ]
Abcorrs['Ab42_regPl_pre_cl'].sort_values()


