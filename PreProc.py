#!/usr/bin/env python
# coding: utf-8

import proteomic_funcs
import importlib
importlib.reload(proteomic_funcs)
from proteomic_funcs import *
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)


# # Preprocessing

# ## Reading in Proteomics
# 
# Proteomics stored in .csv format from assay acquisition, separating Plex4 and Ptau

os.chdir('/home/eduff/biobank/Proteomics/process')


# Reading original pre-release raw data, may require different processing from Biobank

Plex4 = pd.read_csv('4plexE_all_plates_eid.csv')
Plex4.loc[Plex4['spike']=='FALSE ','spike']=False
Plex4.loc[Plex4['spike']=='TRUE ','spike']=True

ptau = pd.read_csv('pTau-181_all_plates_eid.csv')
ptau.loc[ptau['spike']=='FALSE ','spike']=False
ptau.loc[ptau['spike']=='TRUE ','spike']=True


Plex4.loc[Plex4['SampleID'].str[-1:]=='1','Session']='1'
Plex4.loc[Plex4['SampleID'].str[-1:]=='a','Session']='a'

ptau.loc[ptau['SampleID'].str[-1:]=='1','Session']='1'
ptau.loc[ptau['SampleID'].str[-1:]=='a','Session']='a'

Plex4_spike=Plex4[Plex4['spike']==True]
Plex4=Plex4[Plex4['spike']==False]

ptau_spike=ptau[ptau['spike']==True]
ptau=ptau[ptau['spike']==False]


# remove singletons
for a in Plex4['eid']:
    if len(Plex4[(Plex4['eid']==a) & (Plex4['spike']==False)])!=2:
        print(a)
        Plex4=Plex4[Plex4['eid']!=a]
        ptau=ptau[ptau['eid']!=a]
for a in ptau['eid']:
    if len(ptau[(ptau['eid']==a) & (ptau['spike']==False)])!=2:
        Plex4=Plex4[Plex4['eid']!=a]
        ptau=ptau[ptau['eid']!=a]
        
# remove potentially mislabelled samples
for a in [339983033,339983037]:
    Plex4=Plex4[Plex4['VacutainerID']!=a]
    ptau=ptau[ptau['VacutainerID']!=a]
    


# combine proteins
all_prot=pd.concat([ptau,Plex4],axis=1)
all_prot = all_prot.loc[:,~all_prot.columns.duplicated()]

# rename proteins
all_prot =all_prot.rename(columns={'Neurology 4plexE AB42 pg/ml Replicate 1':'Ab42','Neurology 4plexE AB40 pg/ml Replicate 1':'Ab40','Neurology 4plexE NFL pg/ml Replicate 1' : 'NfL' , 'Neurology 4plexE GFAP pg/ml Replicate 1' : 'GFAP', 'pTau-181 pg/ml Replicate 1' : 'pTau-181' })


# ### Preprocess Proteomics: log, regress plateID

# Non A-beta assays have skewed distributions - taking logs 

for a in ['GFAP','NfL','pTau-181']:
    all_prot.loc[:,a+'_orig'] = all_prot.loc[:,a]
    all_prot.loc[:,a] = np.log(all_prot.loc[:,a] )

# remove outliers 
for a in ['Ab42','Ab40','GFAP','NfL','pTau-181']:
    all_prot.loc[:,a]=remove_outliers(all_prot.loc[:,a])
   
for  a in ['GFAP','NfL','pTau-181']:
    all_prot.loc[all_prot.loc[:,a]==np.nan,a+'orig']=np.nan   
    
    
# regress plate ID
for a in ['Ab42','Ab40','GFAP','NfL','pTau-181']:
    mod=smf.ols(formula="  (Q('"+a+"'))  ~  C(PlateID)   ", data=all_prot).fit()
    all_prot.loc[:,a+'_regPl']=all_prot.loc[:,a]-(mod.fittedvalues-np.mean(mod.fittedvalues))
    all_prot.loc[:,a+'_orig_regPl']=np.exp(all_prot.loc[:,a+'_regPl'])
                 
#for a in ['Ab42','Ab40','GFAP','NfL','pTau-181']
# ratio 
all_prot.loc[:,'Ab42/Ab40']=all_prot['Ab42']/all_prot['Ab40']
all_prot.loc[:,'Ab42/Ab40_regPl']=all_prot['Ab42_regPl']/all_prot['Ab40_regPl']
all_prot_pre=(all_prot[all_prot['Session']=='a']).sort_values('eid').add_suffix('_pre')
all_prot_post=(all_prot[all_prot['Session']=='1']).sort_values('eid').add_suffix('_post')
all_prot_pre=all_prot_pre.set_index('eid_pre')
all_prot_post=all_prot_post.set_index('eid_post')
proteomics=pd.concat([all_prot_pre,all_prot_post],axis=1)


# bunch of helper variables
prot_els=['Ab40','Ab42','Ab42/Ab40','pTau-181','NfL','GFAP']
prot_els_pre =  [ x+'_pre' for x in prot_els]
prot_els_post = [ x+'_post' for x in prot_els]
prot_els_all = prot_els_pre + prot_els_post
prot_els_regPl =  [ x+'_regPl' for x in prot_els]
prot_els_regPl_pre =  [ x+'_pre' for x in prot_els_regPl]
prot_els_regPl_post = [ x+'_post' for x in prot_els_regPl]
prot_els_regPl_pre_cl =  [ x+'_pre_cl' for x in prot_els_regPl]
prot_els_regPl_post_cl = [ x+'_post_cl' for x in prot_els_regPl]
prot_els_all_regPl = prot_els_regPl_pre + prot_els_regPl_post
#prot_els_regPl = prot_els_regPl_pre + prot_els_regPl_post


# ## Read in Phenotypic , Genetic Data

# load phenotypic data
# merged_subjsO.csv is a merged file of all subjects with proteomics data across all baskets for project.
pheno= pd.read_csv('../../phenotype/merged_subjsO.csv',low_memory=False)
pheno = pheno.set_index('eid')
# bbd translates from biobank variables to more readable names
bbd={'AlzPRS':'26206','BP_sys':'4080','BP_dia':'4079','Activity':'22034','Weight':'21002','Alcohol':'1558','BMI':'21001',"gSex":'22001-0.0'}


## Load in Genetics
# focusing on APOE for this study.
genetics = pd.read_csv('../../genotyping/APOE_genotyping_ED.raw',delimiter='\t')
genetics=genetics.rename(columns={'FID':'eid'})
genetics=genetics.drop(columns={'PAT','MAT','PHENOTYPE','IID'})
genetics = genetics.set_index('eid')

# Some additional Genetic Variants for future analysis
variants = pd.read_csv('../../genotyping/Genes.csv')
variants = variants.rename(columns={"IID":"eid"})
variants=variants.drop('SEX',axis=1)
variants = variants.set_index('eid')
#variants=variants.drop('SEX',axis=1)
variants_cols = list(variants.columns)
LRRK=pd.read_csv('/home/eduff/biobank/genotyping/LRRK_genotyping_ED2.raw',delimiter='\t')
LRRK=LRRK.set_index('FID')
els=variants.index
variants.loc[els,'LRRK']=LRRK.loc[els,'rs76904798_T']
variants.iloc[:,2:].columns



AD_variants=['APOE_score']+list(variants.iloc[:,2:].columns)


# tidying up APOE for analysis

genetics.loc[:,'APOE']=np.NaN
genetics.loc[(genetics['rs7412_T']==0) & (genetics['rs429358_C']==0),'APOE']='A3A3'
genetics.loc[(genetics['rs7412_T']==0) & (genetics['rs429358_C']==1),'APOE']='A3A4'
genetics.loc[(genetics['rs7412_T']==0) & (genetics['rs429358_C']==2),'APOE']='A4A4'
genetics.loc[(genetics['rs7412_T']==1) & (genetics['rs429358_C']==0),'APOE']='A3A2'
genetics.loc[(genetics['rs7412_T']==2) & (genetics['rs429358_C']==0),'APOE']='A2A2'

# create an APOE var without 'na'
genetics.loc[:,'APOE_missing']=genetics.loc[:,'APOE']
genetics.loc[genetics.loc[:,'APOE'].isna(),'APOE_missing']='missing'

# counts of A4
genetics.loc[:,'A4']=np.NaN
genetics.loc[genetics.loc[:,'APOE']=='A4A4','A4']=2
genetics.loc[genetics.loc[:,'APOE']=='A3A4','A4']=1
genetics.loc[genetics.loc[:,'APOE']=='A3A3','A4']=0
genetics.loc[genetics.loc[:,'APOE']=='A3A2','A4']=0
genetics.loc[genetics.loc[:,'APOE']=='A2A2','A4']=0
# counts of A2
genetics.loc[:,'A2']=np.NaN
genetics.loc[genetics.loc[:,'APOE']=='A4A4','A2']=0
genetics.loc[genetics.loc[:,'APOE']=='A3A4','A2']=0
genetics.loc[genetics.loc[:,'APOE']=='A3A3','A2']=0
genetics.loc[genetics.loc[:,'APOE']=='A3A2','A2']=1
genetics.loc[genetics.loc[:,'APOE']=='A2A2','A2']=2
genetics.loc[:,'APOE_score']=genetics.loc[:,'A4']-genetics.loc[:,'A2']
genetics.loc[:,'A4_bin']=genetics.loc[:,'A4']>0
genetics.loc[:,'A44_bin']=genetics.loc[:,'A4']>1
genetics.loc[:,'A33vA34']=0
genetics.loc[genetics.loc[:,'APOE']=='A3A4','A33vA34']=0.5
genetics.loc[genetics.loc[:,'APOE']=='A3A3','A33vA34']=-0.5
genetics.loc[:,'A33vA44']=0
genetics.loc[genetics.loc[:,'APOE']=='A4A4','A33vA44']=0.5
genetics.loc[genetics.loc[:,'APOE']=='A3A3','A33vA44']=-0.5

genetics['A34']=genetics.loc[:,'APOE']=='A3A4'
genetics.loc[genetics.loc[:,'APOE'].isna(),'A34']=np.nan
genetics['A34']=genetics['A34'].convert_dtypes()
genetics['A44']=genetics.loc[:,'APOE']=='A4A4'
genetics.loc[genetics.loc[:,'APOE'].isna(),'A44']=np.nan
genetics['A44']=genetics['A44'].convert_dtypes()
genetics['A32']=genetics.loc[:,'APOE']=='A3A2'
genetics.loc[genetics.loc[:,'APOE'].isna(),'A32']=np.nan
genetics['A32']=genetics['A32'].convert_dtypes()
genetics['A33']=genetics.loc[:,'APOE']=='A3A3'
genetics.loc[genetics.loc[:,'APOE'].isna(),'A33']=np.nan
genetics['A33']=genetics['A33'].convert_dtypes()


#  APOE contrasts scores
genetics.loc[:,'A33vA34']=np.nan
genetics.loc[genetics.loc[:,'APOE']=='A3A4','A33vA34']=0.5
genetics.loc[genetics.loc[:,'APOE']=='A3A3','A33vA34']=-0.5
genetics['A33vA34']=(genetics['A33vA34']-genetics['A33vA34'].mean())
genetics.loc[genetics['A33vA34'].isna(),'A33vA34']=0
genetics.loc[:,'A33vA44']=np.nan
genetics.loc[genetics.loc[:,'APOE']=='A4A4','A33vA44']=0.5
genetics.loc[genetics.loc[:,'APOE']=='A3A3','A33vA44']=-0.5
genetics['A33vA44']=genetics['A33vA44']-genetics['A33vA44'].mean()
genetics.loc[genetics['A33vA44'].isna(),'A33vA44']=0
genetics.loc[:,'A33vA32']=np.nan
genetics.loc[genetics.loc[:,'APOE']=='A3A2','A33vA32']=0.5
genetics.loc[genetics.loc[:,'APOE']=='A3A3','A33vA32']=-0.5
genetics['A33vA32']=genetics['A33vA32']-genetics['A33vA32'].mean()
genetics.loc[genetics['A33vA32'].isna(),'A33vA32']=0
genetics.loc[:,'A33vA32']=np.nan
genetics.loc[genetics.loc[:,'APOE']=='A3A2','A33vA32']=0.5
genetics.loc[genetics.loc[:,'APOE']=='A3A3','A33vA32']=-0.5
genetics['A33vA32']=genetics['A33vA32']-genetics['A33vA32'].mean()
genetics.loc[genetics['A33vA32'].isna(),'A33vA32']=0


# merge proteomics, pheno, genetics
data= pd.merge(proteomics, pheno, how='right',left_index=True, right_index=True)
data= pd.merge(data, variants,how='left', left_index=True, right_index=True,suffixes=("_x",None))
data= pd.merge(data, genetics,how='left', left_index=True, right_index=True,suffixes=("_x",None))
data=data.copy()


PRS=pd.read_csv('/home/eduff/biobank/genotyping/Enhanced_PRS.csv',index_col=0,header=None) 


# ## Removing people withdrawn from UK Biobank 

withdrawals=pd.read_csv('../../withdrawal-76059_2023-04-25.csv',header=None)
np.array([ a in data.index for a in withdrawals[0]]).sum()


# ## Additional data and helper columns from biobank 

# find_vars fine case / control variable
data_case=data.loc[(np.mod(data.loc[:,'41000-3.0'],10)==1) ,['41000-3.0']]
data_case=data_case-1
data_case=data_case.reset_index()
data_case=data_case.rename(columns={'eid':'case'})
data_control=data.loc[(np.mod(data.loc[:,'41000-3.0'],10)==0) ,['41000-3.0']]
data_control=data_control.reset_index()
data_control=data_control.rename(columns={'eid':'control'})

data.loc[:,'Case']='sars'
data.loc[data_control['control'],'Case']='ctr'

# numeric versions of Case/Ctr
data.loc[data.loc[:,'Case']=='sars','Case_bin']=1
data.loc[data.loc[:,'Case']=='ctr','Case_bin']=0

# Identify matched case/controls

cc_els=pd.merge(data_case,data_control,on='41000-3.0')
data.loc[:,'matched']=False
data.loc[np.concatenate([cc_els['case'].values , cc_els['control'].values]),'matched']=True

#SIMOA_matched=(data.loc[:,'matched']==True) &  )

all_case=data.loc[:,'Case']=='sars'
all_control=data.loc[:,'Case']=='ctr'

# matched excluding missing data for each assay

matched={}

for a in ['Ab42/Ab40','pTau-181','NfL','GFAP','Ab42','Ab40']:
    case_present=data.loc[cc_els['case'],[a+'_pre',a+'_post']].notnull().all(axis=1).values
    control_present=data.loc[cc_els['control'],[a+'_pre',a+'_post']].notnull().all(axis=1).values
    matched[a]=(case_present & control_present)
    # print(a + ' ' + str(np.sum(matched[a])))
    data.loc[:,a+'_matched']=False
    data.loc[pd.concat([cc_els.loc[matched[a],'case'],cc_els.loc[matched[a],'control']]),a+'_matched']=True

data['all_prot_matched']=data.loc[:,[ a + '_matched' for a in ['Ab42/Ab40','pTau-181','NfL','GFAP']]].all(axis=1)
all_matched=data.loc[:,[ a + '_matched' for a in ['Ab42/Ab40','pTau-181','NfL','GFAP']]].any(axis=1)
matched_case=all_matched & all_case

# matched case/controls available
data.loc[cc_els['case'],'matched_eid']=cc_els['control'].values.astype(int)
data.loc[cc_els['control'],'matched_eid']=cc_els['case'].values.astype(int)

groupby = 'Case'
group='sars'
dd=data.loc[all_matched,['matched_eid','Gender_pre','Case']]
ids_case = dd[dd[groupby]==group].index
ids_ctr = dd.loc[ids_case,'matched_eid']


# Date information

data.loc[:,'date_post']=pd.to_datetime(data.loc[:,'53-3.0'])
data.loc[:,'date_pre']=pd.to_datetime(data.loc[:,'53-2.0']) 
data.loc[:,'date_1']=pd.to_datetime(data.loc[:,'53-1.0']) 
data.loc[:,'date_0']=pd.to_datetime(data.loc[:,'53-0.0'])                                    
data.loc[:,'assessment_sep_d']=(data.loc[:,'date_post'] - data.loc[:,'date_pre'])#/ np.timedelta64(1, 'M')
data.loc[:,'assessment_sep']=data.loc[:,'assessment_sep_d'].dt.days
data.loc[:,'assessment_sep_m']=((data.loc[:,'date_post'] - data.loc[:,'date_pre'])/ np.timedelta64(1, 'M')).astype(int)
data.loc[:,'assessment_sep^2']=data.loc[:,'assessment_sep']**2

data.loc[:,'DOB'] = pd.to_datetime(data.loc[:,'52-0.0'].astype('int').astype('str')+'/'+data.loc[:,'34-0.0'].astype('int').astype('str'))
data.loc[:,'Age-3.0_d']=(data.loc[:,'date_post'] - data.loc[:,'DOB'])
data.loc[:,'Age-2.0_d']=(data.loc[:,'date_pre'] - data.loc[:,'DOB'])
data.loc[:,'Age-1.0_d']=(data.loc[:,'date_1'] - data.loc[:,'DOB'])
data.loc[:,'Age-0.0_d']=(data.loc[:,'date_0'] - data.loc[:,'DOB'])
data.loc[:,'Age-3.0']=data.loc[:,'Age-3.0_d'].dt.days/365
data.loc[:,'Age-2.0']=data.loc[:,'Age-2.0_d'].dt.days/365
data.loc[:,'Age-1.0']=data.loc[:,'Age-1.0_d'].dt.days/365
data.loc[:,'Age-0.0']=data.loc[:,'Age-0.0_d'].dt.days/365
data.loc[:,'Age-3.0^2']=(data.loc[:,'Age-3.0']-50)**2
data.loc[:,'Age-2.0^2']=data.loc[:,'Age-2.0']**2

# Age-dependent vulnerability  regressor
data.loc[:,'Age-2.0_fc']=10**(data.loc[:,'Age-2.0']*0.0524-3.27)*data['Case_bin']/(10**(data.loc[:,'Age-3.0']*0.0524-3.27)).max()
data.loc[:,'Age-3.0_fc']=10**(data.loc[:,'Age-3.0']*0.0524-3.27)*data['Case_bin']/(10**(data.loc[:,'Age-3.0']*0.0524-3.27)).max()

# Mean age-matched
data.loc[data.loc[cc_els['case'],'matched_eid'],'matched_age_mean']= 0.5*(data.loc[data.loc[cc_els['case'],'matched_eid'],'Age-3.0'].values + data.loc[cc_els['case'],'Age-3.0'].values)
data.loc[cc_els['case'],'matched_age_mean']= 0.5*(data.loc[data.loc[cc_els['case'],'matched_eid'],'Age-3.0'].values + data.loc[cc_els['case'],'Age-3.0'].values)

# COVID vulnerability modulated regressor
dm=10**(data.loc[:,'Age-2.0']*0.0524-3.27)
data.loc[:,'Age-2.0_f']=dm/np.mean(dm)
dm=10**(data.loc[:,'Age-3.0']*0.0524-3.27)
data.loc[:,'Age-3.0_f']=dm/np.mean(dm)

# Health data

data.loc[:,'Hip/Waist-2.0']=data.loc[:,'49-2.0'] / data.loc[:,'48-2.0']
data.loc[:,'Hip/Waist-3.0']=data.loc[:,'49-3.0'] / data.loc[:,'48-3.0']
data.loc[:,'Smoking-2.0']=data.loc[:,'1239-2.0']
data.loc[:,'Smoking-3.0']=data.loc[:,'1239-3.0']

data['Smoking_bin-2.0']=data['Smoking-2.0']>0
data['Smoking_bin-3.0']=data['Smoking-3.0']>0

data.loc[data.loc[:,'Smoking-2.0']==-3,'Smoking-2.0']=np.nan
data.loc[data.loc[:,'Smoking-3.0']==-3,'Smoking-3.0']=np.nan

data.loc[:,'Alcohol-2.0']=data['1558-2.0']
data.loc[data.loc[:,'Alcohol-2.0']<0,'Alcohol-2.0']=np.nan
data.loc[:,'Alcohol-2.0']=   6-data.loc[:,'Alcohol-2.0']
data.loc[:,'Alcohol-3.0']=data['1558-3.0']
data.loc[data.loc[:,'Alcohol-3.0']<0,'Alcohol-3.0']=np.nan
data.loc[:,'Alcohol-3.0']=   6-data.loc[:,'Alcohol-3.0']
data.loc[:,'Obesity-2.0'] = data['21001-2.0']>30
data.loc[:,'Obesity-3.0'] = data['21001-3.0']>30                                             
                                                 
data.loc[:,'Deprivation'] = data.loc[:,'26410-0.0'].fillna(data.loc[:,'26427-0.0']).fillna(data.loc[:,'26426-0.0'])

data.loc[:,'Ethnicity(White)']=(np.mod(data.loc[:,'21000-0.0'],1000)==1.0).astype(int)


data.loc[:,'KeyWorker'] = ~data.loc[:,'28063-0.0'].isnull()
data.loc[:,'WorkingThroughCOVID'] = ((data.loc[:,'28057-0.0'] == 1) | (data.loc[:,'28057-0.0'] == 1) | (data.loc[:,'28057-0.0'] == 8) )
data.loc[:,'Diabetes'] = np.logical_or.reduce(data.loc[:,find_vars(data,'^2443-')]==1,axis=1)
# Cholestoral and BP 
data.loc[:,'Chol_meds']=np.logical_or.reduce(data.loc[:,find_vars(data,'^6153-')]==1,axis=1)
data.loc[:,'BP_meds']=np.logical_or.reduce(data.loc[:,find_vars(data,'^6153-')]==2,axis=1)
data.loc[:,'Insulin_meds']=np.logical_or.reduce(data.loc[:,find_vars(data,'^6153-')]==3,axis=1)
data.loc[:,'Heart_Cond']=np.logical_or.reduce(data.loc[:,find_vars(data,'^6150-')]>0,axis=1)
data.loc[:,'COPD']=np.logical_or.reduce(data.loc[:,find_vars(data,'^22130-')]>0,axis=1)


data['Wheeze-2.0']=data['2316-2.0']==1
data.loc[data['2316-2.0']<0,'Wheeze-2.0']=np.nan
data['Wheeze-2.0']=data['Wheeze-2.0'].convert_dtypes()
data['Wheeze_bin-2.0']=data['Wheeze-2.0'].astype(float)
data['Wheeze-3.0']=data['2316-3.0']==1
data.loc[data['2316-3.0']<0,'Wheeze-3.0']=np.nan
data['Wheeze_bin-3.0']=data['Wheeze-3.0'].astype(float)

data['HandGrip-2.0']=(data['46-2.0']+data['47-2.0'])/2
data['HandGrip-3.0']=(data['46-3.0']+data['47-3.0'])/2

data['GeneralHealth-2.0']=5-data['2178-2.0']
data.loc[data['2178-2.0']<0,'GeneralHealth-2.0']=np.nan
data['GeneralHealth-3.0']=5-data['2178-3.0']
data.loc[data['2178-3.0']<0,'GeneralHealth-3.0']=np.nan
data['GeneralHealth_pre_cl']=data['GeneralHealth-2.0']
data['GeneralHealth_post_cl']=data['GeneralHealth-3.0']
data['GeneralHealth_diff_cl']=data['GeneralHealth-3.0']-data['GeneralHealth-2.0']

data['Activity-2.0']=data['894-2.0']
data.loc[data['894-2.0']<0,'Activity-2.0']=np.nan
data.loc[data['894-2.0']>400,'Activity-2.0']=np.nan

data['Activity_days-2.0']=data['884-2.0']
data.loc[data['884-2.0']<0,'Activity_days-2.0']=np.nan
#data.loc[data['884-2.0']>400,'Activity_days-2.0']=np.nan

data['Activity_vig-2.0']=data['914-2.0']
data.loc[data['914-2.0']<0,'Activity_vig-2.0']=np.nan
#data.loc[data['914-2.0']>400,'Activity-2.0']=np.nan

data['Activity_vig_days-2.0']=data['904-2.0']
data.loc[data['904-2.0']<0,'Activity_vig-2.0']=np.nan

data['Activity_vig-3.0']=data['914-3.0']
data.loc[data['914-3.0']<0,'Activity_vig-3.0']=np.nan
#data.loc[data['914-2.0']>400,'Activity-2.0']=np.nan


data['TV-2.0']=data['1070-2.0']
data.loc[data['1070-2.0']<0,'TV-2.0']=np.nan
#data.loc[data['914-2.0']>400,'Activity-2.0']=np.nan

data.loc[:,'Education_age']=data.loc[:,find_vars(data,'^845-')].max(axis=1)

data['Activity-3.0']=data['894-3.0']
data.loc[data['894-3.0']<0,'Activity-3.0']=np.nan
data.loc[data['894-3.0']>400,'Activity-3.0']=np.nan


data['Activity-2.0']=data['894-2.0']
data.loc[data['894-2.0']<0,'Activity_levels-2.0']=np.nan
data.loc[data['894-2.0']>0,'Activity_levels-2.0']='Low'
data.loc[data['894-2.0']>40,'Activity_levels-2.0']='Medium'
data.loc[data['894-2.0']>100,'Activity_levels-2.0']='High'

data['Activity-3.0']=data['894-3.0']
data.loc[data['894-3.0']<0,'Activity-3.0']=np.nan
data.loc[data['894-3.0']>400,'Activity-3.0']=np.nan


# Smoking 

data.loc[data['1239-3.0']>0,'Smoking-3.0']=True
data.loc[data['1239-3.0']==0,'Smoking-3.0']=False
data.loc[data['1239-3.0']<0,'Smoking-3.0']=np.nan
data['Smoking-3.0']=data['Smoking-3.0'].convert_dtypes()
data['Smoking_bin-3.0']=data['Smoking-3.0'].astype(float)

    
data.loc[data['1239-2.0']>0,'Smoking-2.0']=True
data.loc[data['1239-2.0']==0,'Smoking-2.0']=False
data.loc[data['1239-2.0']<0,'Smoking-2.0']=np.nan
data['Smoking-2.0']=data['Smoking-2.0'].convert_dtypes()
data['Smoking_bin-2.0']=data['Smoking-2.0'].astype(float)

data['Wheeze-3.0']=data['2316-3.0']==1
data.loc[data['2316-3.0']<0,'Wheeze-3.0']=np.nan



# Plot of Vulnerability score
fig, ax = plt.subplots(figsize=(3,3))
rr=np.arange(50,90,0.5)
sns.lineplot(x=rr,y=10**(rr*0.0524-3.27))
ax=plt.gca()
ax.set_ylabel('Score')
ax.set_xlabel('Age')
ax.set_title('COVID-19 Vulnerability Score')
fig.savefig('Vulnerability.svg',dpi=300, bbox_inches = "tight")


# apply cleaning to proteomics 

prot_els_cl=[ a + '_regPl_diff_cl' for a in prot_els]

for a in prot_els_regPl:
    out = remove_conf(data,a,[],flatten=True,suf_pre=True,remove_out=True,suf=False)
    data.loc[:,a+'_pre_cl'] = out[a+'_pre_cl']
    data.loc[:,a+'_post_cl'] = out[a+'_post_cl']

    data.loc[:,a+'_diff_cl'] = data[a+'_post_cl']-data[a+'_pre_cl']
    data.loc[:,a+'_pc_cl'] = 100*data[a+'_diff_cl']/data[a+'_pre_cl']



# ## Reported Disease codes

# Self reported
disease_codes=(data[(find_vars(data,'20002-2'))].values.flatten())
disease_codes=disease_codes[~pd.isnull(disease_codes)]
codes_counts=(pd.Series((disease_codes))).value_counts()
codes_counts=codes_counts[codes_counts>20].drop(99999.0)
codes_counts.index=codes_counts.index.astype(int)
for a in codes_counts.index:
    data.loc[:,'dc_'+str(int(a))]=(data[find_vars(data,'20002-2')]==a).any(axis=1)

disease_codes_names=pd.read_csv('/home/eduff/biobank/phenotype/disease_codes.tsv',sep='\t',header=None)
disease_codes_names.set_index(0,inplace=True)
disease_codes_names.rename(columns={1:'name'},inplace=True)

data=data.copy()    


# Self reported

data['Diabetes2']=(data[(find_vars(data,'20002-'))]==1223).any(axis=1)
data['cardiovascular']=(data[(find_vars(data,'20002-'))]==1071).any(axis=1)
data['hypertension']=(data[(find_vars(data,'20002-'))]==1065).any(axis=1)
data['heart_cond']=(data[(find_vars(data,'20002-'))]==1066).any(axis=1)
data['depression']=(data[(find_vars(data,'20002-'))]==1312).any(axis=1)
data['cholesterol']=(data[(find_vars(data,'20002-'))]==1536).any(axis=1)
data['headinjury']=(data[(find_vars(data,'20002-'))]==1292).any(axis=1)
data['renal']=(data[(find_vars(data,'20002-'))]==1074).any(axis=1)
data['IBS']=(data[(find_vars(data,'20002-'))]==1154).any(axis=1)
data['COPD_ICD']=(data[(find_vars(data,'20002-'))]==1112).any(axis=1)
data['Bronc']=(data[(find_vars(data,'20002-'))]==1113).any(axis=1)

data['hearing-2.0']=data['2247-2.0']>0
data['hearing-3.0']=data['2247-3.0']>0


data['loneliness-2.0']=data['2020-2.0']>0
data['loneliness-3.0']=data['2020-3.0']>0

data['confide-2.0']=(data['2110-2.0']>-1)&(data['2110-2.0']<2)
data['confide-3.0']=(data['2110-3.0']>-1)&(data['2110-3.0']<2)

data['social_visits-2.0']=(data['1031-2.0']>3)
data['social_visits-3.0']=(data['1031-3.0']>3)

data['isolation-2.0']=data['confide-2.0']&data['social_visits-2.0']
data['isolation-2.0']=data['confide-3.0']&data['social_visits-3.0']

# define hypertension

data['hypertension_emp-2.0']=(data['4080-2.0']>130) & (data['4079-2.0']>80)
data['hypertension_emp-3.0']=(data['4080-3.0']>130) & (data['4079-3.0']>80)


# ## COVID reports

# process records of COVID

data.loc[:,'GP']=((data.loc[all_case,['41001-3.0','41001-3.1','41001-3.2','41001-3.3']])==2).any(axis=1)
data.loc[:,'HES']=((data.loc[all_case,['41001-3.0','41001-3.1','41001-3.2','41001-3.3']])==1).any(axis=1)
data.loc[:,'PCR']=((data.loc[all_case,['41001-3.0','41001-3.1','41001-3.2','41001-3.3']])==3).any(axis=1)
data.loc[:,'Lateral']=((data.loc[all_case,['41001-3.0','41001-3.1','41001-3.2','41001-3.3']])==4).any(axis=1)
data.loc[:,'PublicHealthRecord']=(data.loc[:,'HES'])|(data.loc[:,'GP'])

# characterise individuals only identified by antibody test (no reported symptoms)

data.loc[:,'COVID'] = np.nan
els = all_matched & all_case & (data.loc[:,['HES','PCR','GP']].sum(axis=1)==0)  & (data.loc[:,'Lateral']==True)
data.loc[els,'COVID'] = 'SARS'
data.loc[data.loc[els,'matched_eid'].dropna(),'COVID'] = 'SARS_ctr'
els = all_case & (data.loc[:,['HES','PCR','GP']].sum(axis=1)>0)  
data.loc[els,'COVID'] = 'COVID'
data.loc[data.loc[els,'matched_eid'].dropna(),'COVID'] = 'COVID_ctr'

# load up reginal  COVID reports
England = pd.read_csv('/home/eduff/biobank/phenotype/covid19_result_england.txt', sep='\t')
Scotland = pd.read_csv('/home/eduff/biobank/phenotype/covid19_result_scotland.txt', sep='\t')
Wales = pd.read_csv('/home/eduff/biobank/phenotype/covid19_result_wales.txt', sep='\t')

#print(England[England.eid.isin(data.index)].shape)
#print(Scotland[Scotland.eid.isin(data.index)].shape)
#print(Wales[Wales.eid.isin(data.index)].shape)

# remove subjects not in dataset
del Wales
England=England[England.eid.isin(data.index)]
Scotland=Scotland[Scotland.eid.isin(data.index)]               

Scotland.loc[:,'spectype']=np.nan
Scotland.loc[:,'origin']=np.nan
Scotland.loc[:,'acute']=np.nan
Scotland.loc[:,'hosaq']=np.nan
England.loc[:,'site']=np.nan

COVIDreports=pd.concat([England,Scotland])
COVIDreports.loc[:,'specdate']=pd.to_datetime(COVIDreports['specdate'],infer_datetime_format=True)

# convert to lists for each id

COVIDreports=COVIDreports.groupby(['eid']).agg(list)

# calculate time since first and most recent +ve test when available

for a in COVIDreports.index:
    
    # calculate first and last positive test before 3rd assessment
    results_pos=np.array(COVIDreports.loc[a,'result'])==1
    dates_pos=np.array(COVIDreports.loc[a,'specdate'])[results_pos]
    dates_pos=dates_pos[dates_pos < pd.to_datetime(data.loc[a,'53-3.0'])]
    if len(dates_pos)>0:
        COVIDreports.loc[a,'firstpos']=dates_pos.min()
        COVIDreports.loc[a,'lastpos']=dates_pos.max()
        COVIDreports.loc[a,'time_since_pos_d']=(pd.to_datetime(data.loc[a,'53-3.0']) - dates_pos.max())
        
        COVIDreports.loc[a,'time_since_first_pos_d']=(pd.to_datetime(data.loc[a,'53-3.0'])  - dates_pos.min())
        
    results_neg=np.array(COVIDreports.loc[a,'result'])==0
    dates_neg=np.array(COVIDreports.loc[a,'specdate'])[results_neg]  
    dates_neg=dates_neg[dates_neg < pd.to_datetime(data.loc[a,'53-3.0'])]
    
    if len(dates_neg) > len(dates_pos):
        COVIDreports.loc[a,'neg_reports']=True

COVIDreports.loc[:,'time_since_pos']=COVIDreports.loc[:,'time_since_pos_d'].dt.days
COVIDreports.loc[:,'time_since_first_pos']=COVIDreports.loc[:,'time_since_first_pos_d'].dt.days


data=pd.merge(data, COVIDreports, left_index=True, right_index=True,how='left',suffixes=("_x",None))


# ## Define case / control variables and helper vars

# hospitalisation and non-COVID

hosp=np.sum(np.max(data.loc[all_case,find_vars(data,'41270')]=='U071',axis=1))
main_hosp=(data.loc[:,find_vars(data,'41202')]=='U071').any(axis=1)
secondary_hosp=(data.loc[:,find_vars(data,'41204')]=='U071').any(axis=1)&~main_hosp
data.loc[:,'main_hosp']=np.max(data.loc[:,find_vars(data,'41202')]=='U071',axis=1)

# hospitalisation dates
vars_h=(find_vars(data,'41270'))

[rws,cls] = np.where(data.loc[:,vars_h]=='U071')
for el in range(len(cls)):
    eid=data.index[rws[el]]
    col= '41280-0.'+ vars_h[cls[el]][8:]
    date_pos = pd.to_datetime(data.loc[eid,col])
    if date_pos < pd.to_datetime(data.loc[eid,'53-3.0']):
        if pd.isnull(data.loc[eid,'firstpos']):
            
            data.loc[eid,'firstpos']=date_pos
            data.loc[eid,'lastpos']=date_pos
        else:
            if date_pos < data.loc[eid,'firstpos']:
                
                data.loc[eid,'firstpos']=date_pos
            if date_pos > data.loc[eid,'lastpos']:
                
                data.loc[eid,'lastpos']=date_pos
      
# hospitalisation    
data.loc[data['Case']=='ctr','Case_hosp'] = 'ctr'
data.loc[data['Case']=='sars','Case_hosp'] = 'sars'
data.loc[data['main_hosp'],'Case_hosp'] = 'sars_hosp'

data.loc[:,'Case_hosp_matched']=np.nan
data.loc[data['Case']=='sars','Case_hosp_matched'] = 'sars'
data.loc[data['main_hosp'],'Case_hosp_matched'] = 'sars_hosp'
data.loc[data.loc[data['Case_hosp']=='sars','matched_eid'].dropna(),'Case_hosp_matched'] = 'ctr_sars' 
data.loc[data.loc[data['Case_hosp']=='sars_hosp','matched_eid'].dropna(),'Case_hosp_matched'] = 'ctr_sars_hosp' 

data.loc[:,'Case_hosp_bin']=0
data.loc[data.loc[:,'Case_hosp_matched']=='sars_hosp','Case_hosp_bin']=0.5
data.loc[data.loc[:,'Case_hosp_matched']=='ctr_sars_hosp','Case_hosp_bin']=-0.5

data.loc[:,'Case_hosp_bin_only']=0
data.loc[data.loc[:,'Case_hosp_matched']=='sars_hosp','Case_hosp_bin_only']=0.5
data.loc[data.loc[:,'Case_hosp_matched']=='ctr_sars_hosp','Case_hosp_bin_only']=-0.5

data.loc[:,'Case_nohosp_bin_only']=0
data.loc[data.loc[:,'Case_hosp_matched']=='sars','Case_nohosp_bin_only']=0.5
data.loc[data.loc[:,'Case_hosp_matched']=='ctr_sars','Case_nohosp_bin_only']=-0.5
    
data['Case_hosp_bin_neg']=data['Case_hosp_bin']
els=data['Case_hosp_bin']>0
data.loc[data.loc[els,'matched_eid'].dropna(),'Case_hosp_bin_neg']=-1


data.loc[:,'COVID_bin']=0
data.loc[data['COVID']=='COVID_ctr','COVID_bin']=-1
data.loc[data['COVID']=='COVID','COVID_bin']=1
    
# COVID death (0)
np.sum(np.max(data.loc[all_case,find_vars(data,'40002')]=='U071',axis=1))
np.sum(np.max(data.loc[all_case,find_vars(data,'40001')]=='U071',axis=1))
np.sum(np.max(data.loc[all_case,find_vars(data,'40007')],axis=1))

# Wave 1-6 Hospitalisation report
np.sum(data.loc[:,find_vars(data,'28029')]=='1')

# first vacination prior to pos
data.loc[:,'vac_prior_first_pos']=np.nan  
diff=(pd.to_datetime(data.loc[:,'27984-0.0']) - COVIDreports.loc[:,'firstpos']).dt.days
data.loc[diff<0,'vac_prior_first_pos']='prior'
data.loc[diff>0,'vac_prior_first_pos']='post'

# first vaccination prior to assessment 
data.loc[:,'vac_prior_assess']=np.nan  
diff=(pd.to_datetime(data.loc[:,'27984-0.0']) - pd.to_datetime(data.loc[:,'53-3.0'])).dt.days
data.loc[diff<0,'vac_prior_assess']=True
data.loc[diff>0,'vac_prior_assess']=False

# first vaccination prior to pos
data.loc[:,'vac_prior_first_pos']=np.nan  
diff=(pd.to_datetime(data.loc[:,'27984-0.0']) - COVIDreports.loc[:,'firstpos']).dt.days
data.loc[diff<0,'vac_prior_first_pos']=True
data.loc[diff>0,'vac_prior_first_pos']=False


# vaccination
data.loc[:,'Case_vax_prior_bin']=0
data.loc[data['vac_prior_first_pos']=='prior','Case_vax_prior_bin']=1
els=data['Case_vax_prior_bin']>0
#data.loc[data.loc[els,'matched_eid'].dropna(),'Case_vax_prior_bin']=-1
data.loc[data['vac_prior_first_pos']=='post','Case_vax_prior_bin']=-1

data['Case_vax']=data['Case']
data.loc[data['vac_prior_first_pos']=='prior','Case_vax']='vax'


# first vaccination prior to assessment 
data.loc[:,'vac_prior_assess']=np.nan  
diff=(pd.to_datetime(data.loc[:,'27984-0.0']) - pd.to_datetime(data.loc[:,'53-3.0'])).dt.days
data.loc[diff<0,'vac_prior_assess']=True
data.loc[diff>0,'vac_prior_assess']=False

# first vaccination prior to pos
data.loc[:,'vac_prior_first_pos']=np.nan  
diff=(pd.to_datetime(data.loc[:,'27984-0.0']) - COVIDreports.loc[:,'firstpos']).dt.days
data.loc[diff<0,'vac_prior_first_pos']=True
data.loc[diff>0,'vac_prior_first_pos']=False


# vaccinations
#print(data.loc[all_matched,'vac_prior_first_pos'].sum())
#print(data.loc[all_matched,'vac_prior_first_pos'].notna().sum())


# CHECK THRIVA
## THRIVA - followup antibody test
THRIVA_ids=data.loc[(data.loc[:,'27990-0.0']==1) & (data.loc[:,'Case']=='sars'),'matched_eid'].dropna().index
THRIVA_ctr_ids=data.loc[THRIVA_ids,'matched_eid']
THRIVA_ids=data.loc[(data.loc[:,'27990-0.0']==1) & (data.loc[:,'Case']=='sars'),'matched_eid'].dropna().index

data.loc[:,'THRIVA_pos']=np.NaN
data.loc[THRIVA_ids,'THRIVA_pos']='THRIVA_pos'
data.loc[THRIVA_ctr_ids,'THRIVA_pos']='ctr_THRIVA_pos'


# self tests:  27981, 27990, 28140
data.loc[((data.loc[:,'27990-0.0']!=1) | (data.loc[:,'27981-0.0']!=1) | (data.loc[:,'28140-0.0']!=1) ) & (data.loc[:,'Case']=='ctr')]

#print(len(data.loc[(data.loc[:,'27981-0.0']==0) & (data.loc[:,'Case']=='ctr')]))

#print(len(data.loc[(data.loc[:,'28140-0.0']==0) & (data.loc[:,'Case']=='ctr')]))


# calculate post - pre differences and percentage change for a variety of variables

prot_els_diff=[]    
for IDP in prot_els + prot_els_regPl:
    name=IDP+'_diff'
    prot_els_diff.append(name)
    data.loc[:,name]=(data.loc[:,IDP+'_post']) - (data.loc[:,IDP+'_pre'])
    
# also adjust original verions of proteins log normalised   
for IDP in ['pTau-181','GFAP','NfL']:
    IDP=IDP+'_orig'
    name=IDP+'_diff'
    prot_els_diff.append(name)
    data.loc[:,name]=(data.loc[:,IDP+'_post']) - (data.loc[:,IDP+'_pre'])

for IDP in ['27183-2.0','27004-2.0','26593-2.0','6350-2.0','6348-2.0','27163-2.0','27202-2.0','27200-2.0','27293-2.0','6350-2.0','Age-2.0','4080-2.0','4079-2.0','21002-2.0','21001-2.0','Alcohol-2.0','738-2.0','709-2.0','6350-2.0','20018-2.0','2178-2.0','25781-2.0']:
    IDP=IDP[:-4]
    data.loc[:,IDP+'_diff']=(data.loc[:,IDP+'-3.0']) - (data.loc[:,IDP+'-2.0'])
    
IDP='27183' # lOBC (thick)
IDP='27004'  # lpHg (cont)
IDP='26593' # hipp volume

for IDP in ['27183-2.0','27004-2.0','26593-2.0','6350-2.0','6348-2.0','27163-2.0','27202-2.0','27200-2.0','27293-2.0','6350-2.0','Age-2.0','4080-2.0','4079-2.0','21002-2.0','21001-2.0','Alcohol-2.0','738-2.0','709-2.0','6350-2.0','20018-2.0','2178-2.0','25781-2.0']:
    IDP=IDP[:-4]
    data.loc[:,IDP+'-pc']=100*(data.loc[:,IDP+'-3.0']-data.loc[:,IDP+'-2.0'])/data.loc[:,IDP+'-2.0']


# percentage change
prot_els_pc=[]    
for IDP in prot_els + prot_els_regPl:
    name=IDP+'-pc'
    prot_els_pc.append(name)
    data.loc[:,name]=100*(data.loc[:,IDP+'_post']-data.loc[:,IDP+'_pre'])/data.loc[:,IDP+'_pre']
    
for IDP in ['pTau-181','GFAP','NfL']:
    IDP=IDP+'_orig'
    name=IDP+'-pc'
    prot_els_pc.append(name)
    data.loc[:,name]=100*(data.loc[:,IDP+'_post']-data.loc[:,IDP+'_pre'])/data.loc[:,IDP+'_pre']
    
for IDP in ['pTau-181','GFAP','NfL']:
    IDP=IDP+'_orig_regPl'
    name=IDP+'-pc'
    prot_els_pc.append(name)
    data.loc[:,name]=100*(data.loc[:,IDP+'_post']-data.loc[:,IDP+'_pre'])/data.loc[:,IDP+'_pre']
    
data=data.copy()


# employment information
for a in ['2','3']:
    data['paid_employ-'+a+'.0']=(data.loc[:,find_vars(data,'^6142-'+a+'.')]==1).sum(axis=1)>0
    data.loc[:,'retired-'+a+'.0']=(data.loc[:,find_vars(data,'^6142-'+a+'.')]==2).sum(axis=1)>0
    data.loc[:,'caring-'+a+'.0']=(data.loc[:,find_vars(data,'^6142-'+a+'.')]==3).sum(axis=1)>0
    data.loc[:,'no_employ-sick-'+a+'.0']=(data.loc[:,find_vars(data,'^6142-'+a+'.')]==4).sum(axis=1)>0
    data.loc[:,'unemployed-'+a+'.0']=(data.loc[:,find_vars(data,'^6142-'+a+'.')]==4).sum(axis=1)>0
    data.loc[:,'voluntary-'+a+'.0']=(data.loc[:,find_vars(data,'^6142-'+a+'.')]==6).sum(axis=1)>0
    negs=((data.loc[:,find_vars(data,'^6142-'+a+'.')]).sum(axis=1)<0)
    data.loc[negs,'paid_employ-'+a+'.0']=np.nan
    data.loc[negs,'retired-'+a+'.0']=np.nan
    data.loc[negs,'caring-'+a+'.0']=np.nan
    data.loc[negs,'no_employ-sick-'+a+'.0']=np.nan
    data.loc[negs,'unemployed-'+a+'.0']=np.nan
    data.loc[negs,'voluntary-'+a+'.0']=np.nan
 
    # not converting paid_emply to not flip in table alg
    data['paid_employ-'+a+'.0']=data['paid_employ-'+a+'.0'].convert_dtypes()
    data['paid_employ_bin-'+a+'.0']=data['paid_employ-'+a+'.0'].astype(float)
    data['retired-'+a+'.0']=data['retired-'+a+'.0'].convert_dtypes()
    data['retired_bin-'+a+'.0']=data['retired-'+a+'.0'].astype(float)
    data['caring-'+a+'.0']=data['caring-'+a+'.0'].convert_dtypes()
    data['no_employ-sick-'+a+'.0']=data['no_employ-sick-'+a+'.0'].convert_dtypes()
    data['unemployed-'+a+'.0']=data['unemployed-'+a+'.0'].convert_dtypes()
    data['voluntary-'+a+'.0']=data['voluntary-'+a+'.0'].convert_dtypes()


# # Preprocess Imaging and Cognitive vars
# 

# Imaging variables for AD brain phenotypes
#26562	Volume of Hippocampus (left hemisphere)	Freesurfer ASEG  
#26593	Volume of Hippocampus (right hemisphere)	Freesurfer ASEG  
#26641	Volume of Whole-hippocampus (left hemisphere)	Freesurfer subsegmentation  
#26663	Volume of Whole-hippocampus (right hemisphere)	Freesurfer subsegmentation  
#25886	Volume of grey matter in Hippocampus (left)	Regional grey matter volumes (FAST)  
#25887	Volume of grey matter in Hippocampus (right)	Regional grey matter volumes (FAST)  
#25019	Volume of hippocampus (left)	Subcortical volumes (FIRST)  
#5020	Volume of hippocampus (right)
    
#27432	Mean thickness of G-precuneus (left hemisphere)	Freesurfer a2009s  
#27654	Mean thickness of G-precuneus (right hemisphere)	Freesurfer a2009s  
#26779	Mean thickness of precuneus (left hemisphere)	Freesurfer desikan white  
#27196	Mean thickness of precuneus (left hemisphere)	Freesurfer DKT  
#26880	Mean thickness of precuneus (right hemisphere)	Freesurfer desikan white  
#27289	Mean thickness of precuneus (right hemisphere)	Freesurfer DKT  
# WM hyperintensities 
np.sum(data.loc[all_matched & (data['Case']=='sars'),'25781-2.0']<1000)

class Cols:
    Weight = "21002"
    Hip = "48"
    BMI = "49"
    Alcohol = "1558"
    BP_sys = "4080"
    BP_dia = "4079"
    lOBC='27183' # lOBC (thick)
    lpHg='27004'  # lpHg (cont)
    HippVol='26593' # hipp volume
    EntL_DKT_mt='27177'
    EntR_DKT_mt='27270'


#  RSNs 
RSNs=np.array(['DMN1','Vis1','SMot1','DMN2','Att1','Att2','DMN3','Vis2','Vis3','SMot2','SMot3','SMot4','Temp1','Front1','Cb1','DMN4','Front2','BG1','Vis4','DMN5','Temp2'])
DMNs=['DMN' in a for a in RSNs]
Viss=['Vis' in a for a in RSNs]
SMots=['SMot' in a for a in RSNs]
Atts=['Att' in a for a in RSNs]
Cbs=['Cb' in a for a in RSNs]
BGs=['BG' in a for a in RSNs]
Fronts=['Front' in a for a in RSNs]
Temps=['Temp'in a for a in RSNs]
RSN_groups=[DMNs,Viss,SMots,Atts,Cbs,BGs,Fronts,Temps]
RSN_group_names=['DMN','Vis','SMot','Att','Cb','BG','Front','Temp']

# core confound variables
# imaging vars 

# trail making
for a in ['6350','6348']:
    out=(remove_conf(data,a,pos=True,remove_out=True,log=True,flatten=True))
    data.loc[:,[a+'_cl-2.0',a+'_cl-3.0']]=out.loc[:,[a+'-2.0_cl',a+'-3.0_cl']].values
# RT
data.loc[:,'RT-2.0']=np.log(data.loc[:,'20023-2.0'])
data.loc[:,'RT-3.0']=np.log(data.loc[:,'20023-3.0'])

# pairs matching
data.loc[:,'Pairs-2.0']=np.log(1+data.loc[:,'399-2.2'])
data.loc[:,'Pairs-3.0']=np.log(1+data.loc[:,'399-3.2'])
      
# Prospective memory 20018
data.loc[:,'ProsMem-2.0']=np.mod(data.loc[:,'20018-2.0'],2)
data.loc[:,'ProsMem-3.0']=np.mod(data.loc[:,'20018-3.0'],2)

# Fluid Intell 20016
                 
# Numeric memory 4282
data.loc[:,'NumMem-2.0'] = data.loc[:,'4282-2.0'] 
data.loc[data.loc[:,'NumMem-2.0']==-1,'4282-2.0'] = np.nan
data.loc[:,'NumMem-3.0'] = data.loc[:,'4282-3.0'] 
data.loc[data.loc[:,'NumMem-3.0']==-1,'4282-3.0'] = np.nan

# digit substitution 23324

# trail making b, pairs(2), Reaction Time Test ,

# picture vocab 6364

# paired associate 20197

# tower rearranging 21004

# matrix pattern completion 6373

cog_vars=['6350','6348','4282','20018']
cog_vars_gi=['Pairs','RT','ProsMem','20016','NumMem','6350_cl','23324','6364','20197','21004','6373']
cog_vars_gi_weights=np.array([-0.49,-0.41,0.54,0.66,0.58,-0.79,0.73,0.19,0.49,0.65,0.66])
cog_vars_alz_weights=np.array([0,0,0,0.66,0.58,0,0,0,0,0,0.66])

for a in ['-2.0','-3.0']:
    cols=[col+a for col in cog_vars_gi]
    df=data.loc[:,[]]
    for id in data.index:
        if np.sum(data.loc[id,cols].notna())==11:
            df.loc[id,'cog_vars_gi'+a]=np.dot(data.loc[id,cols],cog_vars_gi_weights)
            df.loc[id,'cog_vars_alz'+a]=np.dot(data.loc[id,cols],cog_vars_alz_weights)
        else:
            df.loc[id,'cog_vars_gi'+a]=np.nan
            df.loc[id,'cog_vars_alz'+a]=np.nan         

    data.loc[:,'cog_vars_gi'+a]=df.loc[:,'cog_vars_gi'+a]
    data.loc[:,'cog_vars_alz'+a]=df.loc[:,'cog_vars_alz'+a]

data.loc[:,'cog_vars_gi_diff'] = data.loc[:,'cog_vars_gi-3.0'] - data.loc[:,'cog_vars_gi-2.0']
data.loc[:,'cog_vars_gi_diff_cl'] = data.loc[:,'cog_vars_gi_diff']
data.loc[:,'cog_vars_alz_diff'] = data.loc[:,'cog_vars_alz-3.0'] - data.loc[:,'cog_vars_alz-2.0']
data.loc[:,'cog_vars_alz_diff_cl'] = data.loc[:,'cog_vars_alz_diff'] 

data.loc[:,'cog_vars_gi_post_cl'] = data.loc[:,'cog_vars_gi-3.0'] 
data.loc[:,'cog_vars_gi_pre_cl']= data.loc[:,'cog_vars_gi-2.0']
data.loc[:,'cog_vars_alz_post_cl'] = data.loc[:,'cog_vars_alz-3.0'] 
data.loc[:,'cog_vars_alz_pre_cl']= data.loc[:,'cog_vars_alz-2.0']

# UKB Pairs Matching Test -0.49 
# UKB Reaction Time Test -0.41
# UKB Prospective Memory Test 0.54 
# UKB Fluid Intelligence Test 0.66
# UKB Numeric Memory Test 0.58 
# UKB Trail Making Test part B -0.79 
# UKB Symbol Digit Test 0.73 
# UKB Picture Vocabulary Test 0.19 
# UKB Paired Associate Learning Test 0.49
# UKB Tower Test 0.65 
# UKB Matrix Pattern Test 0.66 
    
IDP_info=pd.read_excel('/home/eduff/biobank/imaging/media-2.xlsx')
    
imaging_vars= ['26593','27177','27004','27183','25093','25494','26562','25793','27202','27200','27293','27146','27239']
hipp_vars=['26562','26593', '26641', '26663', '25886', '25887', '25019', '25020']
prec_vars=['27432' , '27654'  , '26779', '27196', '26880' , '27289']
ent_vars=['27177','27270']
RSN_ampvars=[a + '_Amp' for a in RSNs ]
RSN_corrvars=['DMN_'+a + '_Corr' for a in RSN_group_names[1:] ]
imaging_vars=imaging_vars+hipp_vars+prec_vars+ent_vars # +RSN_ampvars+RSN_corrvars
#imaging_vars_ro_norm = [ a + '_norm' for a in imaging_vars]

#[head_size_scaling conf_TablePos_COG_Table conf_TablePos_COG_Z conf_YTRANS age];
CONF1_vars = ['25000', '25759',  '25758', '25757', 'Age']
#CONF1_vars_ro = [a +'_ro' for a in ['25000', '25759',  '25758', '25757', 'Age']]

# second set of confound IDPs
# Weight, Waist Circ, Hip Circ, BMI, BPdia, BPsys, SMOKING, Alcohol, Diabetes, Chol_med, BP_meds,'Deprivation'
CONF2_vars = ['21002','48','49','21001','4079', '4080', 'Smoking','Alcohol','Diabetes','Chol_meds','BP_meds','Deprivation']

# remove outliers conf vars
#data.loc[:,CONF1_vars_ro]=remove_outliers(data.loc[all_matched,CONF1_vars],15,repl=0,visits=False).values
#data.loc[:,imaging_vars_ro]=remove_outliers(data.loc[:,imaging_vars],8).values


for a in imaging_vars: 
    # deconfound headsize etc
    data.loc[:,[a+'-2.0_cl',a+'-3.0_cl']]=remove_conf(data,a,CONF1_vars[:-1],remove_out=True)

    # remove difference outliers 
    CONF=['Age-3.0_f','assessment_sep','assessment_sep^2','Ethnicity(White)','Gender_pre']
    naels=remove_outliers(remove_conf(data,a+'-3.0_cl',CONF+[a+'-2.0_cl'],flatten=False,suf=False)).isna()

    for pref in ['-3.0_cl','-2.0_cl']:
        data.loc[naels.squeeze(),a+pref]=np.nan
    
    data.loc[:,a+'_diff_cl'] = data[a+'-3.0_cl']-data[a+'-2.0_cl']
    data.loc[:,a+'_pc_cl'] = 100*data[a+'_diff_cl']/data[a+'-2.0_cl']
    

    


#os.chdir('/home/eduff/biobank/'
# names of IDPs
IDP_info=pd.read_excel('/home/eduff/biobank/imaging/media-2.xlsx',)
IDP_info.set_index('UKB ID')
IDP_info['IDP short name']=IDP_info['IDP short name'].str.replace(" ","_")
IDP_info['IDP short name'][IDP_info['IDP short name'].str.contains('bankssts')]


# ### Extract ADNI AD signature 

features_cort = ['thickness_bankssts_lh',  'thickness_caudalanteriorcingulate_lh', 'thickness_caudalmiddlefrontal_lh',  'thickness_cuneus_lh',  'thickness_entorhinal_lh',  'thickness_fusiform_lh',  'thickness_inferiorparietal_lh',  'thickness_inferiortemporal_lh',  'thickness_isthmuscingulate_lh',  'thickness_lateraloccipital_lh',  'thickness_lateralorbitofrontal_lh',  'thickness_lingual_lh',  'thickness_medialorbitofrontal_lh',  'thickness_middletemporal_lh',  'thickness_parahippocampal_lh',  'thickness_paracentral_lh',  'thickness_parsopercularis_lh',  'thickness_parsorbitalis_lh',  'thickness_parstriangularis_lh',  'thickness_pericalcarine_lh',  'thickness_postcentral_lh',  'thickness_posteriorcingulate_lh',  'thickness_precentral_lh',  'thickness_precuneus_lh',  'thickness_rostralanteriorcingulate_lh', 'thickness_rostralmiddlefrontal_lh',  'thickness_superiorfrontal_lh',  'thickness_superiorparietal_lh',  'thickness_superiortemporal_lh',  'thickness_supramarginal_lh',  'thickness_frontalpole_lh',  'thickness_temporalpole_lh',  'thickness_transversetemporal_lh',  'thickness_insula_lh',  'thickness_bankssts_rh',  'thickness_caudalanteriorcingulate_rh', 'thickness_caudalmiddlefrontal_rh',  'thickness_cuneus_rh',  'thickness_entorhinal_rh',  'thickness_fusiform_rh',  'thickness_inferiorparietal_rh',  'thickness_inferiortemporal_rh',  'thickness_isthmuscingulate_rh',  'thickness_lateraloccipital_rh',  'thickness_lateralorbitofrontal_rh',  'thickness_lingual_rh',  'thickness_medialorbitofrontal_rh',  'thickness_middletemporal_rh',  'thickness_parahippocampal_rh',  'thickness_paracentral_rh',  'thickness_parsopercularis_rh',  'thickness_parsorbitalis_rh',  'thickness_parstriangularis_rh',  'thickness_pericalcarine_rh',  'thickness_postcentral_rh',  'thickness_posteriorcingulate_rh',  'thickness_precentral_rh',  'thickness_precuneus_rh',  'thickness_rostralanteriorcingulate_rh', 'thickness_rostralmiddlefrontal_rh',  'thickness_superiorfrontal_rh',  'thickness_superiorparietal_rh',  'thickness_superiortemporal_rh',  'thickness_supramarginal_rh',  'thickness_frontalpole_rh',  'thickness_temporalpole_rh',  'thickness_transversetemporal_rh',  'thickness_insula_rh']
missing_feats=[]
features_cort_ids=[]
for a in features_cort:
    names=IDP_info['IDP short name'][IDP_info['IDP short name'].str.contains(a[:-3])]
    ukb_name=([s for s in names.values if re.search('Desikan*'+a[-3:],s)])
    if ukb_name == []:
        missing_feats.append(a)
        #print(a)
    else:
        features_cort_ids.append(int(IDP_info.loc[IDP_info['IDP short name']==ukb_name[0],'UKB ID']))

features_vol = [ a.replace('thickness','volume') for a in features_cort]
features_vol_ids=[]
for a in features_vol:
    names=IDP_info['IDP short name'][IDP_info['IDP short name'].str.contains(a[:-3])]
    ukb_name=([s for s in names.values if re.search('Desikan*'+a[-3:],s)])
    if ukb_name == []:
        missing_feats.append(a)
        #print(names)
        #print(a)
    else:
        features_vol_ids.append(int(IDP_info.loc[IDP_info['IDP short name']==ukb_name[0],'UKB ID']))



features_vol_extra = ['volume_Left-Cerebellum-White-Matter', 'volume_Left-Cerebellum-Cortex',
                          'volume_Left-Thalamus-Proper', 'volume_Left-Caudate', 'volume_Left-Putamen',
                          'volume_Left-Pallidum','volume_Brain-Stem', 'volume_Left-Hippocampus', 
                          'volume_Left-Amygdala', 'volume_Left-Accumbens-area', 'volume_Right-Cerebellum-White-Matter',
                          'volume_Right-Cerebellum-Cortex', 'volume_Right-Thalamus-Proper', 'volume_Right-Caudate',
                          'volume_Right-Putamen', 'volume_Right-Pallidum', 'volume_Right-Hippocampus', 
                          'volume_Right-Amygdala', 'volume_Right-Accumbens-area']
#print(len(features_vol_extra))
features_vol_extra_ids = [26556,26557,26558,26559,26560,26561,26526,26562,26563,26564, 26587,26588,26559,26590,26591,26592,26593,26594,26595]
len(features_vol_extra_ids) 


df_ml = pd.DataFrame(columns=features_cort+features_vol+features_vol_extra)
all_ni_ad_feats=features_cort+features_vol+features_vol_extra
pres_ni_ad_feats_ids=features_cort_ids+features_vol_ids+features_vol_extra_ids
pres_ni_ad_feats=([a for a in all_ni_ad_feats if a not in missing_feats ])


ad_feats={}
for session in ['2','3']:
    sess_str='-'+session+'.0'
    ad_feats[session]=data.loc[:,[str(a)+sess_str for a in pres_ni_ad_feats_ids]]
    ad_feats[session]=ad_feats[session].set_axis(pres_ni_ad_feats,axis=1)
    for a in missing_feats:
        ad_feats[session].loc[:,a]=0

    ad_feats[session].loc[:,'gender'] = data.loc[:,'Gender_pre']
    ad_feats[session].loc[:,'age'] = data.loc[:,'Age'+sess_str]
    ad_feats[session].loc[:,'eTIV'] = data.loc[:,'26521'+sess_str]
 


# Preprocess ADNI data

import pandas as pd
import statsmodels.api as smfapi

from sklearn.preprocessing import StandardScaler
import joblib 

def calculate_residuals(data, columns, name_regress):
    for c in columns:
        lin_model = smfapi.load(f'../../imaging/adni_phenotypes/regress_models/{name_regress}_{c}.pkl')
        
        data[c] = data[c] - lin_model.predict(data)
    return data


# Calculate residuals and scale for ukb data (load model and calculate from previously fitted adni data)
scaler = joblib.load('../../imaging/adni_phenotypes/data/adni_train_scaler.joblib')
scaler.mean_

ad_feats_scaled={}

for session in ['2','3']:
    ad_feats[session]=calculate_residuals(data=ad_feats[session], columns=features_cort, name_regress='adni')
    ad_feats[session]=calculate_residuals(data=ad_feats[session], columns=features_vol+features_vol_extra, name_regress='adni')

# Removing demographic columns and putting in the right order for scaler
    ad_feats[session]=ad_feats[session][features_cort + features_vol + features_vol_extra]
    ad_feats_scaled[session] = pd.DataFrame(scaler.transform(ad_feats[session]), columns=ad_feats[session].columns, index=ad_feats[session].index)

    ad_feats_scaled[session].to_csv('ukb_scaled_corrected_'+session+'.csv')


ad_results={}
for session in ['2','3']:
    ad_results[session]=pd.read_csv('../../imaging/adni_phenotypes/results/adni_results_'+session+'.csv')
    ad_results[session]=ad_results[session].set_index('ukb_id')
    for b in ['mean','std']:
        data.loc[:,'adni_'+b+'-'+session+'.0']=ad_results[session][b]


#data.loc[:,'adni-3.0']=
for a in ['adni_mean']:
    out=(remove_conf(data,a,pos=True,remove_out=True,log=True,flatten=True))
    data.loc[:,[a+'_cl-2.0',a+'_cl-3.0']]=out.loc[:,[a+'-2.0_cl',a+'-3.0_cl']].values
data.loc[:,'adni_mean_cl_diff']=data.loc[:,'adni_mean_cl-3.0']-data.loc[:,'adni_mean_cl-2.0']
data.loc[:,'adni_mean_diff']=data.loc[:,'adni_mean-3.0']-data.loc[:,'adni_mean-2.0']
data['adni_mean_pre_cl']=data['adni_mean_cl-2.0']
data['adni_mean_post_cl']=data['adni_mean_cl-3.0']
data['adni_mean_diff_cl']=data['adni_mean_cl_diff']


# # Data for Liz
# liz_vars=[]
# for a in ['21001','4079','4080','102','6152','6150','2443','2453','6142','709','738','1558','20116',
# '21022','31','20118','22189','21000','6138','54']:
#     liz_vars=liz_vars+find_vars(data,"^"+a+"-")
# with open('liz_vars.txt', 'w') as f:
#     for line in liz_vars:
#         f.write(f"{line}\n")
# f.close()
# data.loc[:,liz_vars].to_csv('liz_vars.csv')


