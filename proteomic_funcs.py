import pickle as pkl
import pandas
import pandas as pd
from pandas.api.types import is_numeric_dtype
from pandas.api.types import is_bool_dtype
import numpy as np
import pickle as pkl 
import os,re
import scipy
scipy.__version__
import matplotlib.pyplot as plt

import scipy.stats as stats
import sklearn
import statsmodels
import statsmodels.formula.api as smf
import seaborn as sns
import xarray

from scipy.stats import chi2_contingency
import numbers,numpy.random
from statsmodels.multivariate.pca import PCA
from sklearn.datasets import make_regression

# helper functions
# set of helper functions for biobank analysis


# stats class to store basic outputs from regresion fits
class basic_stats:
  def __init__(self, stats):
    self.pvalues = stats.pvalues
    self.params = stats.params
    self.rsquared_adj = stats.rsquared_adj
    self.conf_int = stats.conf_int()
    self.df_model = stats.df_model
    self.df_resid = stats.df_resid

    
class ols_simp:
  def __init__(self, formula='',data='',pre=[],case='Case_bin'):
    self.formula=formula
    self.data=data
    self.pre=pre
    self.case=case

  def fit(self):
    fit=(smf.ols(formula=self.formula, data=self.data,missing="drop").fit())
    out=basic_stats(fit)
    sd_x=fit.model.data.exog.std(axis=0)
    sd_y=fit.model.data.endog.std(axis=0)
    conf_int = fit.conf_int()
    coefficients=fit.params
    beta_coefficients = pd.Series()
    conf_int_norm=conf_int.copy()
    

    # Iterate through independent variables and calculate beta coefficients
    for col, i in enumerate(fit.model.data.xnames):
        beta = coefficients[i] * (sd_x[col] / sd_y)
        beta_coefficients[i]= beta
        
        conf_int_norm.loc[i,[0,1]] = conf_int.loc[i,[0,1]] * (sd_x[col] / sd_y)
    out.params_norm=beta_coefficients
    out.conf_int_norm=conf_int_norm
    # Print beta coefficients
    #for var, beta in beta_coefficients:
    #print(f' {var}: {beta}')
    
    if len(self.pre)>0:
      #print(out.params)
        if "'" in self.case:
            case=re.search("'.*'",self.case)[0][1:-1]
        else:
            case=self.case
        
        out.ve = (out.params[self.case]*np.var(self.data[case]))/np.var(self.data[self.pre]) 
        out.pc=100*(out.params[self.case]/self.data[self.pre].mean())
    return out



def mean_std(nums,range=False,ext=1):
    # pretty print mean, std
    mm=np.mean(nums)
    ss=np.std(nums)

    if range:
        minn=np.min(nums)
        maxn=np.max(nums)
        if ext==1:
            return f"{mm:.1f}" + '±' + f"{ss:.1f}" + ' (' + f"{minn:.0f}" + '-' + f"{maxn:.0f}" + ')'
        else:
            return f"{mm:.0f}" + '±' + f"{ss:.0f}" + ' (' + f"{minn:.0f}" + '-' + f"{maxn:.0f}" + ')'
    else:
        if ext==1:
            return f"{mm:.1f}" + '±(' + f"{ss:.1f}" + ')'
        else:
            return f"{mm:.0f}" + '±(' + f"{ss:.0f}" + ')'
def find_vars(data,var):
    # return any matching vars in biobank dataset
    #return [s for s in data.columns if var in s]
    return [s for s in data.columns if re.search(var,s)]


def table_means_ext(df, show=[], groupby='Case', cols=[], rows=[], var_names=[], statkeys=['sars']):
    """
    Generate a grouped table of statistics and comparisons.

    Args:
        df (DataFrame): The input DataFrame containing the data.
        show (list, optional): List of columns to show in the table. Defaults to an empty list.
        groupby (str, optional): The column to group the data by. Defaults to 'Case'.
        cols (list, optional): List of columns to include in the table. Defaults to an empty list.
        rows (list, optional): List of rows to include in the table. Defaults to an empty list.
        var_names (list, optional): List of variable names. Defaults to an empty list.
        statkeys (list, optional): List of statistic keys. Defaults to ['sars'].

    Returns:
        DataFrame: The grouped table of statistics and comparisons.
    """
  # draw a pretty grouped table of stats and comparisons 
    
    if show:
        cols=show+[groupby]
        df_c=df.loc[:,cols + ['matched_eid']]
    else:
        df_c=df
    
    if len(var_names)==0:
        var_names=show
    
    # create an empty dictionary
    gp = df_c.groupby(groupby)
    all_keys = [key for key, _ in gp]
   
    #base_gp=gp.get_group(basegp)
    #all_keys.remove(basegps)
    
    stat_results = {}
    for a in statkeys:
        stat_results[a]={}
        
    # loop over column_list and execute code explained above
    for group in statkeys:
        group2 = gp.get_group(group)
        group1 = df_c.loc[group2['matched_eid'],:]
        for column in show:
        # itentify categorical
        #if (~np.isnan(pd.to_numeric(list(filter(lambda v: v==v, df[column])), errors='coerce')).any()):
            ids_case = df_c.loc[df_c[groupby]==group,column].dropna().index
            ids_ctr=df_c.loc[ids_case,'matched_eid'].dropna().values
            ids_ctr=df_c.loc[ids_ctr,column].dropna().index
            ids_case=df_c.loc[ids_ctr,'matched_eid']
              
            #ids_ctr = df_c.loc[ids_case,'matched_eid'].dropna().values
            
            ids = np.concatenate([ids_case,ids_ctr])
            
            if not(is_numeric_dtype(df[column])) or is_bool_dtype(df[column]):
                    contingency=pd.crosstab(index=df_c.loc[ids,groupby],columns=df_c.loc[ids,column])
                    chi2=chi2_contingency(contingency)
                    stat_results[group][column]={'statistic':chi2.statistic, 'pvalue':chi2.pvalue, 'df':len(ids_case)}
            else:
                    # add the output to the dictionary 
                    
                    out=stats.ttest_rel(df_c.loc[ids_case,column],df_c.loc[ids_ctr,column],nan_policy='omit')
                    stat_results[group][column] = {'statistic':out.statistic, 'pvalue':out.pvalue, 'df':len(ids_case)}
                 
    
    idx = pd.IndexSlice
    df_desc_full = df_c.groupby(groupby).describe(include= 'all')
    
    df_desc = df_desc_full.loc[idx[:],idx[:,["mean", "std"]]].T
    df_desc.loc[idx[:,["mean"]],idx[:]] = df_desc.loc[idx[:,["mean"]],idx[:]
                                               ].applymap(
                                               lambda x: "{:.2f}".format(x))
    df_desc.loc[idx[:,["std"]],idx[:]] = df_desc.loc[idx[:,["std"]],idx[:]
                                               ].applymap(
                                               lambda x: " ("+"{:.2f}".format(x)+")")

    df_desc.loc[idx[:,["mean"]],idx[:]] = df_desc.loc[idx[:,["mean"]],idx[:]].values+df_desc.loc[idx[:,["std"]],idx[:]].values
    df_desc=df_desc.loc[idx[:,["mean"]],idx[:]].droplevel(1)
    
    for column in show:
         # flipping to least common as often most interpretable number for bools       
    
        if is_bool_dtype(df[column]):
            #df_desc.loc[column,all_keys] = ((~df_desc_full[column]['top']).astype(str)+":"+(df_desc_full[column]['count']-df_desc_full[column]['freq']).astype(str))
            #df_desc.loc[column,all_keys] = (df_desc_full[column]['count']-df_desc_full[column]['freq']).astype(str)
            df_desc.loc[column,all_keys] = df_c.groupby(groupby)[column].sum().astype(str)
            
            
        elif not(is_numeric_dtype(df[column])) :
            #df_desc.loc[column,all_keys] = df_desc_full[column]['freq'].astype(str)
            df_desc.loc[column,all_keys] = (df_desc_full[column]['count']-df_desc_full[column]['freq']).astype(str)
            if bool(df[column].value_counts().index[0]):
                df_desc.loc[column,all_keys] = (df_desc_full[column]['freq']).astype(str)
                
         
            
       
    for group in statkeys:
        results_df = pd.DataFrame.from_dict(stat_results[group],orient='Index')
        results_df.columns = ['statistic','pvalue','df']

        df_desc.loc[:,'p value (' + group + ')']=results_df.pvalue

        df_desc.loc[:,'count (' + group + ')']=(results_df.df.astype('string'))
    
    df_desc=df_desc.drop('matched_eid',axis=0)
 
    df_desc.index=var_names
    
    df_desc.columns.name = None

    return df_desc#, df_desc_full

# remove outliers from data 
def remove_outliers(data, mult=8, repl=np.nan, axis=0, demean=False, pos=False, log=False):
    """
    Remove outliers from a dataset.

    Parameters:
    - data: numpy array or pandas DataFrame
        The input data.
    - mult: int or float, optional (default=8)
        The multiplier for the median absolute deviation (MAD) threshold.
    - repl: int, float, or np.nan, optional (default=np.nan)
        The value to replace the outliers with.
    - axis: int, optional (default=0)
        The axis along which to compute the MAD and apply the outlier removal.
    - demean: bool, optional (default=False)
        Whether to demean the non-outlier values by subtracting their mean.
    - pos: bool, optional (default=False)
        Whether to set negative values to np.nan before computing the MAD.
    - log: bool, optional (default=False)
        Whether to take the logarithm of the data before computing the MAD.

    Returns:
    - out: numpy array or pandas DataFrame
        The data with outliers removed and replaced with the specified value.
    """
    data = data.copy()
    if pos:
        data[data < 0.00000001] = np.nan
    if log:
        data = np.log(data)

    mad = stats.median_abs_deviation(data, nan_policy='omit', axis=axis)
    out = data.copy()
    els = np.abs(data - np.nanmedian(data, axis=axis)) > (mult * mad)
    out[els] = repl
    if demean:
        out[~els] = out[~els] - np.mean(out[~els], axis=axis)

    return out


# remove confounds from data
# calls remove_outliers,
# remove confounds from data, flattening to combine pre and post points for deconfounding if required
def remove_conf(data, data_var, conf_vars=[], demean=False, return_model=False, remove_out=False, out_pref='_cl', flatten=True, suf=True, suf_pre=False, pos=False, pos_conf=False, log=False, naszero=True, int=False):
    """
    Removes confounding variables from the data and performs outlier removal if specified.

    Parameters:
    - data (DataFrame): The input data.
    - data_var (str): The variable to be analyzed and adjusted.
    - conf_vars (list): List of confounding variables to be removed.
    - demean (bool): If True, the adjusted data will be demeaned.
    - return_model (bool): If True, the fitted model will be returned along with the adjusted data.
    - remove_out (bool): If True, outliers will be removed from the data.
    - out_pref (str): Prefix for the adjusted data variable names.
    - flatten (bool): If True, the data will be flattened out pre and post.
    - suf (bool): If True, suffixes will be added to the confounding variable names.
    - suf_pre (bool): If True, the suffixes for pre and post will be '_pre' and '_post' respectively. Otherwise, they will be '-2.0' and '-3.0'.
    - pos (bool): If True, positive outliers will be removed.
    - pos_conf (bool): If True, positive outliers in the confounding variables will be removed.
    - log (bool): If True, the data will be log-transformed before outlier removal.
    - naszero (bool): If True, missing values in the confounding variables will be replaced with zeros.
    - int (bool): If True, the data will be transformed to ranks and then to standard normal distribution.

    Returns:
    - out (DataFrame): The adjusted data.
    - fitted_model (OLS): The fitted model if return_model is True.
    - df (int): The number of rows in the data if return_model is True.
    """

    vars_all=[data_var]+conf_vars
    
    
    if flatten:
        # flatten out pre and post
        if suf_pre:
           suf_pre='_pre'
           suf_post='_post'
        else:
            suf_pre='-2.0'
            suf_post='-3.0'

         
        #vars_all_pre = [ a + suf_pre for a in vars_all ]
        data_flat_pre=data.loc[:,[]]
        data_flat_pre.loc[:,'eid']=data.index
        data_flat_post=data_flat_pre.copy()
        
        for a in conf_vars:
            if suf:
                data_flat_pre.loc[:,a]=data.loc[:,a+suf_pre]
                data_flat_post.loc[:,a]=data.loc[:,a+suf_post]
            else:
                data_flat_pre.loc[:,a]=data.loc[:,a]
                data_flat_post.loc[:,a]=data.loc[:,a]
                
        data_flat_pre.loc[:,data_var]=data.loc[:,data_var+suf_pre]
        data_flat_post.loc[:,data_var]=data.loc[:,data_var+suf_post]    
        
        data_flat_pre.loc[:,'post']=False
        data_flat_post.loc[:,'post']=True

        df = pd.concat([data_flat_pre,data_flat_post],axis=0,ignore_index=True)
        
        out=df.loc[:,[data_var,'post','eid']]
        
    else:
        df = data.loc[:,vars_all]
        out=df.loc[:,[data_var]]
    # initial removal of outliers   
    if remove_out:
        df[data_var] = remove_outliers(df[data_var],pos=pos,log=log)
        df[conf_vars] = remove_outliers(df[conf_vars],mult=15,repl=0,demean=True,pos=pos_conf)
        
    # remove nas from conf (to retain
    if naszero:
        df[conf_vars]=df[conf_vars].fillna(0)
    
    if int:
        df[data_var] =stats.rankdata(df[data_var],nan_policy='omit')
        df[data_var] =stats.norm.ppf(df[data_var]/(len(df[data_var])+1))
        # scale to mean 0 and std 1
        df[data_var]=df[data_var]-df[data_var].mean()
        df[data_var]=df[data_var]/df[data_var].std() 
    
    # make model
    if len(conf_vars)>0:
        form = ''.join([ "Q('"+aa+"') +" for aa in conf_vars])[:-1]
        fitted_model = smf.ols(formula=" Q('"+data_var+"') ~ " +  form, data=df).fit()
    
    
        if not demean:
            out[data_var]=(df[data_var] - fitted_model.fittedvalues)+df[data_var].mean()
        else:
            out[data_var]=(df[data_var]-fitted_model.fittedvalues)
    else:
        out[data_var]=df[data_var]
    
    # revert to pre/post if flattened
    if flatten:
        out_pre=data_var+suf_pre+out_pref
        out_post=data_var+suf_post+out_pref

        outf=data.loc[:,[]]
        outf.loc[:,out_pre]=out.loc[out['post']==False,data_var].values
        outf.loc[:,out_post]=out.loc[out['post']==True,data_var].values
        out=outf
    else:
        out=df[data_var] 
        fitted_model=[]
        df=df[data_var].shape[0]
    
    # 
    if not return_model:
        return out
    else:
        return fitted_model,df

# print a full table
def print_full(x):
    pd.set_option('display.max_rows', len(x))
    print(x)
    pd.reset_option('display.max_rows')
    

# memory tests    
import sys
def sizeof_fmt(num, suffix='B'):
    ''' by Fred Cirera,  https://stackoverflow.com/a/1094933/1870254, modified'''
    for unit in ['Mi','Gi',]:
        if abs(num) < 1024.0:
            return "%3.1f %s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f %s%s" % (num, 'Yi', suffix)

# print a full table
def model_table(outputs, assays, model_names, vars_show, beta_dirs=[], beta_vars=[], rownames=[], colnames=[], diff_models=False, FPR='row', alpha=0.05, rem_txt=-6, returndata=False,beta_name_add='',beta_name_only=[],show_norms=True):
    """
    Generate a table of model results for proteomics data.

    Args:
        outputs (dict): A dictionary containing the model outputs.
        assays (list): A list of assay names.
        model_names (list): A list of model names.
        vars_show (list): A list of variable names to be shown in the table.
        beta_dirs (list, optional): A list of beta directions. Defaults to an empty list.
        beta_vars (list, optional): A list of beta variables. Defaults to an empty list.
        rownames (list, optional): A list of row names for the table. Defaults to an empty list.
        colnames (list, optional): A list of column names for the table. Defaults to an empty list.
        diff_models (bool, optional): Flag indicating whether different models are used. Defaults to False.
        FPR (str, optional): The false positive rate calculation method. Defaults to 'row'.
        alpha (float, optional): The significance level. Defaults to 0.05.
        rem_txt (int, optional): The number of characters to remove from assay names. Defaults to -6.
        returndata (bool, optional): Flag indicating whether to return additional data. Defaults to False.
        beta_name_add (str): Additional string to pick out beta variable.

    Returns:
        pandas.DataFrame: The model results table.
    """
    
 
    if len(colnames)==0:
        colnames=vars_show
    if len(rownames)==0:
        rownames=[a[:rem_txt] for a in assays]
    if len(beta_dirs)==0:
        beta_dirs=np.zeros(len(assays))
    if len(beta_vars)==0:
        beta_vars=np.zeros(len(assays))
        
    table2=pd.DataFrame({'Assay':rownames})
    table2=table2.set_index('Assay')
    results_tables={}

    
    pvals=np.empty((len(rownames),len(colnames)))
    pvals_corr=np.empty((len(rownames),len(colnames)))
    sigs=np.empty((len(rownames),len(colnames)))
    betas=np.empty((len(rownames),len(colnames)))
    betas_norm=np.empty((len(rownames),len(colnames)))
    conf_int_l=np.empty((len(rownames),len(colnames)))         
    conf_int_l_norm=np.empty((len(rownames),len(colnames)))
    conf_int_u=np.empty((len(rownames),len(colnames)))         
    conf_int_u_norm=np.empty((len(rownames),len(colnames))) 
    pcs=np.empty((len(rownames),len(colnames)))
    resultstables={}
    
    for col in range(len(vars_show)):
        
        resultstable=pd.DataFrame({'Assay':rownames})
        resultstable=resultstable.set_index('Assay',drop=False)

        for row in range(len(assays)): 
        
            if diff_models==True: 
                name=assays[row]+'_'+vars_show[col]+'_'+model_names[0]
            elif diff_models==2:
                name=assays[row]+model_names[0]
            else:
                name=assays[row]+model_names[0]

            ll=list(outputs[name].pvalues.index)
            # show n for model. shows final model if diff_model = True
            n=str(int(outputs[name].df_model + outputs[name].df_resid))
            table2.loc[rownames[row],'n']=n
            resultstable.loc[rownames[row],'n']=n
        
     
            if beta_name_only:
                indices = [i for i, s in enumerate(ll) if (beta_name_only in s)]
            else:  
                indices = [i for i, s in enumerate(ll) if ((vars_show[col] in s)) &  ((beta_name_add in s))]

            #print( [s for i, s in enumerate(ll)]  )
                #print(vars_show[col] )
            #if len(indices)>1:
            #    [ print(ll[i]) for i in indices]
            
            betas[row,col]=(outputs[name].params[indices[0]])
            resultstable.loc[rownames[row],'betas']=betas[row,col]
            betas_norm[row,col]=(outputs[name].params_norm[indices[0]])
            resultstable.loc[rownames[row],'betas_norm']=betas_norm[row,col]

            pcs[row,col]=(outputs[name].pc)
            resultstable.loc[rownames[row],'pc']=pcs[row,col]            
            
            resultstable.loc[rownames[row],'ci_l']=(outputs[name].conf_int.iloc[indices[0],0])
            resultstable.loc[rownames[row],'ci_l_norm']=(outputs[name].conf_int_norm.iloc[indices[0],0])
            resultstable.loc[rownames[row],'ci_u']=(outputs[name].conf_int.iloc[indices[0],1])   
            resultstable.loc[rownames[row],'ci_u_norm']=(outputs[name].conf_int_norm.iloc[indices[0],1])
            pvals[row,col]=(outputs[name].pvalues[indices[0]])
            
            # check if 1-sided
            
            if (beta_dirs[row]*beta_vars[col])==(np.sign(betas[row,col])):
                pvals[row,col]=pvals[row,col]/2

            elif  (beta_dirs[row]*beta_vars[col])==(-np.sign(betas[row,col])): 
                pvals[row,col]=0.5+(1-pvals[row,col])/2
                
            resultstable.loc[rownames[row],'pvals']= pvals[row,col]
        
        resultstables[vars_show[col]]=resultstable
    
    if FPR=='row':
        for row in range(len(assays)):
            (sigs[row,:],pvals_corr[row,:]) = statsmodels.stats.multitest.fdrcorrection(pvals[row,:],alpha=alpha,method='indep')   
    elif FPR=='col':
        for col in range(len(vars_show)):
            (sigs[:,col],pvals_corr[:,col]) = statsmodels.stats.multitest.fdrcorrection(pvals[:,col],alpha=alpha,method='indep')   
    elif FPR=='all':        
        (sigs,pvals_corr) = statsmodels.stats.multitest.fdrcorrection(pvals.flatten(),alpha=alpha,method='indep')
        sigs=sigs.reshape(pvals.shape)

    if show_norms:
        betas=betas_norm
        
    for col in range(len(vars_show)):
        for row in range(len(assays)):
            beta_txt= (' ' + '{0:.4f}'.format(betas[row,col]) ) # '{0:.2f}'.format(outpts_flat[a].params[rowvars[b]]) +
            if (sigs[row,col]==1) :
                table2.loc[rownames[row],colnames[col]]= (beta_txt+' (p=' + '{0:.3f}'.format(pvals[row,col]) + '**)') # '{0:.2f}'.format(outpts_flat[a].params[rowvars[b]]) +
            elif (pvals[row,col]<0.05):
                table2.loc[rownames[row],colnames[col]]= (beta_txt+' (p=' + '{0:.3f}'.format(pvals[row,col]) + '*)') # '{0:.2f}'.format(outpts_flat[a].params[rowvars[b]]) +
            else:
                table2.loc[rownames[row],colnames[col]]= (beta_txt+' (p=' + '{0:.3f})'.format(pvals[row,col]) + ' ') # '{0:.2f}'.format(outpts_flat[a].params[rowvars[b]]) +
                
    
    if returndata:
        return(table2,pvals,pvals_corr,sigs,betas,resultstables)
    else:
        return(table2)

def model_plot(inputs, assays, model_name, var_show, alpha=0.05, plot=False):

    output = pd.DataFrame(index=assays)
    
    for dis in assays:
        tmp = inputs[dis + model_name]
        # hypotheses = "Case_bin = 0, Q('Age-3.0_f'):Case_bin=0"
        disease_output.loc[dis, 'pvalues'] = inputs[dis + model_name].pvalues[var_show]/2
        disease_output.loc[dis, 'params'] = inputs[dis + model_name].params[var_show]
        errs = inputs[dis + model_name].conf_int
        errs_min = errs.loc[var_show, 0]
        errs_max = errs.loc[var_show, 1]
        err = (errs_max - errs_min) / 2
        disease_output.loc[dis, 'err'] = err
    [sigs,pvals] = statsmodels.stats.multitest.fdrcorrection(disease_output.loc[:, 'pvalues'], alpha=alpha, method='poscorr')
    disease_output.loc[:, 'sigs']=sigs
    disease_output.loc[:,  'pvals_corr']=pvals
    
    if plot:
        ax = disease_output.plot(y='params', yerr='err', ms='20', marker='.', ls="", legend=None)
        ax.set_xticks(range(len(disease_output.index)))
        
        ax.set_xticklabels(disease_output.index)
        ax.set_xticklabels([label.replace('_', ' ') for label in disease_output.index])
        plt.xticks(rotation=90)
        ax.grid(visible=True)

        # Add a horizontal line at y=0
        ax.axhline(0, color='black', linestyle='-')

        # Add a star on the plot if sigs is True
        for i, sig in enumerate(disease_output['sigs']):
            if sig:
                ax.text(i, disease_output['params'].max()*1.6, '*', fontsize=14, ha='center', va='bottom',color='red', fontweight='bold')

    return disease_output

def checki(y):
    print(y)

