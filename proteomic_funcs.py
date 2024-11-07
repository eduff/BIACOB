import pickle as pkl
import pandas
import pandas as pd
from pandas.api.types import is_numeric_dtype
from pandas.api.types import is_bool_dtype
import numpy as np
import pickle as pkl 
import os,re
import scipy
import forestplot as fp
scipy.__version__
import matplotlib.pyplot as plt
from pandas.tseries.offsets import DateOffset

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
                                               lambda x: "{:.3f}".format(x))
    df_desc.loc[idx[:,["std"]],idx[:]] = df_desc.loc[idx[:,["std"]],idx[:]
                                               ].applymap(
                                               lambda x: " ("+"{:.3f}".format(x)+")")

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
            pvals_show=pvals[row,col]
            
            if (pvals_show < 0.001) :
                pvals_show="{:.1e}".format(pvals_show)
            else:
                pvals_show="{0:.3f}".format(pvals_show)
            
            if (sigs[row,col]==1) :        
                table2.loc[rownames[row],colnames[col]]= (beta_txt+' (p=' + pvals_show + '**)') # '{0:.2f}'.format(outpts_flat[a].params[rowvars[b]]) +
            elif (pvals[row,col]<0.05):
                table2.loc[rownames[row],colnames[col]]= (beta_txt+' (p=' + pvals_show  + '*)') # '{0:.2f}'.format(outpts_flat[a].params[rowvars[b]]) +
            else:
                table2.loc[rownames[row],colnames[col]]= (beta_txt+' (p=' + pvals_show  + ')') # '{0:.2f}'.format(outpts_flat[a].params[rowvars[b]]) +
                
    
    if returndata:
        return(table2,pvals,pvals_corr,sigs,betas,resultstables)
    else:
        return(table2)

def model_plot(inputs, assays, model_name, var_show, alpha=0.05, plot=False):

    disease_output = pd.DataFrame(index=assays)
    
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
    


def gauss_sm(x):
    return(np.dot(norm.pdf(np.arange(len(x)),loc=5,scale=1),x))


def time_plot(data,IDP,groupby='Case',rows=[],col_names=[],time='Age-3.0_d',rolling="1500d",ylim=[],xlim=[],grps=None,xlabel=[],ylabel=[],title=[],leg=True,ax=[],loc='upper right'):

    if ax == []:
        ax=plt.gca()

    df=data.loc[:,[groupby,time,IDP]]
    df=df.sort_values(by=time,ascending=True)
    df=df.set_index(time)
    dfg=df.groupby(groupby)

    roll = dfg.rolling(rolling,min_periods=1,center=True) 

    #mm = roll.median()
    mm = roll.mean()
    #mm = roll.apply(gauss_sm)
    sem = roll.sem()
    #std = roll.std()/np.sqrt(40)
    qnt_l=roll.quantile(.1)
    qnt_u=roll.quantile(.9)
    
    if grps is None:
        grps=dfg.groups.keys()
    
    for key in grps:
        if type(leg)==dict:
            name=leg[key]
        else:
            name=key
        if type(rolling)==str:
            xs=mm.loc[(key,slice(None)),:].index.get_level_values(1)/ pd.to_timedelta(1, unit='D') /365 
        else:
            xs=mm.loc[(key,slice(None)),:].index.get_level_values(1)   
   
        df= mm.loc[(key,slice(None)),:]

        sns.lineplot(x=xs[~df.index.duplicated()],y=IDP,data=df[~df.index.duplicated()],label=name,ax=ax)
        plt.sca(ax)
        plt.fill_between(xs,mm.loc[(key,slice(None)),IDP]-sem.loc[(key,slice(None)),IDP],mm.loc[(key,slice(None)),IDP]+sem.loc[(key,slice(None)),IDP],alpha=0.2)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    if leg:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[:],labels[:],loc=loc)
    else:
        ax.get_legend().remove()
        
    if ylim:
        ax.set_ylim(ylim)
    if xlim:
        ax.set_xlim(xlim)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    return(ax)

# Function to calculate primary pre->post prediction models for proteomics data.

def calc_pre_post_models(data,els_in,outpts,assays,data_ins=['simp'],models=['modpre'],inter_vars=[],a='post_cl',c='',ext='-3.0'):
    ###
    # Parameters:
    # data (DataFrame): The input data containing proteomics measurements.
    # els : elements of dataframe for analysis
    # outpts (dict): Output dictionary to add results to.
    # assays (list): List of assays to be analyzed, e.g., ['cog_vars_gi', 'adni_mean', 'GeneralHealth'].
    # data_ins (list, optional): List of data reduction methods, default is ['simp'].
    # models (list, optional): List of model types, default is ['modpre'].
    # inter_vars (list, optional): List of variables to add to models for confound/interaction plots.
    # a (str, optional): Suffix for post data, default is 'post_cl'.
    # c (str, optional): Additional suffix, default is ''.
    # ext (str, optional): Suffix for data, default is '-3.0'.
    # # calculate pre and post models for proteomics data
    # outpts - output dictionary to add to
    # assays - e.g.  SIMOA_assays+ OLINK_assays_final ['cog_vars_gi','adni_mean','GeneralHealth']: # diseases: # SIMOA_assays: # vars['park']:
    # ext: suffix for data (0)
    # data_in: data reduction: 'simp' normal, 'red' reduced, 'disease' under age 73, 'PCR' COVID confirmed
    # model: model type: 'modpre' standard, 'modpre_ext' extended, 'modpre_APOE'  APOE type , 'modpre_GFR' GFR 
    # inter_vars - vars to add to models for confound/interaction plots: e.g. 'adni_mean_diff_cl','cog_vars_gi_diff_cl','GeneralHealth_diff_cl'] + vars['pre_comorb']+ ['vac_prior_first_pos']:  #['hearing-2.0']: # ['vac_prior_first_pos'] vars['diseases_pre']:#vars['OLINK_pre']: # vars['pre_comorb']:
    ###
    
    for b in assays:  # ['cog_vars_gi','adni_mean','GeneralHealth']: #+ OLINK_assays_final:  ['cog_vars_gi','adni_mean'] diseases:#  : # [-2:]:#+  ['cog_vars_gi','adni_mean']: #  + diseases SIMOA_assays + OLINK_assays + diseases +  SIMOA_assays  # ['Ab42','Ab40','Ab42/Ab40','GFAP','pTau-181','NfL']:
        
        for data_in in data_ins: # ['PCR']: #,'red']:# 'simp','disease'
            print(b)
            for model in models: # ['modpre','modpre_ext','modpre_APOE','modpre_GFR']: # 'modpre',#_riskfactors']:#,'modpre']: # 'modpre','modpre_ext',
                IDP=b+"_"+c+a # " regPl_
                IDP_pre=b+"_" +c+'pre_cl'#"cl" # " regPl_
                IDP_diff=b+"_" +c+'diff_cl'
                # matched data only (no missing for IDP)
                els=els_in.copy()
                all_case=els&(data.loc[:,'Case_bin']==1)
                all_control=els&(data.loc[:,'Case_bin']==0)
                
                # ensure that matched case and controls exist for IDP
                matched=data.loc[:,IDP_diff].dropna().index
                matched=data.loc[data.loc[matched,'matched_eid'].dropna().values,IDP_diff].notna()
                ids=matched[matched].index
                
                els[:]=False
                els[ids]=True 

                ids=(all_case&els)
                ids=ids[ids].index
                
                ids2=(all_control&els)
                ids2=ids2[ids2].index

                idels=(ids==ids) #(np.abs((data.loc[ids,'Activity-3.0'].values-data.loc[ids2,'Activity-3.0'].values))<2000) #&  (data.loc[ids,'Smoking_bin-2.0'].values==data.loc[ids2,'Smoking_bin-2.0'].values) # & (data.loc[ids,'BP_meds'].values==data.loc[ids2,'BP_meds'].values) # (np.abs((data.loc[ids,'Hip/Waist-2.0'].values-data.loc[ids2,'Hip/Waist-2.0'].values))<20) #& 
                els[:]=False

                if data_in =='red':
                    idels= idels&(np.abs(data.loc[ids,'APOE_score'].values-data.loc[ids2,'APOE_score'].values))<2
                elif data_in == 'disease':
                    # special case for OLINK disease risk models
                    idels = idels & (data.loc[ids,'matched_age_mean']<73)
                elif data_in == 'PCR':
                    idels = idels & (data.loc[ids,'COVID']=='COVID') 
                else:
                    idels=(ids==ids)


                els[ids[idels]]=True
                els[ids2[idels]]=True  
                
                for age_f in ['','age_f_']: 
                    base_model =  " Q('"+ IDP + "')   ~   Q('"+ IDP_pre + "') + Q('Age"+ext+"') + Q('22001-0.0')   + Q('assessment_sep')+ Q('assessment_sep^2')   "

                    if model == 'modpre_ext':
                        base_model = base_model + "+ Q('21002-2.0') + Q('Activity-2.0')+Q('709-3.0')+C(KeyWorker)  "

                    if model == 'modpre_APOE':
                        base_model = base_model + "+ (Q('A33vA34')) + (Q('A33vA44')) +(Q('A33vA32')) +  Q('Hip/Waist-2.0') +  C(Q('Diabetes2')) + Q('GeneralHealth-2.0') +  C(Q('Smoking_bin"+ext+"'))"

                    if model == 'modpre_riskfactors':
                        base_model = base_model + "+  Q('Hip/Waist-2.0') +  C(Q('Diabetes2'))"

                    if model == 'modpre_GFR':
                        base_model = base_model + "+  Q('GFR_Cys_pre_cl') "            

                    if age_f == '':
                        case_var="Case_bin"
                        case_var_hosp=" (Case_hosp_bin_only)  + (Case_nohosp_bin_only) "
                    else:
                        base_model=base_model.replace('Age'+'-3.0','Age-3.0_f')
                        case_var="Case_bin*Q('Age-3.0_f')"
                        case_var_hosp=" (Case_hosp_bin_only)*Q('Age-3.0_f')  + (Case_nohosp_bin_only)*Q('Age-3.0_f') "
                    
          
                    if model=='modpre':
                        # sensitivity analyses for comorbidities
                        
                        for inter in inter_vars:#  
                            if data[inter].dtype==bool :
                                final_model=base_model + " +  C(Q('" + inter + "'))*"+case_var
                                #query_var = "C(Q('" + inter + "'))[T.True]"
                                final_model_conf=base_model + " +  C(Q('" + inter + "')) + "+case_var

                                #query_var_conf = "C(Q('" + inter + "'))[T.True]"
                            else:
                                # remove age when assessing age-related-vulnerability
                                if (inter[-2:]=='_f') & (age_f==''): 
                                    final_model=base_model.replace('Age'+ext,inter)
                                    final_model = final_model + " +  Q('" + inter + "'):"+case_var+ " +  "+case_var
                                    final_model_conf = final_model + " +  "+case_var
                                    
                                    #query_var = "Q('" + inter + "')"
                                    #final_model_conf=base_model + " +  Q('" + inter + "') + "+case_var                               
                                else:    
                                    final_model=base_model + " +  Q('" + inter + "')*"+case_var
                                    #query_var = "Q('" + inter + "')"
                                    final_model_conf=base_model + " +  Q('" + inter + "') + "+case_var
                                    #query_var_conf = "Q('" + inter + "')"
                            
                            #print([IDP+'_'+inter +'_int_'+age_f+data_in + '_' + model])
                            outpts[IDP+'_'+inter +'_int_'+age_f+data_in + '_' + model]=  ols_simp(formula=final_model,data=data[els],pre=IDP_pre).fit() # +  Q('"+bbd['gSex']+"')++ C(PlateID_"+a+")
                            outpts[IDP+'_'+inter +'_conf_'+age_f+data_in + '_' + model]=  ols_simp(formula=final_model_conf,data=data[els],pre=IDP_pre).fit() # +  Q('"+bbd['gSex']+"')++ C(PlateID_"+a+")
                            

                    outpts[IDP+'_'+age_f+data_in + '_' + model]=ols_simp(formula=base_model + " + " + case_var ,data=data[els],pre=IDP_pre).fit() # +  Q('"+bbd['gSex']+"')++ C(PlateID_"+a+")
                    outpts[IDP+'_'+age_f+data_in + '_' + model +'_hosp_only']=ols_simp(formula=base_model  + " + " + case_var_hosp ,data=data[els],pre=IDP_pre,case='Case_hosp_bin_only').fit()
                    outpts[IDP+'_'+age_f+data_in + '_' + model +'_basxtime']=ols_simp(formula= base_model +  " + Case_bin +  Q('"+ IDP_pre + "')*Q('assessment_sep')" ,data=data[els],pre=IDP_pre).fit()
                
                outpts[IDP+'_simp_modpre_f']=ols_simp(formula="   Q('"+ IDP + "')   ~   Case_bin*Q('Age-3.0_f')  + Q('assessment_sep')+ Q('assessment_sep^2') +  Q('22001-0.0') +  Q('"+ IDP_pre + "')    ", data=data[els],pre=IDP_pre).fit() # +  Q('"+bbd['gSex']+"')++ C(PlateID_"+a+")
                
    return outpts


def plot_gp(data,IDP,els=[],els_ctr=[],ax=[],time='Age-3.0_d',input=[],scatterplot=False):
    if len(els)==0:
        els=~np.isnan(data[IDP]) & (data['Case']=='sars')
        els_ctr=~np.isnan(data[IDP]) & (data['Case']=='ctr')
    else:
        els=els&data.loc[:,IDP].notnull()
        
    if len(ax)>0:
        ax=plt.gca()
        
    sigma_f, l = data.loc[els,IDP].std()/np.sqrt(2), 1
    sigma_f, l =  15,15
    kernel = GPy.kern.RBF(1, sigma_f, l)
    #print(data[time].dt.days.values[els,None]/365)
    #print(data[time].dt.days.values[els,None]/365)
    model = GPy.models.GPRegression(data[time].dt.days.values[els,None],data.loc[:,IDP].values[els,None],kernel) 
   
    if len(els_ctr)>0:
        kernel = GPy.kern.RBF(1, sigma_f, l)
        model_ctr = GPy.models.GPRegression(data[time].dt.days.values[els_ctr,None],data.loc[:,IDP].values[els_ctr,None],kernel) 

    #ax=sns.scatterplot(data=data.loc[els,:],y=IDP,x='Age-3.0',hue='Case')
    if len(input)==0:
        input=np.zeros((30,1))
        input[:,0]=np.arange(30)+50
        
    #input[:,1]=False
    ax=plt.plot(input,model.predict(input)[0][:,0],color='red',label='COVID')
    #input[:,1]=True
    qnts=model.predict_quantiles(input,quantiles=[45,55])
    plt.fill_between(input[:,0],qnts[0].flatten(),qnts[1].flatten(),color='red',alpha=0.2)
    
    #ax=sns.lineplot(x=np.arange(30)+50,y=model_ctr.predict(input)[0][:,0],color='blue',)
    if len(els_ctr)>0:
        ax=plt.plot(input,model_ctr.predict(input)[0][:,0],color='blue',label='Control')
    #input[:,1]=True
        qnts=model_ctr.predict_quantiles(input,quantiles=[45,55])
        plt.fill_between(input[:,0],qnts[0].flatten(),qnts[1].flatten(),color='blue',alpha=0.2)

    #ax.sns.f
    #
    ax=plt.gca()

    ax.set_xlim([input.min(),input.max()])
    #ax.set_ylim([1.6,3.8])
    #ax.set_xlabel("Age")
    ax.set_ylabel(IDP)
    #plt.legend()
    
    return(ax,model)