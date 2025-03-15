import numpy as np
import pandas as pd
from scipy.stats import hypergeom
import statsmodels.api as sm
from statsmodels.genmod.families import Binomial
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import multiprocessing
import functools
import warnings

###############################
# Functions
###############################

def filter_parent_reads(x, cutoff=20, baseline=1, test='s'):
    """
    Check if the parent reads pass threshold.
    x: array-like numeric values.
    For test='s' (strict), the maximum value must be >= cutoff and the minimum value >= baseline.
    For test='l' (loose), only the maximum value must be >= cutoff.
    """
    x_sorted = np.sort(x)[::-1]  # sort in descending order
    if test == 's':
        if x_sorted[0] >= cutoff and x_sorted[-1] >= baseline:
            return True
        else:
            return False
    elif test == 'l':
        if x_sorted[0] >= cutoff:
            return True
        else:
            return False
    else:
        raise ValueError("Test values: 's' for strict and 'l' for loose")

def filter_offspring_reads(x, cutoff=20, baseline=1, test='s'):
    """
    Check if offspring reads pass threshold.
    x: array-like with length 48.
    Splits x into 12 replicates, each with 4 elements:
      replicate i = [ x[i], x[i+12], x[i+24], x[i+36] ]  for i=0,...,11.
    For each replicate, the first two elements represent one condition and the last two another.
    For both test modes, at least one replicate must have sum(first two) >= cutoff and
    at least one replicate must have sum(last two) >= cutoff.
    For test 's' the minimum overall x must be >= baseline.
    """
    x = np.array(x)
    if len(x) != 48:
        raise ValueError("Input x must have length 48")
    # Reshape into 12 replicates (rows) with 4 values each (using Fortran order to match R indexing)
    replicates = np.reshape(x, (12, 4), order='F')
    # For each replicate compute sum of first two and sum of last two elements.
    sum_first = replicates[:, 0] + replicates[:, 1]
    sum_last = replicates[:, 2] + replicates[:, 3]
    cond_first = np.any(sum_first >= cutoff)
    cond_last  = np.any(sum_last  >= cutoff)
    if test == 'l':
        return cond_first and cond_last
    elif test == 's':
        if cond_first and cond_last and (np.sort(x)[0] >= baseline):
            return True
        else:
            return False
    else:
        raise ValueError("Test values: 's' for strict and 'l' for loose")

def calc_sim_factor(x):
    """
    Calculate simulation factor based on x which is a two-element array-like.
    If both x[0]>=0 and x[1]>0, returns x[0]/x[1].
    If both are 0, returns 1.
    If x[0]>0 and x[1]==0, returns None.
    """
    x = np.array(x)
    if x[0] >= 0 and x[1] > 0:
        return x[0] / x[1]
    elif x[0] == 0 and x[1] == 0:
        return 1
    elif x[0] > 0 and x[1] == 0:
        return None
    else:
        return None

###############################
# GLM tests on temperature-separated data
###############################

def imprint_GLM_test_on_temp(x, group, temp):
    """
    Test imprinting using a GLM on temperature-specific data.
    x: list or array of counts.
       allele1 is constructed from indices [0:12] and [24:36],
       allele2 from [12:24] and [36:48].
    group: pandas DataFrame containing at least a 'maternal' column.
           Rows of group are selected based on the temperature.
           For temp, the rows used are:
             13: rows [0:3] and [12:15]
             18: rows [3:6] and [15:18]
             23: rows [6:9] and [18:21]
             29: rows [9:12] and [21:24]
    temp: one of 13, 18, 23, or 29.
    Returns the p-value for the maternal effect.
    """
    x = np.array(x)
    if len(x) < 48:
        raise ValueError("Input x must have at least 48 values")
    allele1 = np.concatenate((x[0:12], x[24:36]))
    allele2 = np.concatenate((x[12:24], x[36:48]))
    
    # Determine indices based on temperature (adjusted for 0-indexing)
    if temp == 13:
        idx = np.concatenate((np.arange(0, 3), np.arange(12, 15)))
    elif temp == 18:
        idx = np.concatenate((np.arange(3, 6), np.arange(15, 18)))
    elif temp == 23:
        idx = np.concatenate((np.arange(6, 9), np.arange(18, 21)))
    elif temp == 29:
        idx = np.concatenate((np.arange(9, 12), np.arange(21, 24)))
    else:
        raise ValueError("Temperature must equal one of 13, 18, 23 or 29")
        
    allele1_temp = allele1[idx]
    allele2_temp = allele2[idx]
    # Select corresponding rows from group DataFrame.
    # Here we assume group has enough rows and that the rows correspond to the same order.
    group_temp = group.iloc[idx].copy()
    
    # Prepare endog as a 2-column array for binomial GLM (successes, failures)
    imprint = np.column_stack((allele1_temp, allele2_temp))
    # Add constant to predictor (maternal)
    X = sm.add_constant(group_temp['maternal'])
    # Fit GLM with binomial family. Note: R used quasibinomial, but here we use Binomial.
    try:
        model = sm.GLM(imprint, X, family=Binomial())
        result = model.fit()
    except Exception as e:
        warnings.warn("GLM failed to converge: " + str(e))
        return np.nan
    # Return the p-value for the 'maternal' coefficient.
    # (If you prefer to run an ANOVA test, more work is needed; here we simply report the predictor’s p-value.)
    p_value = result.pvalues.get('maternal', np.nan)
    return p_value

def transReg_GLM_test_on_temp(x, group, eff_lib_size, temp):
    """
    Test for trans‐regulation using a GLM on temperature-specific data.
    x: list/array of counts.
    group: pandas DataFrame containing a 'temperature' column and a 'generation' column.
    eff_lib_size: array-like effective library sizes.
    temp: one of 13, 18, 23, or 29.
    Returns a dict with p-value ('p.transReg') and coefficient ('coef') for the generation effect.
    """
    x = np.array(x)
    eff_lib_size = np.array(eff_lib_size)
    
    if temp == 13:
        allele1_F0_scaled = x[0:3] / eff_lib_size[0:3]
        allele2_F0_scaled = x[12:15] / eff_lib_size[3:6]
        allele1_F1 = np.concatenate((x[24:27], x[48:51]))
        allele2_F1 = np.concatenate((x[36:39], x[60:63]))
    elif temp == 18:
        allele1_F0_scaled = x[3:6] / eff_lib_size[0:3]
        allele2_F0_scaled = x[15:18] / eff_lib_size[3:6]
        allele1_F1 = np.concatenate((x[27:30], x[51:54]))
        allele2_F1 = np.concatenate((x[39:42], x[63:66]))
    elif temp == 23:
        allele1_F0_scaled = x[6:9] / eff_lib_size[0:3]
        allele2_F0_scaled = x[18:21] / eff_lib_size[3:6]
        allele1_F1 = np.concatenate((x[30:33], x[54:57]))
        allele2_F1 = np.concatenate((x[42:45], x[66:69]))
    elif temp == 29:
        allele1_F0_scaled = x[9:12] / eff_lib_size[0:3]
        allele2_F0_scaled = x[21:24] / eff_lib_size[3:6]
        allele1_F1 = np.concatenate((x[33:36], x[57:60]))
        allele2_F1 = np.concatenate((x[45:48], x[69:72]))
    else:
        raise ValueError("Temperature must equal one of 13, 18, 23 or 29")
    
    allele1 = np.concatenate((allele1_F0_scaled, allele1_F1))
    allele2 = np.concatenate((allele2_F0_scaled, allele2_F1))
    # Compute expression proportion
    total = allele1 + allele2
    # Avoid division by zero
    with np.errstate(divide='ignore', invalid='ignore'):
        expr = np.true_divide(allele1, total)
        expr[total == 0] = np.nan

    # Subset group for the given temperature.
    group_temp = group[group['temperature'] == temp].copy()
    # Ensure group_temp and expr have matching lengths.
    if len(expr) != len(group_temp):
        raise ValueError("Mismatch in length between expression data and group data for temp {}".format(temp))
    # Prepare predictors (add constant) for GLM. 
    X = sm.add_constant(group_temp['generation'])
    # For a binomial GLM, we supply the total counts as frequency weights.
    try:
        model = sm.GLM(expr, X, family=Binomial(), freq_weights=total)
        result = model.fit()
    except Exception as e:
        warnings.warn("GLM failed to converge: " + str(e))
        return {'p.transReg': np.nan, 'coef': np.nan}
    p_transReg = result.pvalues.get('generation', np.nan)
    coef_trans = result.params.get('generation', np.nan)
    return {'p.transReg': p_transReg, 'coef': coef_trans}

###############################
# Multipanel plotting
###############################

def multiplot(plots, cols=1, layout=None):
    """
    Arranges multiple matplotlib Axes objects in a grid.
    plots: a list of functions that accept an Axes as argument and plot onto it,
           or already-created Axes objects.
    cols: number of columns in the layout.
    layout: a tuple (nrows, ncols). If None, calculated from cols.
    """
    num_plots = len(plots)
    if num_plots == 0:
        return
    if layout is None:
        ncols = cols
        nrows = int(np.ceil(num_plots / cols))
    else:
        nrows, ncols = layout

    fig, axs = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axs = np.array(axs).flatten() if num_plots > 1 else [axs]
    for ax, plot in zip(axs, plots):
        if callable(plot):
            plot(ax)
        else:
            # Assume plot is an Axes object; re-draw its content into ax.
            # This branch can be adjusted based on how plots are provided.
            pass
    # Hide unused axes
    for ax in axs[num_plots:]:
        ax.axis('off')
    plt.tight_layout()
    plt.show()

###############################
# Enrichment test (Hypergeom)
###############################

def HyperGeomFDR(steps, pool, select, GO, nthreads=4):
    """
    Perform enrichment test based on hypergeometric distribution with simulation for FDR.
    
    Arguments:
      steps: total number of simulation rounds.
      pool: list of all genes.
      select: list of selected genes.
      GO: dictionary mapping GO names to lists of genes.
      nthreads: number of parallel threads.
      
    Returns:
      A pandas DataFrame with columns:
         GO_names, GO_in_select, GO_in_pool, Genes_in_GO, P,
         P_adj_Bonf, P_adj_BH, R_obs, R_exp, FDR.
    """
    GO_names = list(GO.keys())
    num_GO = len(GO_names)
    size_pool = len(pool)
    size_select = len(select)
    
    GO_in_select = np.empty(num_GO, dtype=int)
    GO_in_pool   = np.empty(num_GO, dtype=int)
    Genes_in_GO  = np.empty(num_GO, dtype=int)
    P_val        = np.empty(num_GO, dtype=float)
    R_obs        = np.empty(num_GO, dtype=int)
    
    # Compute hypergeometric p-values for each GO term.
    for i, name in enumerate(GO_names):
        GO_i = set(GO[name])
        GO_in_select[i] = len(set(select) & GO_i)
        GO_in_pool[i]   = len(set(pool) & GO_i)
        Genes_in_GO[i]  = len(GO_i)
        # p-value = 1 - CDF(k-1)
        P_val[i] = 1 - hypergeom.cdf(GO_in_select[i] - 1, size_pool, GO_in_pool[i], size_select)
    
    # Calculate observed rank statistic R_obs for each GO term.
    # R_obs[i] = number of GO terms with p-value <= P_val[i]
    for i in range(num_GO):
        R_obs[i] = np.sum(P_val <= P_val[i])
    
    # Round p-values to high precision (mimicking R's rounding)
    P_val_round = np.round(P_val, decimals=15)
    
    # Simulation: split simulation steps among processes.
    steps_per_thread = [steps // nthreads] * nthreads
    for i in range(steps % nthreads):
        steps_per_thread[i] += 1

    pool_args = []
    for steps_local in steps_per_thread:
        pool_args.append((steps_local, pool, select, GO))
    
    # Worker function: see sim_hyperGeom_worker below.
    with multiprocessing.Pool(nthreads) as pool_proc:
        results = pool_proc.starmap(sim_hyperGeom_worker, pool_args)
    
    # Each result is a numpy array of shape (num_GO, steps_local).
    sim_pvals = np.hstack(results)  # shape (num_GO, steps)
    # Compute q-values: for each simulation round, get the minimum p-value across GO terms.
    min_sim_pvals = np.min(sim_pvals, axis=0)
    q_cutoff = np.quantile(min_sim_pvals, 0.05)
    print('FWER Q-value =', q_cutoff)
    
    # For each GO term, compute expected rank (R_exp) from simulation.
    R_exp = np.array([np.sum(np.round(sim_pvals[i, :], 15) <= P_val_round[i]) for i in range(num_GO)])
    R_exp = R_exp / steps
    # Compute FDR: R_exp / R_obs (avoid division by zero)
    FDR = np.where(R_obs > 0, R_exp / R_obs, np.nan)
    
    # Adjust p-values using Bonferroni and Benjamini–Hochberg procedures.
    P_adj_Bonf, _, _, _ = multipletests(P_val, method='bonferroni')
    P_adj_BH, _, _, _ = multipletests(P_val, method='fdr_bh')
    
    df = pd.DataFrame({
        'GO_names': GO_names,
        'GO_in_select': GO_in_select,
        'GO_in_pool': GO_in_pool,
        'Genes_in_GO': Genes_in_GO,
        'P': P_val,
        'P_adj_Bonf': P_adj_Bonf,
        'P_adj_BH': P_adj_BH,
        'R_obs': R_obs,
        'R_exp': R_exp,
        'FDR': FDR
    })
    return df

def sim_hyperGeom_worker(steps, pool_genes, select, GO):
    """
    Worker function for simulation.
    For each simulation round, sample 'select' from pool_genes and compute hypergeom p-values for each GO.
    Returns: numpy array of shape (num_GO, steps) with p-values.
    """
    GO_names = list(GO.keys())
    num_GO = len(GO_names)
    size_pool = len(pool_genes)
    size_select = len(select)
    sim_mat = np.empty((num_GO, steps), dtype=float)
    
    pool_set = set(pool_genes)
    for j in range(steps):
        # Sample without replacement
        rand_select = np.random.choice(pool_genes, size=size_select, replace=False)
        rand_select_set = set(rand_select)
        for i, name in enumerate(GO_names):
            GO_i = set(GO[name])
            count_in_rand = len(rand_select_set & GO_i)
            count_in_pool = len(pool_set & GO_i)
            sim_mat[i, j] = 1 - hypergeom.cdf(count_in_rand - 1, size_pool, count_in_pool, size_select)
    return sim_mat