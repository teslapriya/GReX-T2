from astropy.io import ascii
import numpy as np
import hdbscan
import pandas as pd 
import json 
import seaborn as sns
import matplotlib.pyplot as plt 


def parse_candsfile(candsfile, selectcols=['itime', 'idm', 'ibox', 'ibeam']):
    """ Takes standard MBHeimdall giants output and returns full table, classifier inputs and snr tables.
    selectcols will take a subset of the standard MBHeimdall output
    (Can add cleaning here, eventually)
    """
    #tab = pd.read_csv(candsfile, names=['snr', 'if', 'itime', 'mjds', 'ibox', 'idm', 'dm', 'ibeam'], delim_whitespace=True) 
    tab = ascii.read(candsfile, names=['snr', 'if', 'itime', 'mjds', 'ibox', 'idm', 'dm', 'ibeam'])
    data = np.lib.recfunctions.structured_to_unstructured(tab[selectcols].as_array())  # ok for single dtype (int)
    snrs = tab['snr']
    # how to use ibeam?

    return tab, data, snrs


def cluster_data(data, min_cluster_size=3, min_samples=5, metric='hamming', returnclusterer=False, allow_single_cluster=True):
    """ Take data from parse_candsfile and identify clusters via hamming metric.
    """

    clusterer = hdbscan.HDBSCAN(metric=metric, min_cluster_size=min_cluster_size,
                                min_samples=min_samples, cluster_selection_method='eom',
                                allow_single_cluster=allow_single_cluster).fit(data) 

    nclustered = np.max(clusterer.labels_ + 1) 
    nunclustered = len(np.where(clusterer.labels_ == -1)[0]) 
    print('Found {0} clustered and {1} unclustered rows'.format(nclustered, nunclustered))
    print('All labels: {0}'.format(np.unique(clusterer.labels_)))
    data_labeled = np.hstack((data, clusterer.labels_[:,None])) 

    if returnclusterer:
        return clusterer
    else:
        return data_labeled


def get_peak(data_labeled, snrs):
    """ Given labeled data, find max snr row per cluster
    """

    clsnr = []
    clusters = data_labeled[:, 3]
    for cluster in np.unique(clusters):
        clusterinds = np.where(cluster == clusters)[0]
        maxsnr = snrs[clusterinds].max()
        imaxsnr = np.where(snrs == maxsnr)[0][0]
        clsnr.append((imaxsnr, maxsnr))

    return clsnr


def dump_cluster_results(tab, clsnr, outputfile, output_cols=['mjds', 'snr', 'ibox', 'dm']):    
    """   
    Takes tab from parse_candsfile and clsnr from get_peak, 
    output columns output_cols into a jason file outputfile. 
    """
    imaxsnr_arr = [clsnr[i][0] for i in range(len(clsnr))] 
    output_dict = {} 
    for i in imaxsnr_arr:
        output_dict[tab['mjds'][i]] = {} 
        for col in output_cols:
            if type(tab[col][i]) == np.int64:
                output_dict[tab['mjds'][i]][col] = int(tab[col][i]) 
            else: 
                output_dict[tab['mjds'][i]][col] = tab[col][i] 

    with open(outputfile, 'w') as f: #encoding='utf-8'
        json.dump(output_dict, f, ensure_ascii=False, indent=4) 
    

def plot_cluster(clusterer, clsnr, snrs, data, tab, cols, plot_dir="./"):
    """
    Plot cluster probabilities in color and greyness. 

    Parameters
    ----------
    clusterer : from cluster_data(). 
    clsnr: from get_peak(). 
    tab, data, snrs: from parse_candsfile(). 
    cols : parameters to plot. 
    plot_dir : where to save output plots. 
    """
    # grey for unclustered noise, pure color for well clustered points. 
    palette = sns.color_palette() 
    cluster_colors = [sns.desaturate(palette[col], sat) if col >= 0 else (0.5, 0.5, 0.5) for col, sat in zip(clusterer.labels_, clusterer.probabilities_)] 
    
    for i in range(len(cols)): 
        fig, ax = plt.subplots() 
        ax.cla() 
        ax.scatter(data[:,i], snrs, s=3, c=cluster_colors) 
        ax.set_xlabel(cols[i]) 
        ax.set_ylabel('snr') 
        ax.set_title('cluster cols:'+str(cols))
        fig.savefig(plot_dir+'snr_'+str(cols[i])+'.pdf') 
        fig.clf()
        plt.close('all')
        
        for j in range(len(cols)): 
            if j > i: 
                fig, ax = plt.subplots() 
                ax.cla()
                ax.scatter(data[:,i],data[:,j], c=cluster_colors)
                ax.set_xlabel(cols[i]) 
                ax.set_ylabel(cols[j]) 
                
                for k in range(len(clsnr)): 
                    imaxsnr = clsnr[k][0]
                    maxsnr = int(clsnr[k][1])
                    ax.scatter(data[:,i][imaxsnr],data[:,j][imaxsnr], s=maxsnr, c='k', marker='*') 
                    ax.text(data[:,i][imaxsnr],data[:,j][imaxsnr], str(maxsnr)) 
                
                ax.set_title('cluster cols:'+str(cols))
                fig.savefig(plot_dir+'cluster_prob_'+cols[i]+'_'+cols[j]+'.pdf') 
                fig.clf()
                plt.close('all') 


def parse_socket():
    pass 