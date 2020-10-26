#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# dsahead python 3.7
from astropy.io import ascii
import numpy as np
import hdbscan
import json 
import dsautils.dsa_syslog as dsl
logger = dsl.DsaSyslogger()
logger.subsystem('software')
logger.app('T2')


def parse_candsfile(candsfile, selectcols=['itime', 'idm', 'ibox', 'ibeam']):
    """ Takes standard MBHeimdall giants output and returns full table, classifier inputs and snr tables.
    selectcols will take a subset of the standard MBHeimdall output
    (Can add cleaning here, eventually)
    """
    #tab = pd.read_csv(candsfile, names=['snr', 'if', 'itime', 'mjds', 'ibox', 'idm', 'dm', 'ibeam'], delim_whitespace=True) 
    tab = ascii.read(candsfile, names=['snr', 'if', 'itime', 'mjds', 'ibox', 'idm', 'dm', 'ibeam'], guess=True, fast_reader=False)
    tab['ibeam'] = tab['ibeam'].astype(int)
    data = np.lib.recfunctions.structured_to_unstructured(tab[selectcols].as_array())  # ok for single dtype (int)
    snrs = tab['snr']
    # how to use ibeam?
    
    print("table has", len(tab), "rows")

    logger.info("Parsed candsfile")
    return tab, data, snrs


def cluster_data(data, min_cluster_size=3, min_samples=5, metric='hamming', return_clusterer=False, allow_single_cluster=True):
    """ Take data from parse_candsfile and identify clusters via hamming metric.
    """

    clusterer = hdbscan.HDBSCAN(metric=metric, min_cluster_size=min_cluster_size,
                                min_samples=min_samples, cluster_selection_method='eom',
                                allow_single_cluster=allow_single_cluster).fit(data) 

    nclustered = np.max(clusterer.labels_ + 1) 
    nunclustered = len(np.where(clusterer.labels_ == -1)[0]) 
    logger.info('Found {0} clustered and {1} unclustered rows'.format(nclustered, nunclustered))
    logger.info('All labels: {0}'.format(np.unique(clusterer.labels_)))
    data_labeled = np.hstack((data, clusterer.labels_[:,None])) 
    
    if return_clusterer:
        return clusterer, data_labeled
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
    

