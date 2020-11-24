#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# dsahead python 3.7
from astropy.io import ascii
from astropy.io.ascii.core import InconsistentTableError
import numpy as np
import hdbscan
import json
import os.path
from dsautils import dsa_store
ds = dsa_store.DsaStore()
import dsautils.dsa_syslog as dsl
logger = dsl.DsaSyslogger()
logger.subsystem('software')
logger.app('T2')


def parse_candsfile(candsfile):
    """ Takes standard MBHeimdall giants output and returns full table, classifier inputs and snr tables.
    (Can add cleaning here, eventually)
    """

    if os.path.exists(candsfile):
        print(f'Candsfile {candsfile} is path, so opening it')
        candsfile = open(candsfile, 'r').read()
    else:
        ncands = len(candsfile.split('\n'))-1
        print(f'Received {ncands} candidates')
#    candsfile = '\n'.join([line for line in candsfile.split('\n') if line.count(' ') == 7])
#    print(f'Received {ncands0} candidates, removed {ncands0-ncands} lines.')

    try:
        tab = ascii.read(candsfile, names=['snr', 'if', 'itime', 'mjds', 'ibox', 'idm', 'dm', 'ibeam'], guess=True, fast_reader=False, format='no_header')
    except InconsistentTableError:
        print('Inconsistent table. Skipping...')
        return ([], [], [])

    tab['ibeam'] = tab['ibeam'].astype(int)
#
#    snrs = tab['snr']
    # how to use ibeam?
   
#    return tab, data, snrs
    return tab

def cluster_data(tab, selectcols=['itime', 'idm', 'ibox', 'ibeam'], min_cluster_size=3, min_samples=5, metric='hamming', return_clusterer=False, allow_single_cluster=True):
    """ Take data from parse_candsfile and identify clusters via hamming metric.
    selectcols will take a subset of the standard MBHeimdall output
    """

    data = np.lib.recfunctions.structured_to_unstructured(tab[selectcols].as_array())  # ok for single dtype (int)
    try:
        clusterer = hdbscan.HDBSCAN(metric=metric, min_cluster_size=min_cluster_size,
                                    min_samples=min_samples, cluster_selection_method='eom',
                                    allow_single_cluster=allow_single_cluster).fit(data) 

        nclustered = np.max(clusterer.labels_ + 1) 
        nunclustered = len(np.where(clusterer.labels_ == -1)[0]) 
        cl = clusterer.labels_
    except ValueError:
        logger.info("Clustering did not run. Each point assigned to unique cluster.")
        cl = np.arange(len(data))
        nclustered = 0
        nunclustered = len(cl)

    print(f'Found {nclustered} clustered and {nunclustered} unclustered rows')
    logger.info('Found {0} clustered and {1} unclustered rows'.format(nclustered, nunclustered))

    # hack assumes fixed columns
    bl = data[:, 3]
    cntb, cntc = np.zeros((len(data), 1), dtype=int), np.zeros((len(data), 1), dtype=int)
    ucl = np.unique(cl)

    for i in ucl:
        ww = np.where(i == cl)
        cntc[ww] = len(ww[0]) 
        ubl = np.unique(bl[ww])
        for j in ubl:
            wwb = np.where(j == bl)
            cntb[wwb] = len(wwb[0]) 

    # append useful metastats to original data
#    data_labeled = np.hstack((data, cl[:,None], cntb, cntc))
    # modifies tab in place
    tab['cl'] = cl.tolist()
    tab['cntc'] = cntc.flatten().tolist()
    tab['cntb'] = cntb.flatten().tolist()

    if return_clusterer:
#        return clusterer, data_labeled
        return clusterer
#    else:
#        return data_labeled



def get_peak(tab):
    """ Given labeled data, find max snr row per cluster
    Adds in count of candidates in same beam and same cluster.
    Note that unclustered candidates are ignored.
    """

#    clsnr = []
#    cl = datal[:, 4].astype(int)   # hack. should really use table.
#    cnt_beam = datal[:, 5].astype(int)
#    cnt_cl = datal[:, 6].astype(int)
    cl = tab['cl'].astype(int)
#    cnt_beam = tab['cntb'].astype(int)
#    cnt_cl = tab['cntc'].astype(int)
    snrs = tab['snr']
    ipeak = []
    for i in np.unique(cl):
        if i == -1:
            continue
        clusterinds = np.where(i == cl)[0]
        maxsnr = snrs[clusterinds].max()
        imaxsnr = np.where(snrs == maxsnr)[0][0]
        ipeak.append(imaxsnr)
#        clsnr.append((imaxsnr, maxsnr, cnt_beam[imaxsnr], cnt_cl[imaxsnr]))
    print(f"Cluster peaks at {ipeak}:\n{tab[ipeak]}")

    return tab[ipeak]


def filter_clustered(tab, min_snr=None, min_dm=None, max_ibox=None, min_cntb=None, max_cntb=None, min_cntc=None,
                     max_cntc=None):
    """ Function to select a subset of clustered output.
    Can set minimum SNR, min/max number of beams in cluster, min/max total count in cluster.
    """

#    for (imaxsnr, snr, cntb, cntc) in clsnr:
    good = [True] * len(tab)

    if min_snr is not None:
        good *= tab['snr'] > min_snr
    if min_dm is not None:
        good *= tab['dm'] > min_dm
    if max_ibox is not None:
        good *= tab['ibox'] < max_ibox
    if min_cntb is not None:
        good *= tab['cntb'] > min_cntb
    if max_cntb is not None:
        good *= tab['cntb'] < max_cntb
    if min_cntc is not None:
        good *= tab['cntc'] > min_cntc
    if max_cntc is not None:
        good *= tab['cntc'] < max_cntc

    #    clsnr_out.append((imaxsnr, snr, cntb, cntc))
    tab_out = tab[good]

    print(f'Filtering from {len(tab)} to {len(tab_out)} candidates.')
    return tab_out


def dump_cluster_results_json(tab, clsnr, outputfile, output_cols=['mjds', 'snr', 'ibox', 'dm', 'ibeam'], trigger=False):
    """   
    Takes tab from parse_candsfile and clsnr from get_peak, 
    output columns output_cols into a jason file outputfile. 
    trigger is bool to update DsaStore to trigger data dump.
    """

    output_dict = {}
    for i, (imaxsnr, maxsnr, cnt_beam, cnt_cl) in enumerate(clsnr):
        output_dict[str(tab['itime'][imaxsnr])] = {}
        for col in output_cols:
            if type(tab[col][imaxsnr]) == np.int64:
                output_dict[str(tab['itime'][imaxsnr])][col] = int(tab[col][imaxsnr])
            else:
                output_dict[str(tab['itime'][imaxsnr])][col] = tab[col][imaxsnr]

        # append fields from clsnr for now
        output_dict[str(tab['itime'][imaxsnr])]['nbeam'] = int(cnt_beam)
        output_dict[str(tab['itime'][imaxsnr])]['ncluster'] = int(cnt_cl)

    with open(outputfile, 'w') as f: #encoding='utf-8'
        print('Writing candidates to json')
        json.dump(output_dict, f, ensure_ascii=False, indent=4) 

    if trigger and len(output_dict):
        print('*Triggering on first candidate*')
        itimes = list(output_dict.keys())
#        for i, dd in enumerate(output_dict.values()):
#            send = True
#            for kk, vv in dd.items():
#                if kk == 'dm':
#                    if vv < dmmin:
#                        send = False
#                elif kk == 'snr':
#                    if vv < snrmin:
#                        send = False
#                elif kk == 'ibox':
#                    if vv > boxmax:
#                        send = False
        itime = (int(itimes[0])-477)*16  # just take first
        ds.put_dict('/cmd/corr/0', {'cmd': 'trigger', 'val': f'{itime}'})


def dump_cluster_results_heimdall(tab, clsnr, outputfile): 
    """   
    Takes tab from parse_candsfile and clsnr from get_peak, 
    output T2-clustered results with the same columns as heimdall.cand into a file outputfile.
    The output is in pandas format with column names in the 1st row.
    """

    imaxsnr = [clsnr[i][0] for i in range(len(clsnr))] 
    cnt_cl =  [clsnr[i][3] for i in range(len(clsnr))] 
    
    output = tab['snr','if','itime', 'mjds','ibox','idm', 'dm'][imaxsnr] 
    output['members'] = cnt_cl 
    output['ibeam'] = tab['ibeam'][imaxsnr]

    output.write(outputfile, format='ascii')
