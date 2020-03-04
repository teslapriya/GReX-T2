from astropy.io import ascii
import numpy as np
import hdbscan


def parse_candsfile(candsfile, selectcols=['itime', 'idm', 'ibox']):
    """ Takes standard MBHeimdall output and returns classifier inputs and snr tables.
    selectcols will take a subset of the standard MBHeimdall output
    (Can add cleaning here, eventually)
    """

    tab = ascii.read(candsfile, names=['snr', 'if', 'itime', 'mjds', 'ibox', 'idm', 'dm', 'ibeam'])
    data = np.lib.recfunctions.structured_to_unstructured(tab[selectcols].as_array())  # ok for single dtype (int)
    snrs = data['snr']
    # how to use ibeam?

    return data, snrs


def cluster_data(data, min_cluster_size, min_samples):
    """ Take data from parse_candsfile and identify clusters via hamming metric.
    """

    clusterer = hdbscan.HDBSCAN(metric='hamming', min_cluster_size=min_cluster_size, min_samples=min_samples, cluster_selection_method='eom', allow_single_cluster=True).fit(data) 

    nclustered = np.max(clusterer.labels_ + 1) 
    nunclustered = len(np.where(clusterer.labels_ == -1)[0]) 
    print('Found {0} clustered and {1} unclustered rows'.format(nclustered, nunclustered))
    print('All labels: {0}'.format(np.unique(clusterer.labels_)))
    data_labeled = np.hstack((data, clusterer.labels_[:,None])) 

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
