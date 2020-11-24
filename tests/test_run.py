import T2
import pytest
import os.path

_install_dir = os.path.abspath(os.path.dirname(__file__))

def test_T2():
    """
    Test T2 functions. Read in heimdall giants, cluster them, send results to outputfile and save plots. 

    Parameters
    ----------
    mode : String. "file" or "socket". 
    outputfile : default is "T2_output.txt".
    candsfile : default is "giants.cand.".
    host, ports : use the same host, port with heimdall -coincidencer HOST:PORT 
    plot : default is True.
    plot_dir : default is "./".
    """

    outputfile = os.path.join(_install_dir, 'data/T2_output.txt')
    candsfile = os.path.join(_install_dir, 'data/giants_1.cand')

    # read in giants 
    tab, data, snrs = T2.cluster_heimdall.parse_candsfile(candsfile, selectcols=['itime', 'dm', 'ibox', 'ibeam'])

    # T2 cluster 
    clusterer, data_labeled = T2.cluster_heimdall.cluster_data(data, min_cluster_size=10, min_samples=10,
                                                               metric='euclidean', allow_single_cluster=True,
                                                               return_clusterer=True)
    clsnr = T2.cluster_heimdall.get_peak(data_labeled, snrs) 
    
    # send T2 cluster results to outputfile
    T2.cluster_heimdall.dump_cluster_results_json(tab, clsnr, outputfile, output_cols=['mjds', 'snr', 'ibox', 'dm'])
    
#    if plot: 
#        T2.cluster_heimdall.plot_giants(tab, plot_dir=plot_dir) # plot giants      
#        T2.cluster_heimdall.plot_clustered(clusterer, clsnr, snrs, data, tab, cols=['itime', 'idm', 'ibox'], plot_dir=plot_dir) # plot cluster results  
