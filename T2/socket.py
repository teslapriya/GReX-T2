import numpy as np
import socket 
from T2 import cluster_heimdall, plotting

import dsautils.dsa_syslog as dsl
logger = dsl.DsaSyslogger()
logger.subsystem('software')
logger.app('T2')


def parse_socket(host, port, selectcols=['itime', 'idm', 'ibox', 'ibeam'], outputfile=None, plot_dir=None):
    """ 
    Takes standard MBHeimdall giants socket output and returns full table, classifier inputs and snr tables.
    selectcols will take a subset of the standard MBHeimdall output for cluster. 
    
    host, port: same with heimdall -coincidencer host:port 
    selectcol: list of str.  Select columns for clustering. 
    """

    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM) 
    s.bind((host,port))    # assigns the socket with an address
    s.listen(5)             # accept no. of incoming connections
    
    gulp_i = 0 # track heimdall gulp number 
    
    while True:
        try:
            clientsocket, address = s.accept() # stores the socket details in 2 variables
        except KeyboardInterrupt:
            logger.info("Escaping socket connection")
            break

        logger.info(f"Connection from {address} has been established")
        gulp_i += 1 
        logger.info(f"gulp {gulp_i}")

        # read in heimdall socket output  
        ascii_letter = clientsocket.recv(1)           # recieves an alphabet whose ASCII value is the size of the message 
        
        if len(ascii_letter) == 0:
            logger.info("This gulp has no heimdall giants output.")
        else: 
            logger.info("Reading candsfile...")
            candsfile = clientsocket.recv(10000000).decode('utf-8')
            clientsocket.close()   
            tab, data, snrs = cluster_heimdall.parse_candsfile(candsfile)
            logger.info(f"Table has {len(tab)} rows")
            try:
                cluster_and_plot(tab, data, snrs, gulp_i, selectcols=selectcols, outputfile=outputfile, plot_dir=plot_dir)
            except KeyboardInterrupt:
                logger.info("Escaping plotting")
                break


def cluster_and_plot(tab, data, snrs, gulp_i, selectcols=['itime', 'idm', 'ibox', 'ibeam'], outputfile=None, plot_dir=None):
    """ 
    Run clustering and plotting on read data.
    """

    # cluster 
    clusterer, data_labeled = cluster_heimdall.cluster_data(data, min_cluster_size=10, min_samples=10, metric='euclidean', allow_single_cluster=True, return_clusterer=True)
    clsnr = cluster_heimdall.get_peak(data_labeled, snrs) 
    clsnr = cluster_heimdall.filter_clustered(clsnr)

    # send T2 cluster results to outputfile
    if outputfile is not None:
        cluster_heimdall.dump_cluster_results(tab, clsnr, outputfile+str(gulp_i)+".txt", output_cols=['mjds', 'snr', 'ibox', 'dm'])
            
    if plot_dir is not None: 
         plotting.plot_giants(tab, plot_dir=plot_dir+str(gulp_i)+"_") # plot giants      
         plotting.plot_clustered(clusterer, clsnr, snrs, data, tab, cols=['itime', 'idm', 'ibox'], plot_dir=plot_dir+str(gulp_i)+"_") # plot cluster results  
