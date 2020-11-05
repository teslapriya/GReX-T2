import numpy as np
import socket 
from T2 import cluster_heimdall, plotting

import dsautils.dsa_syslog as dsl
logger = dsl.DsaSyslogger()
logger.subsystem('software')
logger.app('T2')


def parse_socket(host, ports, selectcols=['itime', 'idm', 'ibox', 'ibeam'], outputfile=None, plot_dir=None):
    """ 
    Takes standard MBHeimdall giants socket output and returns full table, classifier inputs and snr tables.
    selectcols will take a subset of the standard MBHeimdall output for cluster. 
    
    host, ports: same with heimdall -coincidencer host:port
    ports can be list of integers.
    selectcol: list of str.  Select columns for clustering. 
    """

    if isinstance(ports, int):
        ports = [ports]

    assert isinstance(ports, list)

    ss = []
    for port in ports:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM) 
        s.bind((host, port))    # assigns the socket with an address
        s.listen(5)             # accept no. of incoming connections
        ss.append(s)

    gulp_i = 0 # track heimdall gulp number 
    
    while True:
        cls = []
        for s in ss:
            try:
                clientsocket, address = s.accept() # stores the socket details in 2 variables
                logger.info(f"Connection from {address} has been established")
                cls.append(clientsocket)
            except KeyboardInterrupt:
                logger.info("Escaping socket connection")
                break

        gulp_i += 1 
        logger.info(f"gulp {gulp_i}")

        # read in heimdall socket output  
        receiving = all([len(cl.recv(1)) for cl in cls])           # could also set header as code for gulp number?
        
        if not receiving:
            logger.info(f"not all clients are gulping gulp {gulp_i}.")
            print(f"not all clients are gulping gulp {gulp_i}.")
        else:
            print(f"receiving from all ports")
            logger.info("Reading candsfile...")
            candsfile = ''
            for cl in cls:
                ascii_letter = cl.recv(1)           # recieves an alphabet whose ASCII value is the size of the message 
                if len(ascii_letter) == 0:
                    logger.info(f"gulp {gulp_i} has no giant.")
                else: 
                    candsfile += cl.recv(10000000).decode('utf-8')
                    cl.close()
                    candsfile += '\n'
                    print(candsfile)

            try:
                tab, data, snrs = cluster_heimdall.parse_candsfile(candsfile)
                logger.info(f"Table has {len(tab)} rows")
                cluster_and_plot(tab, data, snrs, gulp_i, selectcols=selectcols, outputfile=outputfile, plot_dir=plot_dir)
            except KeyboardInterrupt:
                logger.info("Escaping parsing and plotting")
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
