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

        # read in heimdall socket output  
        logger.info(f"Reading candsfile from {len(cls)} sockets...")
        candsfile = ''
        gulps = []
        for i, cl in enumerate(cls):
            cf = cl.recv(10000000).decode('utf-8')
            cl.close()
            print("received cf:")
            print(cf)
            gulps.append(cf.split('\n')[0])
            cf = '\n'.join(cf.split('\n')[1:])
            candsfile += cf
            candsfile += '\n'

        print(f'Received gulp_i {gulps}')
        if len(gulps) != len(cls):
            logger.info(f"not all clients are gulping gulp {gulp_i}.")
            print(f"not all clients are gulping gulp {gulp_i}. Skipping...")
            continue
        else:
            print(f"receiving from all ports")

        if len(set(gulps)) > 1:
            logger.info(f"not all clients received from same gulp: {set(gulps)}")
            print(f"not all clients received from same gulp: {set(gulps)}. Skipping")
            continue

        print("Created candsfile from all clients:")
        print(candsfile)
        if candsfile == '\n':  # skip empty candsfile
            continue

        try:
            tab, data, snrs = cluster_heimdall.parse_candsfile(candsfile)
            logger.info(f"Table has {len(tab)} rows")
            if len(tab) == 0 and len(data) == 0 and len(snrs) == 0:
                continue
            assert len(set(gulps)) == 1
            cluster_and_plot(tab, data, snrs, set(gulps).pop(), selectcols=selectcols, outputfile=outputfile, plot_dir=plot_dir)
        except KeyboardInterrupt:
            logger.info("Escaping parsing and plotting")
            break
        except OverflowError:
            logger.warning("overflowing value. Skipping this gulp...")
            continue


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
        cluster_heimdall.dump_cluster_results_heimdall(tab, clsnr, outputfile+str(gulp_i)+".cand")
            
    if plot_dir is not None: 
         plotting.plot_giants(tab, plot_dir=plot_dir+str(gulp_i)+"_") # plot giants      
         plotting.plot_clustered(clusterer, clsnr, snrs, data, tab, cols=['itime', 'idm', 'ibox'], plot_dir=plot_dir+str(gulp_i)+"_") # plot cluster results  
