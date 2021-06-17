import numpy as np
import socket 
from T2 import cluster_heimdall, plotting
import time
from astropy.time import Time
import datetime
from event import names

from dsautils import dsa_store, dsa_syslog, cnf
ds = dsa_store.DsaStore()
logger = dsa_syslog.DsaSyslogger()
logger.subsystem('software')
logger.app('T2')
my_cnf = cnf.Conf()
t2_cnf = my_cnf.get('t2')

def parse_socket(host, ports, selectcols=['itime', 'idm', 'ibox', 'ibeam'], outroot=None, plot_dir=None, trigger=False):
    """ 
    Takes standard MBHeimdall giants socket output and returns full table, classifier inputs and snr tables.
    selectcols will take a subset of the standard MBHeimdall output for cluster. 
    
    host, ports: same with heimdall -coincidencer host:port
    ports can be list of integers.
    selectcol: list of str.  Select columns for clustering. 
    """

    # count of output - separate from gulps
    globct = 0
    
    if isinstance(ports, int):
        ports = [ports]

    assert isinstance(ports, list)

    lastname = names.get_lastname()

    ss = []

    while True:
        if len(ss) != len(ports):
            for port in ports:
                s = socket.socket(socket.AF_INET, socket.SOCK_STREAM) 
                s.bind((host, port))    # assigns the socket with an address
                s.listen(1)             # accept no. of incoming connections
                ss.append(s)

        ds.put_dict('/mon/service/T2service',
                    {"cadence": 60, "time": Time(datetime.datetime.utcnow()).mjd})

        cls = []
        try:
            for s in ss:
                clientsocket, address = s.accept() # stores the socket details in 2 variables
                logger.info(f"Connection from {address} has been established")
                cls.append(clientsocket)
        except KeyboardInterrupt:
            logger.info("Escaping socket connection")
            break

        # read in heimdall socket output  
        logger.info(f"Reading candsfile from {len(cls)} sockets...")
        print(f"Reading candsfile from {len(cls)} sockets...")
        candsfile = ''
        gulps = []
        for cl in cls:
            cf = recvall(cl, 100000000).decode('utf-8')

            gulp, *lines = cf.split('\n')
            try:
                gulp = int(gulp)
#                print(f"received gulp {gulp} with {len(lines)-1} lines")
            except ValueError:
                print(f'Could not get int from this read ({gulp}). Skipping this client.')
                continue

            gulps.append(gulp)
            cl.close()

            if len(lines) > 1:
                if len(lines[0]) > 0:
                    candsfile += '\n'.join(lines)

        print(f'Received gulp_i {gulps}')
        if len(gulps) != len(cls):
            print(f"not all clients are gulping gulp {gulps}. Skipping...")
            continue
                
        if len(set(gulps)) > 1:
            logger.info(f"not all clients received from same gulp: {set(gulps)}. Restarting socket connections.")
            print(f"not all clients received from same gulp: {set(gulps)}. Restarting socket connections.")

            for s in ss:
                s.close()
            ss = []
            for port in ports:
                s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                try:
                    s.bind((host, port))    # assigns the socket with an address
                except OSError:
                    print('socket bind failed.')
                    continue
                s.listen(1)             # accept no. of incoming connections
                ss.append(s)
            continue
        else:
            ds.put_dict('/mon/service/T2gulp',
                        {"cadence": 60, "time": Time(datetime.datetime.utcnow()).mjd})

        if candsfile == '\n' or candsfile == '':  # skip empty candsfile
            continue

        try:
            tab = cluster_heimdall.parse_candsfile(candsfile)
            logger.info(f"Table has {len(tab)} rows")
            if len(tab) == 0:
                continue
            lastname = cluster_and_plot(tab, globct, selectcols=selectcols, outroot=outroot,
                                        plot_dir=plot_dir, trigger=trigger, lastname=lastname)
            globct += 1
        except KeyboardInterrupt:
            logger.info("Escaping parsing and plotting")
            break
        except OverflowError:
            logger.warning("overflowing value. Skipping this gulp...")
            continue


def cluster_and_plot(tab, gulp_i, selectcols=['itime', 'idm', 'ibox', 'ibeam'], outroot=None, plot_dir=None,
                     trigger=False, lastname=None, max_ncl=None):
    """ 
    Run clustering and plotting on read data.
    Can optionally save clusters as heimdall candidate table before filtering and json version of buffer trigger.
    lastname is name of previously triggered/named candidate.
    """

    # TODO: put these in json config file
    min_dm = t2_cnf['min_dm']  # smallest dm in filtering
    max_ibox = t2_cnf['max_ibox']  # largest ibox in filtering
    min_snr = t2_cnf['min_snr']  # smallest snr in filtering
    if max_ncl is None:
        max_ncl = t2_cnf['max_ncl']  # largest number of clusters allowed in triggering
#    target_params = (405., 420., 6.5)  # R67?

    # cluster
    cluster_heimdall.cluster_data(tab, metric='euclidean', allow_single_cluster=True, return_clusterer=False)
    tab2 = cluster_heimdall.get_peak(tab)
#    tab3 = cluster_heimdall.filter_clustered(tab2, min_snr=min_snr, min_dm=min_dm, max_ibox=max_ibox, target_params=target_params)
    tab3 = cluster_heimdall.filter_clustered(tab2, min_snr=min_snr, min_dm=min_dm, max_ibox=max_ibox)

    col_trigger = np.zeros(len(tab2), dtype=int)
    if outroot is not None and len(tab3):
        outputfile = outroot+str(gulp_i)+".json"
        tab4, lastname = cluster_heimdall.dump_cluster_results_json(tab3, outputfile, trigger=trigger,
                                                                    max_ncl=max_ncl, lastname=lastname)
        if tab4 is not None and trigger:
            col_trigger = np.where(tab4 == tab2, 1, 0)  # if trigger, then overload

    # write T2 cluster results
    if outroot is not None and len(tab2):
        tab2['trigger'] = col_trigger
        cluster_heimdall.dump_cluster_results_heimdall(tab2, outroot+str(np.floor(time.time()).astype('int'))+".cand")
        
    return lastname


def recvall(sock, n):
    """
    helper function to receive all bytes from a socket
    sock: open socket
    n: maximum number of bytes to expect. you can make this ginormous!
    """

    data = bytearray()
    while len(data) < n:
        packet = sock.recv(n - len(data))
        if not packet:
            return data
        data.extend(packet)

    return data
