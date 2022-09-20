#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 11:20:32 2020

@author: liamconnor
"""

import socket 
from astropy.io import ascii 
import numpy as np
import  pandas as pd

from T2 import cluster_heimdall, socket_grex, names

HOST = "127.0.0.1"
PORT = 12346

gulpsize=16384

s=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
#s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
# Create a UDP socket
s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
# Bind the socket to the port
server_address = (HOST, PORT)
#s.close()
s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
s.bind(server_address)
print("Do Ctrl+c to exit the program !!")

n = 0


candsfile = ''
gulp = 0



while True:
    data, address = s.recvfrom(4096)
    candstr = data.decode('utf-8')
    itime = int(candstr.split('\t')[2])
    iff = int(candstr.split('\t')[1])
    candsfile += candstr

    if itime//gulpsize != gulp:
        print("GULP", gulp)
        socket_grex.filter_candidates(candsfile)
        gulp = itime//gulpsize
        continue
        # tab = ascii.read(candsfile, names=col_heimdall,
        #                  guess=True, fast_reader=False,
        #                  format='no_header')

        # cluster_heimdall.cluster_data(tab, metric='euclidean', 
        #                               allow_single_cluster=True, 
        #                               return_clusterer=False)

        # tab2 = cluster_heimdall.get_peak(tab)
        # tab3 = cluster_heimdall.filter_clustered(tab2, 
        #                                         min_snr=min_snr, 
        #                                         min_dm=min_dm, 
        #                                         max_ibox=max_ibox, 
        #                                         max_cntb=max_cntb,
        #                                         max_cntb0=max_cntb0, 
        #                                         max_ncl=max_ncl, 
        #                                         target_params=target_params)

        # itimes = tab3['itime']
        # maxsnr = tab3['snr'].max()
        # imaxsnr = np.where(tab3['snr'] == maxsnr)[0][0]        
        # itime_imax = str(itimes[imaxsnr])
        # mjd = tab3['mjds'][imaxsnr]
        
        # gulp = itime//gulpsize
        # trigger = False
        # outroot = '/home/liam/data/grex/candidates/T2/'
        # lastname = names.get_lastname_grex(outroot)
        # cat = None
        # beam_model = None
        # coords = None
        # snrs = None
        # nbeams_queue = 0
        # prev_trig_time = None
        # min_timedelt = 60.
        # tab3['mjds'] = 59000.00
        # X = cluster_heimdall.dump_cluster_results_json(
        #                                                 tab3,
        #                                                 trigger=trigger,
        #                                                 lastname=lastname,
        #                                                 cat=cat,
        #                                                 coords=coords,
        #                                                 snrs=snrs,
        #                                                 outroot=outroot,
        #                                                 frac_wide=0.0,
        #                                                )
        
exit()
