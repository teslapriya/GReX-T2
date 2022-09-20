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

from T2 import socket_grex

HOST = "127.0.0.1"
PORT = 12345

# Use roughly 4 seconds as a gulp size
gulpsize=16384

s=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
# Create a UDP socket
s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
# Ensure that you can reconnect 
s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
# Bind the socket to the port
server_address = (HOST, PORT)
s.bind(server_address)

print("Connected to socket %s:%d" % (HOST, PORT))

gulp = 0

while True:
    data, address = s.recvfrom(4096)
    candstr = data.decode('utf-8')
    # Read time sample to keep track of gulp number
    itime = int(candstr.split('\t')[2])
    candsfile += candstr

    # If the gulp number has changed, cluster the current gulp
    if itime//gulpsize != gulp:
        print("GULP", gulp)
        socket_grex.filter_candidates(candsfile)
        gulp = itime//gulpsize
        candsfile = ''
        continue
        
exit()
