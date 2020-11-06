#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# dsahead python 3.7
import T2 
import sys
import argparse


def main(argv):
    parser = argparse.ArgumentParser(description='Parse input to T2 socket clients')
    parser.add_argument('ip', type=str, default='10.41.0.109', help='ip address of heimdall')
    parser.add_argument('ports', type=str, default='12345', help='ports address of heimdall (comma-delimited list)')

    args = parser.parse_args()
    ip = args.ip
    ports = [int(port) for port in args.ports.split(',')]

    # "10.41.0.109" is ip of corr21
    # ports are usually 12345,12346
    print(f'Running parse_socket to ip {ip} and ports {ports}')
    T2.socket.parse_socket(host=ip, ports=ports, selectcols=['itime', 'idm', 'ibox', 'ibeam'], outputfile="cluster_output", plot_dir=None)

if __name__ == '__main__':
    main(sys.argv)
