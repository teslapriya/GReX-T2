#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# dsahead python 3.7
import T2 
import sys
import argparse


def main(argv):
    parser = argparse.ArgumentParser(description='Parse input to T2 socket clients')
    parser.add_argument('--ip', type=str, default='10.41.0.126', help='ip address of heimdall', required=False)
    parser.add_argument('--ports', type=str, default='12345,12346,12347,12348', help='ports address of heimdall (comma-delimited list)', required=False)
    parser.add_argument('--trigger', type=bool, default=True, help='send trigger to dump buffer', required=False)
    args = parser.parse_args()
    ip = args.ip
    ports = [int(port) for port in args.ports.split(',')]
    trigger = args.trigger

    print(f'Running parse_socket to ip {ip} and ports {ports} with voltage trigger={trigger}')
    T2.socket.parse_socket(host=ip, ports=ports, selectcols=['itime', 'idm', 'ibox', 'ibeam'], outputfile="output/cluster_output", plot_dir=None, trigger=trigger)

if __name__ == '__main__':
    main(sys.argv)
