#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# dsahead python 3.7
import T2 
import sys
import argparse


def main(argv):
    parser = argparse.ArgumentParser(description='Parse input to T2 socket clients')
    parser.add_argument('ip', type=str, help='ip address of heimdall')
    parser.add_argument('ports', type=str, help='ports address of heimdall')

    args = parser.parse_args()
    print(args)

if __name__ == '__main__':
    T2.socket.parse_socket(host="10.41.0.109", ports=[12345], selectcols=['itime', 'idm', 'ibox', 'ibeam'], outputfile="cluster_output", plot_dir=None)
