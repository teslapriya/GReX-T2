import logging as logger
import grex_t2
import sys
import argparse


def main(argv):
    parser = argparse.ArgumentParser(description="Parse input to T2 socket clients")
    parser.add_argument(
        "--ip",
        type=str,
        default="127.0.0.1",
        help="ip address of heimdall",
        required=False,
    )
    parser.add_argument(
        "--ports",
        type=str,
        default="12345",
        help="ports address of heimdall (comma-delimited list)",
        required=False,
    )
    parser.add_argument(
        "--trigger",
        type=bool,
        default=True,
        help="send trigger to dump buffer",
        required=False,
    )
    parser.add_argument(
        "--source_catalog",
        type=str,
        default="/home/liam/software/GReX-T2/data/catalog.txt",
        help="set to identify triggers from sources",
        required=False,
    )
    args = parser.parse_args()
    ip = args.ip
    ports = [int(port) for port in args.ports.split(",")]
    trigger = args.trigger
    source_catalog = args.source_catalog

    print(
        f"Running parse_socket to ip {ip} and ports {ports} with voltage trigger={trigger}"
    )
    logger.info(
        f"Running parse_socket to ip {ip} and ports {ports} with voltage trigger={trigger}"
    )
    grex_t2.socket.parse_socket(
        host=ip,
        ports=ports,
        selectcols=["itime", "idm", "ibox", "ibeam"],
        outroot="/home/ubuntu/T2_output/cluster_output",
        plot_dir=None,
        trigger=trigger,
        source_catalog=source_catalog,
    )


if __name__ == "__main__":
    main(sys.argv)
