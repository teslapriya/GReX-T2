from grex_t2 import cluster_heimdall, socket_grex
import os.path

_install_dir = os.path.abspath(os.path.dirname(__file__))


def test_T2():
    """
    Test T2 functions. Read in heimdall giants, cluster them, send results to outputfile and save plots.

    Parameters
    ----------
    mode : String. "file" or "socket".
    outputfile : default is "T2_output.txt".
    candsfile : default is "giants.cand.".
    host, ports : use the same host, port with heimdall -coincidencer HOST:PORT
    plot : default is True.
    plot_dir : default is "./".
    """

    outputfile = os.path.join(_install_dir, "data/T2_output.txt")
    candsfile = os.path.join(_install_dir, "data/giants_1.cand")

    # read in giants
    tab = cluster_heimdall.parse_candsfile(candsfile)

    # T2 cluster
    cluster_heimdall.cluster_data(
        tab,
        min_cluster_size=10,
        min_samples=10,
        selectcols=["itime", "dm", "ibox", "ibeam"],
        metric="euclidean",
        allow_single_cluster=True,
        return_clusterer=True,
    )
    cluster_heimdall.get_peak(tab)

    # send T2 cluster results to outputfile
    row, candname = cluster_heimdall.dump_cluster_results_json(
        tab, outputfile=outputfile, output_cols=["mjds", "snr", "ibox", "ibeam", "dm"]
    )

    assert candname is not None
