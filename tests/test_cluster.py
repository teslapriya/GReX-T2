import pytest
import os.path
from T2 import cluster_heimdall, plotting

_install_dir = os.path.abspath(os.path.dirname(__file__))

@pytest.fixture(scope="module")
def datasnrs():
    candsfile = os.path.join(_install_dir, 'data/giants.cand')
    tab, data, snrs = cluster_heimdall.parse_candsfile(candsfile, selectcols=['itime', 'idm', 'ibox'])
    return tab, data, snrs


def test_parse(datasnrs):
    tab, data, snrs = datasnrs
    assert len(data)
    assert len(data) == len(snrs)


def test_cluster1(datasnrs):
    tab, data, snrs = datasnrs
    datal = cluster_heimdall.cluster_data(data, return_clusterer=False)

    assert len(datal) == len(data)
    assert len(datal[0]) == 4


def test_peak(datasnrs):
    tab, data, snrs = datasnrs
    datal = cluster_heimdall.cluster_data(data, return_clusterer=False)

    clsnr = cluster_heimdall.get_peak(datal, snrs)
    assert len(clsnr) == 1

    i, snr = clsnr[0]

    assert i == 1380
    assert snr == 117.613


def test_json(datasnrs):
    tab, data, snrs = datasnrs
    datal = cluster_heimdall.cluster_data(data, return_clusterer=False)

    clsnr = cluster_heimdall.get_peak(datal, snrs)
    cluster_heimdall.dump_cluster_results(tab, clsnr, 'test.json')

    assert os.path.exists('test.json')


def test_plot_dmhist(datasnrs):
    tab, data, snrs = datasnrs
    plotting.plot_dm_hist(tab)
