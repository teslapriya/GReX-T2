import pytest
import os.path
from T2 import cluster_heimdall

_install_dir = os.path.abspath(os.path.dirname(__file__))

@pytest.fixture(scope="module")
def datasnrs():
    candsfile = os.path.join(_install_dir, 'data/giants.cand')
    data, snrs = cluster_heimdall.parse_candsfile(candsfile, selectcols=['itime', 'idm', 'ibox'])
    return data, snrs


def test_parse(datasnrs):
    data, snrs = datasnrs
    assert len(data)
    assert len(data) == len(snrs)


def test_cluster1(datasnrs):
    data, snrs = datasnrs
    datal = cluster_heimdall.cluster_data(data)

    assert len(datal) == len(data)
    assert len(datal[0]) == 4


def test_peak(datasnrs):
    data, snrs = datasnrs
    datal = cluster_heimdall.cluster_data(data)

    clsnr = cluster_heimdall.get_peak(datal, snrs)
    assert len(clsnr) == 1

    i, snr = clsnr[0]

    assert i == 1376
    assert snr == 117.613
