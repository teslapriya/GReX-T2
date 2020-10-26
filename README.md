[![Build Status](https://travis-ci.com/dsa110/dsa110-T2.svg?branch=master)](https://travis-ci.com/dsa110/dsa110-T2)
[![codecov](https://codecov.io/gh/dsa110/dsa110-T2/branch/master/graph/badge.svg)](https://codecov.io/gh/dsa110/dsa110-T2)

# dsa110-T2

Module for parsing MBHeimdall output and clustering candidates to trigger voltage recording.

## Install

`python setup.py install`

## Dependencies
- python 3.6+
- astropy
- numpy
- hdbscan

## Running it (temporary test setup)

1. Move `giants_1.cand` to local directory
2. install gc/dev version
3. checkout development version
4. python scripts/run_T2_v4.py 
5. run heimdall
`./heimdall -f ~gcchen/fake_data/test_32beams_10p_5w.fil -output_dir ~/ -dm 10 2000 -dm_tol 1.25 -nsamps_gulp 4096 -detect_thresh 8 -boxcar_max 12 -nbeams 32 -v -coincidencer 127.0.0.1:12345`

(Should check that socket output is actually the T2 socket client name/port)

## Test
`pytest`
