import os 
import string
import datetime
import glob 
from astropy import time


def get_lastname():
    """ Look at etcd to get name of last triggered candidate
    Return of None means that the name generation should start anew.
    """

    from dsautils import dsa_store
    ds = dsa_store.DsaStore()

    try:
        lastname, vv = ds.get_dict('/mon/corr/1/trigger').popitem()
    except:
        lastname = None
        
    return lastname

def get_lastname_grex(outroot):
    files = glob.glob(outroot + '/*.json')

    if len(files):
        return max(fl, key = os.path.getctime)
    else:
        return None

def increment_name(mjd, lastname=None, suffixlength=4):
    """ Use mjd to create unique name for event.
    """

    dt = time.Time(mjd, format='mjd', scale='utc').to_datetime()
    if lastname is None:  # generate new name for this yymmdd
        suffix = string.ascii_lowercase[0]*suffixlength
    else:
        yymmdd = lastname.split('_inj')[0][:-suffixlength]
        dt0 = datetime.datetime(int('20'+yymmdd[0:2]), int(yymmdd[2:4]), int(yymmdd[4:6]))
        if dt.year > dt0.year or dt.month > dt0.month or dt.day > dt0.day:
            # new day, so name starts over
            suffix = string.ascii_lowercase[0]*suffixlength
        else:
            # same day, so increment name
            lastsuffix = lastname.split('_inj')[0][-suffixlength:]
            lastnumber = suffixtonumber(lastsuffix)
            suffix = f'{numbertosuffix(lastnumber+1):a>4}'  # increment name

    newname = f'{str(dt.year)[2:]}{dt.month:02d}{dt.day:02d}{suffix}'
    print(f'Incrementing name from "{lastname}" to "{newname}".')

    return newname

def suffixtonumber(suffix):
    """ Given a set of ascii_lowercase values, get a base 26 number.
    a = 0, ... z = 25, aa = 26, ...
    """

    # int(base=26) doesn't quite work, since first ten ints are actually ints!
    base36 = '0123456789abcdefghijklmnopqrstuvwxyz'
    return int(''.join([base36[base36.index(b)-10] for b in suffix]), 26)


def numbertosuffix(num, base=26, numerals=string.ascii_lowercase):
    """ Given a base=26 number, convert to ascii_lowercase-based name.
    Taken from https://stackoverflow.com/questions/60039572/how-to-increment-alphanumeric-number-in-python
    """

    return ((num == 0) and numerals[0]) or (numbertosuffix(num // base, base, numerals).lstrip(numerals[0]) + numerals[num % base])