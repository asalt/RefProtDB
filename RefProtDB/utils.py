import os
import re
from datetime import datetime
from functools import wraps

def print_msg(*msg, end='...'):
    def deco(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            print(datetime.now(), ':', *msg, end=end, flush=True)
            result = func(*args, **kwargs)
            print('done.')
            return result
        return wrapper
    return deco

# heavily inspired and borrowed
# by https://github.com/jorvis/biocode/blob/master/lib/biocode/utils.py#L149
# Joshua Orvis
def fasta_dict_from_file(fasta_file):
    """
    Reads a file of FASTA entries and returns a dict for each entry.
    It parsers the headers such that `>gi|12345|ref|NP_XXXXX| Description`
    returns as `{gi: 12345, ref: 'NP_XXXXX', description : 'Description'}`
    The sequence has all whitespace removed.
    """
    current_id = dict()
    current_seq = ''
    current_header = None
    pat = re.compile('>(\S+)\s*(.*)')
    header_pat = re.compile(r'(\w+)\|(\w+\.?\w*)')

    def parse_header(header):
        keys = header_pat.findall(header)
        header_data = dict()
        for key in keys:
            header_data[key[0]] = key[1]
            # gi -> ProteinGI #, ref -> NP_XXXX
        return header_data

    for line in open(fasta_file):
        line = line.rstrip()
        m = pat.search(line)
        if m:
            ## new residue line matched, purge the existing one, if not the first
            if current_id:
                ## remove all whitespace and save
                current_seq = ''.join(current_seq.split())
                current_id['sequence'] = current_seq
                yield current_id
                current_id.clear()

            current_seq = ''
            header = m.group(1)
            current_id = parse_header(header)
            current_id['description'] = m.group(2)

        else:
            ## python 2.6+ makes string concatenation amortized O(n)
            ##  http://stackoverflow.com/a/4435752/1368079
            current_seq += str(line)

    ## don't forget the last one
    current_seq = ''.join(current_seq.split())
    current_id['sequence'] = current_seq
    yield current_id
