import re
from collections import defaultdict

digit = re.compile(r'(\d+)')

def digit_flag(text):
    if text.isdigit():
        return int(text)
    if text.startswith('NP'):
        return 1
    elif text.startswith('XP'):
        return 0
    return text

def sort_key(text):
    return [ digit_flag(x) for x in digit.split(text) ]

class Record:
    __slots__ = ('sequence', 'geneid', 'taxon', 'ref')
    def __init__(self, sequence, geneid=None, taxon=None, ref=None, gi=None, homologene=None,
                 description=None, symbol=None, **kwargs):
        self.sequence = sequence
        self.geneid = geneid
        self.taxon = set([taxon])
        self.ref = defaultdict(lambda  : { 'gi' : {gi},
                                           'homologene' : {homologene},
                                           'description': {description},
                                           'symbol'     : {symbol}}
        )
        self.ref[ref]['gi'].add(gi)
        self.ref[ref]['homologene'].add(homologene)
        self.ref[ref]['description'].add(description)
        self.ref[ref]['symbol'].add(symbol)

    def __repr__(self):
        seq_len = len(self.sequence)
        len_ = min( seq_len, 12 )
        ret = '{}|{}'.format(self.geneid, self.sequence[0:len_])
        if len_ < seq_len:
            ret += '...'
        return ret

    def __hash__(self):
        return hash((self.sequence, self.geneid))

    def __eq__(self, other):
        return (self.sequence, self.geneid) == (other.sequence, other.geneid)

    def update(self, other):
        if not hash(self) == hash(other):
            raise ValueError('Not the same sequences or invalid object')
        self.taxon |= other.taxon
        # self.ref['gi'] |= other.ref['gi']
        # self.ref['homologene'] |= other.ref['homologene']
        # self.ref['description'] |= other.ref['description']
        # self.ref['symbol'] |= other.ref['symbol']

        # self.symbol |= other.symbol
        # self.homologene |= other.homologene
        # self.description |= other.description

        self.ref.update(other.ref)

    def to_fasta(self):
        top_ref_key = sorted(self.ref, key=sort_key, reverse=True)[0]
        gis = self.ref[top_ref_key]['gi']
        homologenes = self.ref[top_ref_key]['homologene']
        descriptions = self.ref[top_ref_key]['description']
        symbols = self.ref[top_ref_key]['symbol']
        out = {
            'top_ref' : top_ref_key,
            'gis' : ';'.join(map(str, gis)),
            'taxons' : ';'.join(map(str, self.taxon)),  # should only be 1
            'homologene' : ';'.join(map(str, homologenes)),
            'gid' : self.geneid,
            'seq' : self.sequence,
            'description' : ';'.join(descriptions)
        }
        fmt = """>geneid|{gid}|ref|{top_ref}|taxon|{taxons}|gi|{gis}|homologene|{homologene}| {description}\n{seq}"""
        return fmt.format(**out)

class Records(dict):
    def __iadd__(self, record):
        if not isinstance(record, Record):
            raise ValueError('Input must be of instance `Record`')
        h = hash(record)
        if h not in self:
            self[h] = record
            return self
        else:
            self[h].update(record)
        return self
