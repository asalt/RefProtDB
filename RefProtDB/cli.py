import os
import sys
import re
from itertools import zip_longest
from getpass import getuser
from datetime import datetime

import click

from . import _version
from .download_utils import download_all


__author__ = 'Alexander B. Saltzman, Bhoomi Bhatt, Anna Malovannaya'
__copyright__ = _version.__copyright__
__credits__ = ['Alexander B. Saltzman', 'Bhoomi Bhatt', 'Anna Malovannaya']
__license__ = 'MIT'
__version__ = _version.__version__
__maintainer__ = 'Alexander B. Saltzman'
__email__ = 'saltzman@bcm.edu'


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

ORGANISMS = {
    9606  : 'Homo_sapiens',
    10090 : 'Mus_musculus',
}

@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-o', '--outdir', type=click.Path(exists=True),
              default='.', show_default=True)
@click.option('--split/--no-split', default=True,
              help=('After downloading multiple organisms, split the resulting combined fasta'
              'file into a separate file for each organism'),
              show_default=True)
@click.version_option(__version__)
@click.argument('taxonids', nargs=-1, type=int)
def cli(outdir, split, taxonids):
    orgs = dict()
    for taxon in taxonids:
        if taxon not in ORGANISMS:
            msg = ('The name for {} is not known\n'
                   'Please enter it now : (e.g. Homo_sapiens)').format(taxon)
            name = click.prompt(msg, type=str, err=True)
            orgs[taxon] = name
            continue
        orgs[taxon] = ORGANISMS[taxon]
    download_all(orgs, outdir=outdir,split=split)
