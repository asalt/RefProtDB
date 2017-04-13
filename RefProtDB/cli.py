import os
import sys
import re
from itertools import zip_longest
from getpass import getuser
from datetime import datetime

import click

from . import _version
from .download_utils import download_all
from .utils import convert_tab_to_fasta, print_msg, sniff_fasta


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

@click.group(name='main')
@click.version_option(__version__)
def cli():
    pass


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.argument('tabfile', type=click.File('r'))
def convert(tabfile):
    """Convert tab file to fasta"""
    outname, _ = os.path.splitext(tabfile.name)
    outname += '.fa'
    with open(outname, 'w') as outfile:
        for record in convert_tab_to_fasta(tabfile):
            outfile.write(record)
            outfile.write('\n')
    print('Successfully created', outname, 'from', tabfile.name)



@cli.command(context_settings=CONTEXT_SETTINGS)
@click.option('--input-fasta', type=click.Path(exists=True),
              default=None, show_default=True,
              help='Pass in the file of a previously downloaded fasta file'
              'Must have protein accession (ref) )'
)
@click.option('-o', '--outdir', type=click.Path(exists=True),
              default='.', show_default=True)
@click.option('--split/--no-split', default=True,
              help=('After downloading multiple organisms, split the resulting combined fasta'
              'file into a separate file for each organism'),
              show_default=True)
@click.version_option(__version__)
@click.argument('taxonids', nargs=-1, type=int)
def download(outdir, split, taxonids, input_fasta):
    if input_fasta:
        sniff_fasta(input_fasta)  # raises ValueError if not valid input

    """Download and format a fasta file from NCBI"""
    orgs = dict()
    for taxon in taxonids:
        if taxon not in ORGANISMS:
            msg = ('The name for {} is not known\n'
                   'Please enter it now : (e.g. Homo_sapiens)').format(taxon)
            name = click.prompt(msg, type=str, err=True)
            orgs[taxon] = name
            continue
        orgs[taxon] = ORGANISMS[taxon]
    if input_fasta:
        taxonid = click.prompt("Enter the approprate taxon id for your input file", err=True)
        orgs[taxonid] = None

    download_all(orgs, outdir=outdir, split=split, input_fasta=input_fasta)
