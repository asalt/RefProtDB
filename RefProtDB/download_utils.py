import os
import re
import csv
import tempfile
import ftplib
import gzip as gz
from ftplib import FTP
from datetime import datetime
from warnings import warn

import numpy as np
import pandas as pd

from .utils import fasta_dict_from_file, print_msg
from .containers import Record, Records

CONN='ftp.ncbi.nlm.nih.gov'

@print_msg("Looking for files at", CONN)
def find_file(target, conn='ftp.ncbi.nlm.nih.gov'):
    """Tries the find the correct location of the protein refseq file
    for a given organism (target).
    Input should be the organism name (e.g. Homo_sapiens).
    Returns the location within the ftp connection at
    /genomes/target/protein.fa.gz
    or None if cannot be found.

    Raises ValueError if the search name is ambiguous and a single location
    cannot be found

    """

    def search(dir_):
        "Search function"
        search = [x for x in ftp.nlst(dir_) if target in x]
        if search and len(search) == 1:
            return search[0]
        elif len(search) > 1:
            raise ValueError("Multiple folders found for {}".format(target))
        else:
            return None

    def get_file(dir_):
        """Try to find the protein file"""
        target = 'protein.fa.gz'
        targetdir = os.path.join(dir_, 'protein')
        search = [x for x in ftp.nlst(targetdir) if target in x]
        if search and len(search) == 1:
            return search[0]
        elif len(search) > 1:
            raise ValueError("Multiple folders found for {}".format(target))
        else:
            return None

    ftp = ftplib.FTP(conn)
    ftp.login()
    found = search('/genomes')
    if found:
        file_ = get_file(found)
        if file_:
            return file_
        return None
    return None

@print_msg("Downloading files at", CONN)
def ftp_download(full_file, outdir='.', conn=None):
    """Download necessary files"""
    #NCBI downloads: (1) RefSeq (*protein.faa.gz); (2) homologene.data; (3) gene directory files.
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    _, outname = os.path.split(full_file)
    outname_full = os.path.join(outdir, outname)
    with open(outname_full, 'wb') as outfile:
        ftp.retrbinary('RETR {}'.format(full_file), outfile.write)
    ftp.quit()
    if not os.path.exists(outname_full):
        raise FileNotFoundError('Could not download', outname_full)
    return outname_full

# Unzip Files:
def unzip(inputf, outputf, append=False):
    "unzip a file"
    mode = 'wb'
    if append:
        mode = 'ab'
    with gz.open(inputf, 'rb') as zipf, open(outputf, mode) as outf:
        for line in zipf:
            outf.write(line)

# def get_unzipped_fname(gzfile, identifier):
#     unzipf, _ = os.path.splitext(gzfile)
#     path, name = os.path.split(unzipf)
#     return os.path.join(path, '{}_{}'.format(identifier, name))

def download_refseq(orgname, tmpdir, append=False):
    target = find_file(orgname, conn=CONN)
    if target is None:
        raise ValueError("Could not find RefSeq for {} at {}/genomes".format(orgname, CONN))
    gzfile = ftp_download(target, outdir=tmpdir, conn=CONN)
    # unzippedf = get_unzipped_fname(gzfile, orgname)
    unzippedf, _ = os.path.splitext(gzfile)
    unzip(gzfile, unzippedf, append=append)
    return

@print_msg("Downloading gene2accession and homologene data")
def download_mappings(conn, tmpdir='.'):
    """Download gene2accession and homologene tables"""
    ftp = FTP(conn)
    ftp.login()
    g2a_file = os.path.join(tmpdir, 'gene2accession.gz')
    hgene_file = os.path.join(tmpdir, 'homologene.tab')
    with open(g2a_file, 'wb') as g2a, open(hgene_file, 'wb') as hgene:
        ftp.retrbinary('RETR /gene/DATA/gene2accession.gz', g2a.write)
        ftp.retrbinary('RETR /pub/HomoloGene/current/homologene.data', hgene.write)

    g2a_unzipped = os.path.splitext(g2a_file)[0] + '.tab'
    unzip(g2a_file, g2a_unzipped)
    os.remove(g2a_file)

    return g2a_unzipped, hgene_file

@print_msg("Formatting gene2accession file")
def gene2accession_formatter(tokeep : list, g2afile):

    tokeep = [str(x) for x in tokeep]
    path, _ = os.path.split(g2afile)
    outfile = os.path.join(path, 'g2a_filtered.tab')
    fieldnames = ['TaxonID', 'GeneID', 'ProteinAccession', 'ProteinGI']

    with open(g2afile) as csvfile, open(outfile, 'w') as outf:
        reader = csv.DictReader(csvfile, delimiter='\t')
        writer = csv.DictWriter(outf, delimiter='\t',
                                fieldnames=fieldnames)
        writer.writeheader()
        for row in reader:
            taxonid = row['#tax_id']
            protacc = row['protein_accession.version']
            if taxonid not in tokeep or protacc.startswith('-'):
                continue
            outrow = dict(
                TaxonID          = taxonid,
                GeneID           = row['GeneID'],
                ProteinAccession = protacc,
                ProteinGI        = row['protein_gi'],
                )
            writer.writerow(outrow)
    return outfile

@print_msg("Preparing reference data")
def load_reference_data(g2a_file, hgene_file):
    g2a_df = (pd.read_table(g2a_file)
              .assign(acc_short = lambda x: x.ProteinAccession.str.extract('(\w+)\.',
                                                                           expand=False))
              .sort_values(by=['GeneID', 'ProteinAccession'], ascending=[1,0])
              .drop_duplicates(subset=['GeneID', 'ProteinAccession'])
    )
    hid_df = (pd.read_table(hgene_file, header=None,
                            names=('Homologene', 'TaxonID', 'GeneID',
                                   'Symbol', 'ProteinGI', 'ProteinAccession'))
              .assign(acc_short = lambda x: x.ProteinAccession.str.extract('(\w+)\.', expand=False))
              .filter(regex='[^ProteinAccession]')
    )
    merge = pd.merge(g2a_df, hid_df, on=['TaxonID', 'GeneID', 'ProteinGI', 'acc_short'], how='left')
    return merge

def download_all(orgs : dict, outdir='.', split=True, *args, **kwargs):
    """ orgs is the subdict with the organisms desired
    e.g.: 9606: 'Homo_Sapiens', 10090: 'Mus_musculus'
    """
    today = datetime.now().strftime("%Y_%m_%d")
    fasta_output = 'PyGrouper_{}_refseq_{}.fa'.format('_'.join(map(str, orgs)),
                                                      today)
    fasta_output = os.path.join(outdir, fasta_output)
    with tempfile.TemporaryDirectory() as tmpdir:
        for org in orgs.values():
            download_refseq(org, tmpdir, append=True)  # keep appending to the resulting protein.fa file
        g2a_file, hgene_file = download_mappings(CONN, tmpdir=tmpdir)
        filtered_g2a = gene2accession_formatter(tokeep=orgs.keys(), g2afile=g2a_file)  # a much smaller file

        ref_df = load_reference_data(filtered_g2a, hgene_file)  # protein_accession ->

        fasta_file = os.path.join(tmpdir, 'protein.fa')
        fasta = fasta_dict_from_file(fasta_file)  # a generator that yields a dictionary of records
        gid_hid = (ref_df[['GeneID', 'Homologene']]
                   .dropna(subset=['Homologene'])
                   .assign(Homologene = lambda x: x['Homologene'].astype('int'))
                   .drop_duplicates()
                   .set_index('GeneID')
                   .pipe(lambda x: pd.Series(data=x.Homologene, index=x.index))
        ).to_dict()

        acc_short_regex = re.compile('(\w+)\.')
        mapping = ref_df.groupby('acc_short')
        records = Records()
        for record in fasta:
            acc = record['ref']
            m = acc_short_regex.match(record['ref'])
            if m:
                acc = m.group(1)
            else:
                print("Invalid accession number :", acc)
                continue
            try:
                g2a_recs = mapping.get_group(acc).to_dict(orient='records')
            except KeyError:
                print('No gene 2 accession info for {}'.format(record['ref']))
                continue
            if len(g2a_recs) > 1:
                break
            for g2a_rec in g2a_recs:
                hgene = gid_hid.get(g2a_rec['GeneID'], '-')
                record['geneid'] = g2a_rec['GeneID']
                record['taxon'] = g2a_rec['TaxonID']
                record['ref'] = g2a_rec['ProteinAccession']  # includes the `.X` after the accession
                record['homologene'] = hgene
                r = Record(**record)
                records += r
        saved_records = 0
        print(datetime.now(), ": Writing {} ...".format(fasta_output), end='')
        with open(fasta_output, 'w') as fa:
            for record in records.values():
                fa.write(record.to_fasta())
                fa.write('\n')
                saved_records += 1
        print('done.')
        print(datetime.now(), ': Saved {} records.'.format(saved_records))
        if not split or len(orgs) == 1:
            return
        for org in orgs:
            saved_records = 0
            org_ = str(org)
            org_output = 'PyGrouper_{}_refseq_{}.fa'.format(org_, today)
            org_output = os.path.join(outdir, org_output)
            print(datetime.now(), ": Writing {} ...".format(org_output), end='')
            with open(org_output, 'w') as f:
                for record in fasta_dict_from_file(fasta_output):
                    if record['taxon'] == org_:
                        r = Record(**record)
                        f.write(r.to_fasta())
                        f.write('\n')
                        saved_records += 1
            print('done.')
            print(datetime.now(), ': Saved {} records.'.format(saved_records))
