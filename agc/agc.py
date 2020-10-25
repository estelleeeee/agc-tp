#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Estelle Mariaux"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Estelle Mariaux"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Estelle Mariaux"
__email__ = "estelle.mariaux@hotmail.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

#1-De-duplication en sequence 'complète' ======================

def read_fasta(amplicon_file, minseqlen):
    """Gets the sequence with length superior minseqlen
    """
    if amplicon_file.endswith(".gz"):
        file = gzip.open(amplicon_file, 'rb')
        sequence = b""
        for line in file:
            if line.startswith(b'>'):
                if len(sequence) >= minseqlen:
                    yield sequence.decode('ascii')
                sequence = b""
            else:
                sequence += line.strip()
        yield sequence.decode('ascii')
    else:
        file = open(amplicon_file, "r")
        sequence = ""
        for line in file:
            if line.startswith(">"):
                if len(sequence) >= minseqlen:
                    yield sequence
                sequence = ""
            else:
                sequence += line.strip()
        yield sequence
    file.close()

def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """Gets the sequences repeated at least mincount times
    """
    dico_temp = {}
    for seq in read_fasta(amplicon_file, minseqlen):
        if seq in dico_temp.keys():
            dico_temp[seq] += 1
        else:
            dico_temp[seq] = 1
    for key, value in sorted(dico_temp.items(), key = lambda x: x[1], reverse = True):
        if value >= mincount:
            yield [key, value]

#2-Recherche de séquences chimeriques par approche "de novo"===

def get_chunks(sequence, chunk_size):
    """Creates at least 4 segments of 1 sequence
    """
    seq_len = len(sequence)
    chunks = []
    for k in range(0, seq_len//chunk_size):
        chunks.append(sequence[chunk_size*k: chunk_size*(k+1)])
    if len(chunks) < 4:
        raise ValueError("ERROR : get_chunks, less than 4 segments of sequence")
    return chunks

def cut_kmer(sequence, kmer_size):
    """Creates a generator of kmer
    """
    for k in range(len(sequence) - kmer_size+1):
        yield sequence[k : k+kmer_size]

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    """Add kmer to dictionnary
    """
    for kmer in cut_kmer(sequence, kmer_size):
        if kmer in kmer_dict:
            kmer_dict[kmer].append(id_seq)
        else:
            kmer_dict[kmer] = [id_seq]
    return kmer_dict

def search_mates(kmer_dict, sequence, kmer_size):
    """doc
    """
    return [i[0] for i in Counter([ids for kmer in cut_kmer(sequence, kmer_size) \
        if kmer in kmer_dict for ids in kmer_dict[kmer]]).most_common(8)]

def get_identity(alignment_list):
    """Computes percentage of identity
    """
    nucl_id = 0
    len_seq = len(alignment_list[0])
    for k in range(0, len_seq):
        if alignment_list[0][k] == alignment_list[1][k]:
            nucl_id +=1
    return round(nucl_id/len_seq*100,1)

def get_unique(ids):
    """Creates dict of unique ID
    """
    return {}.fromkeys(ids).keys()

def detect_chimera(perc_identity_matrix):
    """Creates boolean
    """
    id1, id2 = [], []
    sd = 0
    for segment in perc_identity_matrix:
        sd += statistics.stdev(segment)
        id1.append(segment[0])
        id2.append(segment[1])
    sd_mean = sd /len(perc_identity_matrix)
    if (len(get_unique(id1))>2 or len(get_unique(id2))>2) and sd_mean > 5:
        chimera = True
    else:
        chimera = False
    return chimera

def common(lst1, lst2):
    """doc
    """
    return list(set(lst1) & set(lst2))

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """Removes chimera
    """
    kmer_dict = {}
    for sequence, count_seq in dereplication_fulllength(amplicon_file, minseqlen, mincount):
        chunks = get_chunks(sequence, chunk_size)
        chunks_mates = []
        for chunk in chunks:
            chunks_mates.append(search_mates(kmer_dict, chunk, kmer_size))
        #on ne choisit que les 2 premiers chunks
        #common
        #get_identity
    #yield [sequence, count_seq]

#3-Regroupement glouton=======================================

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """Creates OTU list
    """
    otu = []
    for element in chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
        otu.append(element)
    return otu

def write_OTU(OTU_list, output_file):
    """Creates fasta file
    """
    with open(output_file, 'w') as file:
        for i in range(0,len(OTU_list)):
            file.write(">OTU_{} occurrence:{}\n".format(i+1, OTU_list[i][1]))
            for j in range (0, len(OTU_list[i][0]), 80):
                file.write("{}\n".format(OTU_list[i][0][j:j+80]))

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    otu = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount, \
       args.chunk_size, args.kmer_size)

    write_OTU(otu, "../results/otu.fasta")


if __name__ == '__main__':
    main()
