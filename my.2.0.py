#!/usr/bin/env python

from Bio.Blast.Applications import NcbiblastnCommandline as blast
from Bio.Blast import NCBIXML
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA

from StringIO import StringIO
from subprocess import Popen, PIPE, STDOUT
import sys
import os
import os.path

from optparse import OptionParser

def main( options ):
    align( options.query_fasta, options.subject_fasta, options.ignore_ambig )


# Wether or not ambiguous bases should be called different or not
ignore_ambig = True


def pairwise_align( seq1, seq2 ):
    alignments = pairwise2.align.globalms( seq1, seq2, 1, -2, -5, -2 )
    count = 1
    get_muts( alignments[0][0], alignments[0][1] )

def tcoffee_align( seq1, seq2 ):
    seqs = StringIO()
    seqs.write( seq1.format( 'fasta' ) )
    seqs.write( seq2.format( 'fasta' ) )
    cmd = "t_coffee -in stdin -output fasta -outfile stdout -quiet"
    p = Popen( cmd.split(), stdout=PIPE, stdin=PIPE, stderr=STDOUT )
    stdeo, stdin = p.communicate( input = seqs.getvalue() )
    s = StringIO( stdeo )

    try:
        p = SeqIO.parse( s, 'fasta' )
        ref = p.next()
        s = p.next()
    except Exception:
        return (stdeo, -1)

    mutations = get_muts( ref.seq, s.seq )

    return (mutations,1)

def is_ambig( base ):
    """
        If a given base is ambiguous or not
    """
    return base.upper() not in IUPACUnambiguousDNA.letters

def get_muts( sub_align, query_align ):
    mutations = []
    count = 1
    for q,s in zip( query_align, sub_align ):
        tmp = "Q: %s S: %s Pos: %s" % (q, s, count)

        if not ignore_ambig and (is_ambig( q ) or is_ambig( s )):
            sys.stderr.write( "Skipping ambiguous base\n" )
            sys.stderr.write( tmp + "\n" )
        elif q != s:
            mutations.append( tmp )
        count += 1
    return mutations

def count_seqs( file_name ):
    total_seq = 0
    fh = open( file_name )
    for l in fh:
        if l[0] == '>':
            total_seq += 1
    fh.close()

    return total_seq

def align( query_fasta, subject_fasta, ignore_ambiguous ):
    global ignore_ambig
    ignore_ambig = ignore_ambiguous
    if not ignore_ambig:
        sys.stderr.write( "Ignoring ambiguous bases" )

    refp = SeqIO.parse( subject_fasta, 'fasta' )
    ref = refp.next()

    fhq = SeqIO.parse( query_fasta, 'fasta' )

    total_seq = count_seqs( query_fasta )

    tcount = 1
    for seq in fhq:
        mutations = []
        sys.stderr.write( "Gathering mutations for sequence %s of %s\n" % (tcount, total_seq) )
        mutations, r1 = tcoffee_align( ref, seq )
        if r1 == -1:
            sys.stderr.write( "%s failed to align with tcoffee. Tcoffee output:\n%s" % (seq.id,r2) )
        print "%s: Total mutations: %s" % (seq.description, len( mutations ))
        for m in mutations:
            print m
        tcount += 1

def check_options( options ):
    if options.query_fasta and options.subject_fasta:
        return True
    return False

def get_options( ):
    o = OptionParser()
    o.add_option( '--query', dest='query_fasta', help='Query fasta file' )
    o.add_option( '--subject', dest='subject_fasta', help='Subject fasta file that has only 1 sequence to align all the query sequences against. Most likely a reference gene fasta file.' )
    o.add_option( '--skip-ambiguous', dest='ignore_ambig', action='store_false', default=True, help='Specify this option to skip ambiguous bases so they are not counted' )
    options, args = o.parse_args()

    if not check_options( options ):
        o.print_help()
        sys.exit( -1 )

    return options

if __name__ == '__main__':
    options = get_options()
    main( options )
