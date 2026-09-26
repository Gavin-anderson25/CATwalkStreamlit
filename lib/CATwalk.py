import numpy as np
import re
from Bio import SeqIO
import streamlit as st

def reverse_complement(seq):
    rc = []
    for i in range(len(seq)):
        if seq[i] == "A":
            rc.append("T")
        if seq[i] == "C":
            rc.append("G")
        if seq[i] == "G":
            rc.append("C")
        if seq[i] == "T":
            rc.append("A")
    rc = rc[::-1]
    return str(rc)

def spot(a, b):
    if a == b:
        return True
    else:
        return False
    
def compare_seqs(seq, cand, match, iter_count, error_rate=0.001, cutoff=1.05):
    ref = seq[match.seqstart:match.seqend]
    subj = cand[match.candstart:match.candend]
    error_expected = len(ref) * error_rate
    i=0
    found_errors = 0
    while i < len(ref):
        if spot(ref[i], subj[i]) == False:
            found_errors += 1
        i += 1
    ratio = found_errors / error_expected
    if ratio <= cutoff:
        st.write(f"Candidate {iter_count} accepted!")
        return True
    else:
        return False

def spoolgen(obj):
    for i in obj:
        yield i.seq

class matchobj:
    def __init__(self, seqstart, seqend, candstart, candend, rev=False):
        self.seqstart = seqstart
        self.seqend = seqend
        self.candstart = candstart
        self.candend = candend
        self.rev = rev

def hook(seq, cand, barblength=15):
    if seq == cand:
        return None
    barb = seq[0:barblength]
    patt = re.compile(barb)
    match = re.search(patt, cand)
    if match:
        found = matchobj(0, len(cand[match.start():-1]) + 1, match.start(), len(cand), rev=False)
        return found
    cand = cand[::-1]
    match = re.search(patt, cand)
    if match:
        found = matchobj(0, len(cand[match.start():-1]) + 1, match.start(), len(cand), rev=True)
        return found
    return None

def hook_back(seq, cand, barblength=15):
    if seq == cand:
        return None
    barb = seq[-barblength:]
    patt = re.compile(barb)
    match = re.search(patt, cand)
    if match:
        found = matchobj(-match.end(), len(seq), 0, match.end(), rev=False)
        return found
    cand = cand[::-1]
    match = re.search(patt, cand)
    if match:
        found = matchobj(-match.end(), len(seq), 0, match.end(), rev=True)
        return found
    return None

def connect(match, seq, cand):
    if match.seqstart == 0:
        product = cand[0:match.candstart]
        product = product + seq
        return product
    else:
        product = cand[match.candend:]
        product = seq + product
        return product

def CATwalk(center_sequence, fasta, product_length, barbsize, error_variable, cutoff_variable):
    seq = str(center_sequence)
    final_length = int(product_length)
    for i in range(3):
        fasta.seek(0)
        gen = SeqIO.parse(fasta, "fasta")
        spool = spoolgen(gen)
        iter_count = 0
        while len(seq) < final_length:
            iter_count += 1
            if iter_count%250000 == 0:
                st.write(f"{iter_count} entries searched...")
            try:
                candidate = next(spool)
            except StopIteration:
                break
            candidate = str(candidate)
            match_front = hook(seq, candidate, barblength=barbsize)
            if match_front != None:
                if compare_seqs(seq, candidate, match_front, iter_count, error_rate=error_variable, cutoff=cutoff_variable) == True:
                    seq = connect(match_front, seq, candidate)
            match_back = hook_back(seq, candidate, barblength=barbsize)
            if match_back != None:
                if compare_seqs(seq, candidate, match_back, iter_count, error_rate=error_variable, cutoff=cutoff_variable) == True:
                    seq = connect(match_back, seq, candidate)
    if len(seq) < final_length:
        st.write(f"Unable to extend sequence to specified length. Returning Product of length {len(seq)}.")
        return seq
    else:
        st.write(f"Extension successful. Returning Product of length {len(seq)}.")
        return seq
