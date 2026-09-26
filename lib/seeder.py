from Bio import SeqIO
import random

def randomseq(fasta, length):
     genome = ""
     gen = SeqIO.parse(fasta, "fasta")
     for i in gen:
          genome = i.seq
     bound = len(genome)-length
     randomcut = random.randint(1, bound)
     randomslice = genome[randomcut:randomcut+length]
     return randomslice
