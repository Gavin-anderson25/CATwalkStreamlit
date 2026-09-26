from Bio import SeqIO
import re
import streamlit as st

def compare(location, ref, subj):
    if ref[location] == subj[location]:
        return True
    else:
        return False

@st.dialog("Validating Sequence Against Reference...")
def validatesequence(fasta, subject):
    genome = ""
    gen = SeqIO.parse(fasta, "fasta")
    for i in gen:
        genome = i.seq
    genome = str(genome)
    match = False
    i = int(len(subject)*(3/4))
    begin = 0
    while match == False:
        if i%100 == 0:
            st.write(f"Unwinding ({i})...")
        result = re.search(subject[:-i], genome)
        if result != None:
            match = True
            begin = int(result.start())
            st.write("Match found!")
            st.write(f"i-value: {i}")
        else:
            i += 1
    reference = genome[begin:begin+len(subject)]
    count = 0
    for i in range(len(subject)):
        if compare(i, reference, subject) == True:
            count += 1
    st.session_state["validator_output"] = count/len(subject)
