from Bio import SeqIO
import io
import random
import lib.seeder as sd
import lib.CATwalk as ct
import lib.validator as vd
import streamlit as st

def getseedbutton():
    if file_entry:
        generator_c.empty()
        fastatext = io.TextIOWrapper(file_entry, encoding="utf-8")
        length = int(seedsize_entry)
        generated_seed = sd.randomseq(fastatext, length)
        st.session_state["seed"] = generated_seed
    else:
        st.toast("Please upload a reference genome to the generator first")

def performcatwalk():
    if library_entry:
        dialog_box()
    else:
        st.toast("Please upload a short-read sequencing library first")

@st.dialog("Running CATwalk...")
def dialog_box():
    with st.spinner("Extending Seed...", show_time=True):
        librarytext = io.TextIOWrapper(library_entry, encoding="utf-8")
        seedsequence = seedsequence_entry
        productlength = int(productsize)
        barb_size = int(barbsize)
        errorvariable = float(errorrate)
        cutoffvariable = float(errorcutoff)
        dialog_box = st.container(border=True, height=300, gap="xxsmall", autoscroll=True)
        with dialog_box:
            CATwalk_product = ct.CATwalk(seedsequence, librarytext, productlength, barb_size, errorvariable, cutoffvariable)
    st.session_state["cwproduct"] = CATwalk_product
    st.success("Finished Extension")
    
def performvalidation():
    if product_entry:
        referencetext = io.TextIOWrapper(reference_entry, encoding="utf-8")
        vd.validatesequence(referencetext, product_entry)
        st.session_state["validator_score"] = round(float(st.session_state["validator_output"]), 5)
    else:
        st.toast("Please upload a reference genome to the validator first")

def transfer_to_field(state, load):
    st.session_state[state] = load

st.html("""
  <style>
    [class*="st-key-generator"] {
        background-color: #ffd9b3;
        border-radius: 8px;
        padding: 10px;
    }
    [class*="st-key-extender"] {
        background-color: #e6e6fa;
        border-radius: 8px;
        padding: 10px;
    }
    [class*="st-key-validator"] {
        background-color: #c4ffc4;
        border-radius: 8px;
        padding: 10px;
    }
  </style>
""")

st.title("CATwalk")
st.subheader("A Naïve Gene Assembler by Gavin Anderson")

st.divider()

seed_generator, seed_extender, product_validator = st.columns(3)

with seed_generator:
    with st.container(key="generator"):
        st.header("Seed Generator")
        file_entry = st.file_uploader("**Upload Reference Genome for Seed Generation (.fasta/.fna):**", type=["fasta", "fna"], max_upload_size=200)
        seedsize_entry = st.number_input("**Specify Seed Length (nt):**", value=250)
        st.divider()
        st.button("Generate Random Sequence", on_click=getseedbutton)
        generator_c = st.container(border=True)
        if "seed" in st.session_state:
            generator_c.write(st.session_state["seed"])
            st.button("Transfer to Extender", on_click=transfer_to_field, args=["transfer1", st.session_state["seed"]])

with seed_extender:
    with st.container(key="extender"):
        st.header("Seed Extender")
        library_entry = st.file_uploader("**Upload Library (.fasta/.fna):**", type=["fasta", "fna"], max_upload_size=1000)
        if "transfer1" not in st.session_state:
            st.session_state["transfer1"] = ""
        seedsequence_entry = st.text_input("**Input Seed Sequence (nt):**", value=st.session_state["transfer1"])
        productsize = st.number_input("**Specify Desired Product Length (nt):**", value=2000)
        barbsize = st.number_input("**Barb Length (nt):**", value=15)
        errorrate = st.number_input("**Error Rate (decimal):**", format="%0.3f", value=0.001, step=0.001)
        errorcutoff = st.number_input("**Error Cutoff (decimal):**", format="%0.2f", step=0.01, value=1.05)
        st.divider()
        st.button("Extend Seed with Library Entries", on_click=performcatwalk)
        extender_c = st.container(border=True, height=200, gap="xxsmall")
        if "cwproduct" in st.session_state:
            extender_c.write(st.session_state["cwproduct"])
            st.button("Transfer to Validator", on_click=transfer_to_field, args=["transfer2", st.session_state["cwproduct"]])
            
with product_validator:
    with st.container(key="validator"):
        st.header("Product Validator")
        reference_entry = st.file_uploader("**Upload Reference Genome for Validation(.fasta/.fna):**", type=["fasta", "fna"], max_upload_size=200)
        if "transfer2" not in st.session_state:
            st.session_state["transfer2"] = ""
        product_entry = st.text_area("**Input Extended Product for Validation:**", value=st.session_state["transfer2"])
        st.button("Validate Extended Product Against Reference", on_click=performvalidation)
        st.divider()
        with st.container():
            if "validator_score" in st.session_state:
                st.write(f"Percent Similarity between Extended Product and Reference: {st.session_state["validator_score"]}")

st.space("small")
st.text("\u00A9 Gavin Anderson 2026", text_alignment="center", width="stretch")
