Welcome to CATwalk Streamlit Edition!

This is the official repository for the Streamlit app for CATwalk, a gene assembly software made from scratch.
You can find a publicly hosted at [catwalk.streamlit.app](https://catwalk.streamlit.app/). Take note that this host is very lightweight, meaning that functionality is limited by
small amounts of memory. This problem is not present in any locally hosted instances of the app, nor would it be present in a more robust public
hosting of it.

This project was the result of my research for my Master's of Bioinformatics at San Jose State University. It solves the problem of quickly assembling small portions of
a genome when an entire assembly is either infeasible or unnecessary. Uploading a short read library (intended for Illumina short-read libraries) allows the program to 
build out a scaffold from fragment entries in the file. The app includes a random sequence generator and a reference genome validator for testing purposes. 

While it can produce accurate sequences very quickly, it contains little error correction and should not be used as a definitive source for reference sequences. CATwalk is considered
naïve, which means that it leverages an one-at-a-time decision-making process rather than storing all sequences in memory and assembling from there. This means that CATwalk trades occasional
low accuracy for speed.

CATwalk is created and maintained solely by Gavin Anderson.

© 2026 Gavin Anderson
