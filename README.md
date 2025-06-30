# ConforPSSP

**Generative model for protein secondary structure prediction**






ConforPSSP is a protein secondary structure prediction (PSSP) which is able to produce multiple outputs corrsponding to alternative protein conformers. ConforPSSP demonstrated high Q8 accuracy on CASP datasets (80-85%).
\
\
ConforPSSP is available under the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).



# Installation on Linux
1. Make sure curl, git, and wget commands are already installed on your PC. If not present, you need install them at first. For Ubuntu, type sudo apt -y install curl git wget.

2. Make sure you have installed Cuda compiler driver is 12.2 or later. The model is light-weighted and can be operated on a general-purpose workstation too.

3. Download ConforPSSP repository using git.
   git clone git://project.url.here

4. Unzip models in "ConforPSSP/models/" folder.

5. Merge and unzip splitted "ConforPSSP/db/freq_table.json" file: cat "freq_table.a?" > freq_table.7z; 7z x freq_table.7z

6. Install miniconda, see: https://docs.conda.io/projects/miniconda/en/latest/miniconda-other-installer-links.html

7. Create enviroment and activate it.
   conda env create -f environment.yml
   wait
   conda activate conforpssp

Installation is finished. The enviroment includes the following software and python packages:
1. python v3.11
2. Tensorflow v2.15.0
3. numpy==1.26
4. h5py==3.13.0
5. cuda-12.2 toolkit and related packages 




# Run the test case
```
cd ConforPSSP
python run_confor-pssp.py --fasta_file test.fasta --output_dir ./ --model 1
```


# Flags
The model acceptsthe following flags:
--sequence: Input amino acid sequence.
--fasta_file: Name of a file with input amino acid sequences.
--output_dir: Directory where to save the output files.
--N: Number of secondary strucure predictons per input amino acid sequence. Default is 15.
--Y: Size of a set of PSS tokens to constract a frequency table. Default is 15.
--model: Select a trained model among 3 versions. Models 1 and 2 are prained as descriped in the paper and model 3 is trained with less noise. Default is 1.


