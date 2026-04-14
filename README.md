# GC-depth visualization
GC depth visualization is one of the robust methods to identify contaminants in the genome. Since every taxon/species has different GC content, visualizing GC content with sequencing depth information can identify whether a genome has a contaminant or not. Contaminant are always present in low sequencing depth and has different GC content with the host genome. Many published paper have already implemented GC-depth visualization, but none of them have published the script to visualized it. So, here I present a Python script to compute and visualize GC content vs sequencing depth per genomic window.

## 1. Install
This script use 3 external python libraries, including `numpy`, `matplotlib`, and `spicy`. Make sure all of those libraries already installed in your system.
```
pip install matplotlib numpy scipy
```

### (1) Clone this repository
```
git clone https://github.com/dedee95/GC-depth.git
cd GC-depth/
python gc-depth-plot.py -h
```
### (2) Direct download the script from Github
You can also just download the script `gc-depth-plot.py` alone and directly use it as you want. 

After download the script, then you can type the help message (`-h` or `--help`) to see if the script is working.
```
$ python gc-depth-plot.py -h
Compute and visualize GC content vs sequencing depth per genomic window. One of rebust way to identify contamination in the genome.

Usage: gc-depth-plot.py <fasta> <pandepth_output> [options]

Positional arguments:
  fasta                Genome FASTA file (gzipped is also fine)
  pandepth             Pandepth windowed depth file (.win.stat.gz)

Options:
  -h, --help           Show this help message and exit
  -w, --window WINDOW  Window size, must match pandepth -w value (default: 1000)
  -o, --output OUTPUT  Output plot file (.png or .pdf, default: gc-depth.png)
  --log-depth          Use logarithmic scale for the depth axis
  --plot-only TSV      Skip processing, re-plot from an existing combined TSV (from --output-data)
  --output-data FILE   Save merged GC and depth data to this TSV file (can be reused with --plot-only)

Author: Dede Kurniawan
Email: dedekurniawan@genomics.cn or dedearkun2710@gmail.com
```

## 2. Usage step by step 
There are several upstream steps that you must do before running `gc-depth-plot.py`. The main purpose of the initial step is to generate sequencing depth information in a specific window size.
### 2.1 Map raw reads to the genome
Updating...

## 3. Example
Updating...
