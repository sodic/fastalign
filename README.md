MashMap
========================================================================

This repository provides an implentation of an approximate algorithm for computing glogal genome alignments. The algorithm utilizes several heuristic techniques. Most important are [Winnowing](http://www.cs.princeton.edu/courses/archive/spr05/cos598E/bib/p76-schleimer.pdf), [Jaccard index](https://en.wikipedia.org/wiki/Jaccard_index) and [MinHash estimation](https://en.wikipedia.org/wiki/MinHash).
The algorithm is taken from the original paper found [here](https://www.biorxiv.org/content/early/2018/02/18/259986.1). The main building block of the algorithm are described [here](https://www.biorxiv.org/content/early/2017/01/27/103812).

## Dependencies

### Linux
1. gcc 4.8+
2. cmake 3.5+
3. Boost Library ( http://www.boost.org )

## Installation
1. Clone the repository:
   ```sh
   $ git clone https://github.com/sodic/fastalign.git
   ```  
2. Recursively update the submodules:
   ```sh
   $ git submodule update --init --recursive
   ```  
3. Run the building script from the **project root directory** (this will create the `build` directory containing the executable and the rest of the binaries)
   ```sh
   $ scripts/build.sh
   ```  
Scripts for cleaning and rebuilding the project can be found in the folder `scripts` as well.
## Usage

* The program expects two arguments: a query file and a reference file (in that order):
  ```sh
  $ build/fastalign <query_file_name> <reference_file_name>
  ```
  The output is space-delimited with each line consisting of query name, length,
  0-based start, end, strand, target name, length, start, end and mapping nucleotide
  identity.

## Plot
You can use the script `generateDotPlot` for the visualization of the found mappings. The script can be found in the `scripts` directory.

## Example
The folder `examples` provides two example genomes for testing the program. This can be done with the following commands:
```sh
$ build/fastalign examples/klebsiella_pneumoniae_reference.fasta examples/escherichia_coli_reference.fasta
```
This will generate the output file. To plot the graph, you can run:
```sh
$ scripts/generateDotPlot png large fastalign.out
```
