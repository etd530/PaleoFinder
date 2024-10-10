# PaleoFinder

## Installation
`PaleoFinder` has the following depenencies:
| Software/Package | Version |
| ---------------- | ------- |
| `Python`         | 3.12.6  |
| `docopt`         | 0.6.2   |
| `bioptyhon`      | 1.84    |
| `pandas`         | 2.2.3   |
| `numpy`          | 2.1.1   |
| `taxopy`         | 0.13.0  |
| EMBOSS           | 6.6.0   |
| BLAST            | 2.16.0+ |
| Diamond          | 2.1.9   |
| lalign36         | 36.3.8i |

To install the pipeline, simply clone this repository and install the required Python packages through `pip`. You will also need either `BLAST` or `diamond` (or both), the `EMBOSS` suite, and `lalign36` from the [`fasta36` package](https://github.com/wrpearson/fasta36).

In addition, to run the program you will need some blast or diamond protein databases, as well as a taxonomy database. We recommend building a blast or diamond database from the NR database of NCBI and using the [taxdump files](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip) form NCBI for the taxonomy, but you can use other ones if you want to.

Note that for `diamond`, you will also need the [prot.accession2taxid](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz) to be able to use taxonomy IDs with the database. You can build the database as follows:

```
diamond makedb --in ../nr/nr.gz --db nr -p 20 --taxonmap prot.accession2taxid.FULL.gz --taxonnodes nodes.dmp --taxonnames names.dmp
```

## How to run the program
In order to run the pipeline you will need:
- The proteome of an organism you want to screen for, in FASTA format.
- The genome where you want to look for horizontal transfers of the organisms to which the proteome belongs, in FASTA format.

Then the pipeline can be run as:
```
pseudogene_finder.py runall --proteins=<your_proteme_file> --genome=<your_genome_file> --blastp_db=<your_blastp_database> --parent_taxid=<taxid_of_your_target>
```
The `--parent_taxid` flag specifies the taxonomy ID of your organism of interest, which you can find in the [NCBI Taxonomy database](https://www.ncbi.nlm.nih.gov/taxonomy). For example, if you are looking for inerstions of a Nudivirus in a genome, you may use the taxonomy ID of Nudiviridae.
