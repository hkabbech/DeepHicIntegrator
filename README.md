![DeepIntegrativeHiC release](https://img.shields.io/badge/DeepIntegrativeHiC-v0.1-blue.svg)
[![License: GNU](https://img.shields.io/badge/License-GNU-yellow.svg)](https://opensource.org/licenses/gpl-license)
![Python version](https://img.shields.io/badge/python-3-brightgreen.svg)
[![Documentation Status](https://readthedocs.org/projects/deephicintegrator/badge/?version=latest)](https://deephicintegrator.readthedocs.io/en/latest/?badge=latest)


<br>

# DeepHicIntegrator

This tool permits the integration of a Hi-C matrix with one or several histone marks by interpolating in the latent space of an Autoencoder.

## Installation

### Clone the repository
```
git clone https://github.com/kabhel/DeepHicIntegrator.git
cd DeepHicIntegrator
```

### Requirements

1. A **linux** distribution.

2. **Python3** and the following python packages : **tensorflow**, **keras**, **docopt**, **schema**, **pandas**, **numpy**, **scipy**, **matplotlib**, **sklearn**, **cooler**, **hic2cool** and **m2r** (for Sphinx).

```
pip3 install -r requirements.txt
```

3. A Hi-C matrix in `.hic` file format.

Please, download the **GSE63525 HUVEC** genome in order to run the toy example.

```
wget -i ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_HUVEC_combined_30.hic.gz
gunzip GSE63525_HUVEC_combined_30.hic.gz
```

4. One or several histone marks in 2D dimension.

## Run the program

### Toy example

```
./deep_hic_integrator data/hic_matrix/GSE63525_HUVEC_combined_30.hic data/histone_marks/100K/
```

### Get help
```
    Usage:
        ./deep_hic_integrator <HIC_FILE> <HM_PATH> [--resolution INT]
                                                   [--chr_train INT]
                                                   [--chr_test INT]
                                                   [--hist_mark_train STR]
                                                   [--square_side INT]
                                                   [--epochs INT]
                                                   [--batch_size INT]
                                                   [--encoder STR]
                                                   [--decoder STR]
                                                   [--output PATH]
                                                   [--help]

    Arguments:
        <HIC_FILE>                          Path of the Hi-C matrix file (.hic format)
        <HM_PATH>                           Path of the repository containing the histone mark files

    Options:
        -r, INT, --resolution INT           Resolution representing the number of pair-ended reads
                                            spanning between a pair of bins. [default: 25000]
        -a INT, --chr_train INT             Chromosome used to train the autoencoder [default: 1]
        -t INT, --chr_test INT              Chromosome used to test the autoencoder [default: 20]
        -m STR, --hist_mark_train STR       Name of the histone mark used to train the autoencoder
                                            [default: h3k4me3]
        -n INT, --square_side INT           Size N*N of a sub-matrix [default: 20]
        -p INT, --epochs INT                Number of epochs for the training [default: 50]
        -b INT, --batch_size INT            Batch size for the training [default: 64]
        -e STR, --encoder STR               Trained encoder model (H5 format) [default: None]
        -d STR, --decoder STR               Trained decoder model (H5 format) [default: None]
        -o PATH, --output PATH              Output path [default: results/]
        -h, --help                          Show this
```

## Documentation

The documentation is generated with Sphinx and built on ReadTheDocs.


## Author

[Hélène Kabbech](https://github.com/kabhel) : Bioinformatics master student intern at the Medical Center University of Goettingen (Germany)

## License

This project is licensed under the GNU License.