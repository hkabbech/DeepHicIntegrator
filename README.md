![DeepIntegrativeHiC release](https://img.shields.io/badge/DeepIntegrativeHiC-v0.1-blue.svg)
[![License: GNU](https://img.shields.io/badge/License-GNU-yellow.svg)](https://opensource.org/licenses/gpl-license)
![Python version](https://img.shields.io/badge/python-3-brightgreen.svg)
[![Documentation Status](https://readthedocs.org/projects/deephicintegrator/badge/?version=latest)](https://deephicintegrator.readthedocs.io/en/latest/?badge=latest)


<br>

# DeepHicIntegrator

## Installation

### Clone the repository
```
git clone https://github.com/kabhel/DeepHicIntegrator.git
cd DeepHicIntegrator
```

### Requirements

1. A **linux** distribution.

2. **Python3** and the following python packages :

**tensorflow-gpu**, **keras**, **docopt**, **schema**, **pandas**, **numpy**, **scipy**, **matplotlib**, **sklearn**, **cooler**, **hic2cool** and **m2r** (for Sphinx).

The next command will install all the required packages. Before running the command, make sure that the first line is uncommented. If you do not have GPUs (or do not want to use them) simply replace **tensorflow-gpu** by **tensorflow**.

```
pip install -r requirements.txt
```

3. A Hi-C matrix in `.hic` file format.

Please, download the **GSE63525 HUVEC** genome in order to run the toy example.

```
wget -i ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_HUVEC_combined_30.hic.gz
gunzip GSE63525_HUVEC_combined_30.hic.gz
mv GSE63525_HUVEC_combined_30.hic data/1_binaries/hic/
rm wget-log
```

## Run the program

### Get help
```
    Usage:
        ./deep_hic_integrator <HIC_FILE> [--resolution INT] [--train INT] [--test INT]
                                         [--square_side INT] [--epochs INT] [--batch_size INT]
                                         [--output PATH]

    Arguments:
        <HIC_FILE>                      Path to the Hi-C matrix file
                                        (.hic format)

    Options:
        -r, INT, --resolution INT       Hi-C matrice resolution to use. [default: 25000]
        -a INT, --train INT             Chromosome for training [default: 1]
        -t INT, --test INT              Chromosome for test [default: 20]
        -n INT, --square_side INT       Size n*n of a sub-matrix [default: 60]
        -e INT, --epochs INT            Number of epochs [default: 50]
        -b INT, --batch_size INT        Size of a batch [default: 128]
        -o PATH, --output PATH          Output path [default: output/]
        -h, --help  
```

### Toy example

```
./deep_hic_integrator data/1_binaries/hic/GSE63525_HUVEC_combined_30.hic -a 10 -t 20 -o results/
```

## Documentation

The documentation is generated with Sphinx and built on ReadTheDocs.


## Authors

- [Hélène Kabbech](https://github.com/kabhel)
- [Eduardo Gade Gusmao](https://github.com/eggduzao)

## License

This project is licensed under the GNU License.