# rpExtractSink

RetroPath2 sink generator. From a given SBML file, the tool will extract all the sink molecules and generate a csv file with the InChI structures. For each molecule, the InChI will be get from:

1. The local cache (rrCache), if available
2. The MetaNetX database from MIRIAM URLs (MetaNetX first)

## Input

Required:
* **input_sbml**: (string) Path to the input SBML file

Optional:
* **--remove-dead-end**: (boolean, default: True) Perform FVA (Flux Variability Analysis) evaluation to remove dead end metabolites
* **--compartment-id**: (string, default: 'c') Specify the compartment from which to extract the sink molecules. The default are for MetaNetX files
* **--standalone**: (boolean, default: False) If True, do not retrieve InChI from Internet
* **--cache-dir**: (string, default: None) Path to the cache directory

## Output

* **output_sbml**: (string) Path to the output csv file


## Installation Guide

### Overview

`rpextractsink` depends on `rplibs`, which depends on `cobra`, which requires `python-libsbml`.  
On Apple Silicon (`arm64`) macOS, `python-libsbml` is not available as a native Conda package.

Therefore, installation must be done using an **Intel (`osx-64`) Conda environment under Rosetta**.

---

### General case
```bash
conda install -c conda-forge rpextractsink
```

---

### Apple Silicon macOS (M1/M2/M3)

#### 1. Install Rosetta 2

```bash
softwareupdate --install-rosetta --agree-to-license
```

#### 2. Install

```bash
CONDA_SUBDIR=osx-64 conda install -c conda-forge rpfba
```

Or with mamba:

```bash
CONDA_SUBDIR=osx-64 mamba install -c conda-forge rpfba
```

#### 3. Persist platform setting

```bash
conda config --env --set subdir osx-64
```

#### 4. Verify installation

```bash
python -c "import rpfba; print('rpfba installed successfully')"
```

#### 5. (Optional) Dev installation

```bash
CONDA_SUBDIR=osx-64 conda env create -f environment.yaml
```

---

### Troubleshooting

#### Solver fails on Apple Silicon

Make sure you are using:

```bash
CONDA_SUBDIR=osx-64
```

#### Wrong architecture environment

Check:

```bash
conda config --show subdir
```

Expected output:

```bash
subdir: osx-64
```

## Use

### Function call from Python code
```python
from rr_cache import rrCache
from rpextractsink import genSink

cache = rrCache()
sink = genSink(cache, args.input_sbml)
```


### Run from CLI
```sh
python -m rpextractsink --help
```

## Tests
Test can be run with the following commands:

### Natively
```bash
python -m pytest
```

## Authors

* **Joan Hérisson**
* Thomas Duigou, Melchior du Lac

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
