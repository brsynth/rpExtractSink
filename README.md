# rpExtractSink


RetroPath2 sink generator

## Input

Required:
* **input_sbml**: (string) Path to the input SBML file

Optional:
* **--remove_dead_end**: (boolean, default: True) Perform FVA evaluation to remove dead end metabolites
* **--compartment_id**: (string, default: 'c') Specify the compartment from which to extract the sink molecules. The default are for MetaNetX files

## Output

* **output_sbml**: (string) Path to the output csv file


## Install
```sh
[sudo] conda install -c conda-forge rpextractsink
```

## Use

### Function call from Python code
```python
from rptools.rpextractsink import rpextractsink

sink = rpExtractSink(input_sbml, output_sink)
sink.genSink()
```

If parameters from CLI have to be parsed, the function `build_args_parser` is available:
```python
from rptools.pextractsink import build_args_parser

parser = buildparser()
params = parser.parse_args()
```

### Run from CLI
```sh
python -m pextractsink \
    <input_sbml> \
    <output_sink> \
    [--compartment_id COMPARTMENT_ID] \
    [--remove_dead_end REMOVE_DEAD_END]
```

## Tests
Test can be run with the following commands:

### Natively
```bash
python -m pytest tests
```

## Authors

* **Melchior du Lac**
* Thomas Duigou, Joan Hérisson

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
