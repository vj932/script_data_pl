# script_data_pl
This is the code and data accompanying the Model of plasmidless generation by asymmetrical cell division detailed in the Materials and Methods section of "Efficient Method for Estimating Plasmidless Cell Generation Frequency", by Christina Egami, Viraj Arora, and Manel Camps.

Command-line usage for `get_loss.py`, which computes plasmid loss rates (per division) for Carbenicillin and/or Kanamycin treatments.

## Setup

The following files must be in the same folder as `get_loss.py` as the script is written:

- `OD_values_growth_curve_uv_vs_ck_070126.txt` — OD growth data, script currently pulls growht curve data for the delta 6 plasmid. 
- `NORMALIZED_CFU_vs_time.txt` — viability data of WT bacteria under exposure of antibiotic, for both carbenicillin and kanamycin 

Dependencies: `numpy`, `pandas`, `scipy`, `matplotlib`.

## Basic form

```
python get_loss.py <prop_dark> [-g] [-o OUTPUT] [--graph-out GRAPH_OUT]
```

`<prop_dark>` can be one of:

- a **decimal** like `0.12`
- a **ratio** like `3/25` (spaces around `/` are fine)
- a **csv filename** (in this folder) with replicate colony-count data, e.g. `new_data.csv`

Append `k` (case-insensitive) to a numeric input to mark it as **Kanamycin** data. Without a trailing `k`, the input is assumed to be **Carbenicillin**. Whitespace inside the argument is ignored. A csv file input covers both treatments.

By default, numerical results are written to `results.txt` in this folder. Adding `-g` shows the bar chart and also saves it to `loss_rates.png` in this folder.

## Examples of usage

Vavigate to the project folder (cd [path to folder])

Csv file input (both treatments), no graph:
```
python3 get_loss.py new_data.csv
```

Csv file input, with graph:
```
```

```
python3 get_loss.py new_data.csv -g
```

Single decimal proportion, Carbenicillin (default):
```
python3 get_loss.py 0.12
```

Single decimal proportion, Kanamycin:
```
python3 get_loss.py 0.12k
```

Ratio form (dark / total), Carbenicillin:
```
python3 get_loss.py 3/25
```

Ratio form, Kanamycin:
```
python3 get_loss.py "3 / 25 k"
```

Write results to a custom text file:
```
python3 get_loss.py new_data.csv -o results.txt
```

Save the graph to a custom image path:
```
python3 get_loss.py new_data.csv -g --graph-out plot.png
```

Custom text output and custom graph output together:
```
python3 get_loss.py 0.08k -g -o kan_results.txt --graph-out kan_plot.png
```

## Options

| Flag | Description | Default |
|---|---|---|
| `prop_dark` | First argument slot takes decimal(s)/ratio(s) or csv file (optional trailing `k` for Kanamycin when inputting ratios, otherwise assumes ratio is from carb exposure data).
| `-g`, `--graph` | Display the bar plot and save it as an image. | off |
| `-o`, `--output` | Path for the numerical results text file. | `results.txt` in this folder |
| `--graph-out` | Path for the saved graph image (only used with `-g`). | `loss_rates.png` in this folder |

