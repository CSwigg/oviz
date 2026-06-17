# Taurus/Sigma Oviz Notebook Handoff

This is the easiest way to share the notebook with someone else: push this repo to GitHub, have them clone it, install `oviz`, put the data files in one folder, then run the notebook.

## Clone And Install

```bash
git clone https://github.com/CSwigg/oviz.git
cd oviz
python -m pip install -e .
python -m pip install jupyterlab scipy matplotlib ipython
```

If they want a clean environment:

```bash
conda create -n oviz-taurus python=3.11
conda activate oviz-taurus
python -m pip install -e .
python -m pip install jupyterlab scipy matplotlib ipython
```

Then open the notebook:

```bash
jupyter lab examples/taurus_core_sigma/taurus_core_sigma_figure_tutorial.ipynb
```

## Data Folder

Make one folder with these files:

```text
data/
  taurus_core_sigma_age-Feb-2025-v2.csv
  cluster_velocities_jan2026.csv
  hunt_sample_chronos_ages_multiprocessing_feb_2026.csv
  mean_and_std_xyz.fits
```

The three CSV files are small enough to share directly through Google Drive, Dropbox, Box, or similar. Do not put the Edenhofer FITS cube in Git.

The Edenhofer Cartesian dust map is large, about 15.7 GB. They should download it themselves:

- Zenodo record: https://zenodo.org/records/8212070
- Direct FITS download: https://zenodo.org/records/8212070/files/mean_and_std_xyz.fits?download=1

## Set Paths

In the notebook, set:

```python
DATA_DIR = Path("/path/to/data")
```

That is enough if all five data files are in the same folder.

They can also set paths from the shell before opening Jupyter:

```bash
export OVIZ_DATA_DIR=/path/to/data
export OVIZ_EDENHOFER_FITS=/path/to/data/mean_and_std_xyz.fits
```

Use the individual `OVIZ_*` variables in the notebook if any file lives somewhere else.

## Output

By default, the notebook writes:

```text
examples/taurus_core_sigma/taurus_core_sigma_figure_tutorial.html
```

That HTML file is generated output and does not need to be committed unless you want to share a pre-rendered figure.
