# dash_gaia_orbits

## Description

This project is a work in progress. It is meant to provide a tool to visualize
the 3D backwards motions of young stellar clusters backwards in time using oribtal integration methods. This type of visual analysis can help reveal the 3D motions of star-forming regions, uncover sites of previous star-forming complexes which have since expanded or dispersed, and it can potentially reveal information about the Galaxy's spiral arms. 

The `StarClusterData` class stores user input data of a star cluster and provides various methods for data analysis.

The `StarClusterCollection` class stores a collection of `StarClusterData` instances.

The `StarClusters3DPlotter` class is used for generating 3D plots of star clusters. It provides a method `generate_3d_plot(collection)` for generating a 3D plot for a `StarClusterCollection`.

The `dash_test.py` file contains functions for creating a Dash application for interactive data visualization. The `layout_2d_plot` function creates the layout for a 2D plot, and the `create_dash_app` function creates the Dash application.

## Installation

To install this project, follow these steps:

1. Clone the repository to your local machine.
2. Navigate into the repo folder and run `python -m pip install .`
3. If still necessary, install the required dependencies, which are listed below, using `pip install -r requirements.txt`

## Dependencies
    - numpy >= 1.26.2
    - pandas >= 1.3.5
    - dash >= 2.14.2
    - plotly >= 5.13.0
    - astropy >= 5.2.1
    - galpy >= 1.9.1

## Usage

To use the `StarClusterData`, `StarClusters3DPlotter`, and Dash application, you need to provide a pandas DataFrame with the following columns: x, y, z, U, V, W, name, age_myr. Here's an example:

```python
import pandas as pd
from cluster_traces import StarClusterData
from StarClusters3DPlotter import StarClusters3DPlotter
from dash_test import create_dash_app




# Create a DataFrame with your data
df = pd.DataFrame({
    'x': [...],
    'y': [...],
    'z': [...],
    'U': [...],
    'V': [...],
    'W': [...],
    'name': [...],
    'age_myr': [...],
})

# Specify the times for the orbits to be integrated
time_array = np.arange(0, -60.5, -0.5)


# Create a StarClusterData instance
cluster_data = StarClusterData(df, 'My Star Cluster')
collection = StarClusterCollection([cluster_data])

# Create a StarClusters3DPlotter instance
plotter = StarClusters3DPlotter()

# Create a Dash application
app = create_dash_app(cluster_data, plotter)
