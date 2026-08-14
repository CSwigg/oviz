"""July 29 source-scene script (plain Python, no notebook).

Running this module builds the July 29 source ThreeJSFigure as
``fig3d`` without writing any file. It was materialized from the
retired notebook pipeline and is now the authoritative source —
edit it directly.
"""

SCENE_SCRIPT_CLUSTER_VELOCITIES = '/Users/swiggumc/Desktop/astro_research/cfa/velocity analysis/outputs/release/cluster_velocities_sdssv_covariance_audited.csv'

import numpy as np 
import pandas as pd
import sys
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import plotly.graph_objects as go


from oviz import Trace, TraceCollection, Animate3D
#from oviz.app import run_dash_app_in_notebook

import matplotlib.pyplot as plt
greys_cmap = plt.get_cmap("Greys")

column_renaming_dict = {'n_stars_hunt' : 'n_stars', 'U_new' : 'U', 'V_new' : 'V', 'W_new' : 'W', 'x_hunt_50' : 'x', 'y_hunt_50' : 'y', 'z_hunt_50' : 'z', 'U_err_new' : 'U_err', 'V_err_new' : 'V_err', 'W_err_new' : 'W_err'}

time_int = np.round(np.arange(0, -121, -1), 1)
static_traces = []
static_traces_times = []

# # read in Gordian scatter map
# dfv = pd.read_csv('/Users/swiggumc/Desktop/astro_research/radcliffe/mustache_work/figures/cubes/edenhofer_scatterized_4pc.csv')
# dfv = dfv.loc[dfv['e'].between(0.003, 1)]
# #dfv = dfv.loc[dfv['extinction'].between(0.0007, 1)]
# ds_index = 1
# scatter_edenhofer = go.Scatter3d(
#     x=dfv['x'].values[::ds_index],
#     y=dfv['y'].values[::ds_index],
#     z=dfv['z'].values[::ds_index],
#     mode='markers',
#     marker=dict(size=2,
#                 color='gray',
#                 symbol='circle',
#                 opacity=.05),
#     line = dict(color = 'gray', width = 0.),
#     name='Edenhofer Dust',
#     visible = True,
#     hovertext='Edenhofer Dust',
#     hoverinfo='skip'
#     )

# df_verg = pd.read_csv('/Users/swiggumc/Downloads/vergely_dust_scatter_large.csv')
# scatter_verg = go.Scatter3d(
#     x=df_verg['X'].values[::ds_index],
#     y=df_verg['Y'].values[::ds_index],
#     z=df_verg['Z'].values[::ds_index],
#     mode='markers',
#     marker=dict(size=2,
#                 color='gray',
#                 symbol='circle',
#                 opacity=.2),
#     line = dict(color = 'gray', width = 0.),
#     name='Vergely dust',
#     visible = True,
#     hovertext='Edenhofer Dust',
#     hoverinfo='skip'
#     )

# # Created in: /Users/swiggumc/Desktop/astro_research/radcliffe/mcallum_maps_scatterized.ipynb
# ds_ha_index = 1
# dfne = pd.read_csv('/Users/swiggumc/Downloads/ne_grid_scatterized.csv')
# scatter_mccallum = go.Scatter3d(
#     x=dfne['x'].values[::ds_ha_index],
#     y=dfne['y'].values[::ds_ha_index],
#     z=dfne['z'].values[::ds_ha_index],
#     mode='markers',
#     marker=dict(size=2,
#                 #color='red',
#                 color='#00BFFF',
#                 symbol='circle',
#                 opacity=.05),
#     line = dict(color = 'red', width = 0.),
#     name='McCallum NE',
#     visible = True,
#     hovertext='McCallum NE',
#     hoverinfo='skip'
#     )

# # append edenhofer
# static_traces.append(scatter_edenhofer)
# static_traces_times.append([0.]) # only show at t=0

# # append vergely
# static_traces.append(scatter_verg)
# static_traces_times.append([0]) # only show at t=0

# # append mccallum
# static_traces.append(scatter_mccallum)
# static_traces_times.append([0]) # only show at t=0

df_hunt_full = pd.read_csv('/Users/swiggumc/Desktop/astro_research/supernovae_map/outputs/paper_variants/chronos_parsec_ybc_flatav_posterior_draws_20260605/parsec/parsec_chronos_mf_mass_velocity_catalog_nomwm_bulkfit.csv')
df_hunt_full['hunt_age_myr'] = pd.to_numeric(
    df_hunt_full['age_myr'] if 'age_myr' in df_hunt_full.columns else pd.Series(np.nan, index=df_hunt_full.index),
    errors='coerce',
)
df_hunt_full['hunt_initial_mass_msun'] = pd.to_numeric(
    df_hunt_full['mass_all_previous'] if 'mass_all_previous' in df_hunt_full.columns else (
        df_hunt_full['mass_all'] if 'mass_all' in df_hunt_full.columns else pd.Series(np.nan, index=df_hunt_full.index)
    ),
    errors='coerce',
)
df_hunt_full['chronos_initial_mass_msun'] = pd.to_numeric(
    df_hunt_full['mass_all'] if 'mass_all' in df_hunt_full.columns else pd.Series(np.nan, index=df_hunt_full.index),
    errors='coerce',
)
df_hunt_full['display_initial_mass_msun'] = df_hunt_full['chronos_initial_mass_msun'].combine_first(
    df_hunt_full['hunt_initial_mass_msun']
)
df_hunt_full['mass_all'] = df_hunt_full['display_initial_mass_msun']
df_hunt_full['mass_source'] = np.where(
    df_hunt_full['chronos_initial_mass_msun'].notnull(),
    'Chronos',
    np.where(df_hunt_full['hunt_initial_mass_msun'].notnull(), 'Hunt', pd.NA),
)
jun6_velocity_cols = ['name', *['x_2026', 'y_2026', 'z_2026', 'x_err_2026', 'y_err_2026', 'z_err_2026', 'U_2026', 'V_2026', 'W_2026', 'U_err_2026', 'V_err_2026', 'W_err_2026', 'n_rvs_2026', 'velocity_fit_status', 'velocity_fit_quality', 'velocity_fit_method']]
jun6_velocity_source = pd.read_csv(
    '/Users/swiggumc/Desktop/astro_research/cfa/velocity analysis/outputs/release/cluster_velocities_sdssv_covariance_audited.csv',
    usecols=lambda _col: _col in jun6_velocity_cols,
)
if jun6_velocity_source['name'].duplicated().any():
    raise RuntimeError('Jun 6 velocity source has duplicate cluster names.')
df_hunt_full = df_hunt_full.drop(
    columns=[_col for _col in jun6_velocity_cols if _col != 'name' and _col in df_hunt_full.columns],
    errors='ignore',
)
df_hunt_full = df_hunt_full.merge(
    jun6_velocity_source,
    on='name',
    how='left',
    validate='1:1',
)
ages_chronos = pd.read_csv('/Users/swiggumc/Desktop/astro_research/supernovae_map/outputs/paper_variants/chronos_parsec_ybc_flatav_posterior_draws_20260605/parsec/parsec_chronos_mode_ages.csv')
if ages_chronos['name'].duplicated().any():
    raise RuntimeError('Jun 6 Chronos mode age source has duplicate cluster names.')
df_hunt_full = df_hunt_full.drop(
    columns=['age_chronos_lo', 'age_chronos_mode', 'age_chronos_hi'],
    errors='ignore',
)
df_hunt_full = df_hunt_full.merge(
    ages_chronos[['name', 'age_chronos_lo', 'age_chronos_mode', 'age_chronos_hi']],
    on='name',
    how='left',
    validate='1:1',
)
for _jun6_col in [
    'age_chronos_lo',
    'age_chronos_mode',
    'age_chronos_hi',
    'x_2026',
    'y_2026',
    'z_2026',
    'U_2026',
    'V_2026',
    'W_2026',
    'U_err_2026',
    'V_err_2026',
    'W_err_2026',
    'n_rvs_2026',
]:
    if _jun6_col in df_hunt_full.columns:
        df_hunt_full[_jun6_col] = pd.to_numeric(df_hunt_full[_jun6_col], errors='coerce')
df_hunt_full['age_myr'] = df_hunt_full['age_chronos_mode'].combine_first(df_hunt_full['hunt_age_myr'])
df_hunt_full['age_source'] = np.where(
    df_hunt_full['age_chronos_mode'].notnull(),
    'Chronos',
    np.where(df_hunt_full['hunt_age_myr'].notnull(), 'Hunt', pd.NA),
)
df_hunt_full = df_hunt_full.drop(
    columns=['x', 'y', 'z', 'U', 'V', 'W', 'U_err', 'V_err', 'W_err'],
    errors='ignore',
)
df_hunt_full = df_hunt_full.rename(columns={
    'U_2026': 'U',
    'V_2026': 'V',
    'W_2026': 'W',
    'U_err_2026': 'U_err',
    'V_err_2026': 'V_err',
    'W_err_2026': 'W_err',
    'x_2026': 'x',
    'y_2026': 'y',
    'z_2026': 'z',
})
print(
    'Using cluster velocities from /Users/swiggumc/Desktop/astro_research/cfa/velocity analysis/outputs/release/cluster_velocities_sdssv_covariance_audited.csv: '
    f"{df_hunt_full['age_myr'].notnull().sum()} clusters with Chronos mode ages, "
    f"{df_hunt_full[['x', 'y', 'z', 'U', 'V', 'W']].notnull().all(axis=1).sum()} with complete positions/velocities."
)
df_hunt_good = df_hunt_full.loc[
    (df_hunt_full['U_err'] < 10) &
    (df_hunt_full['V_err'] < 10) &
    (df_hunt_full['W_err'] < 10) &
    (df_hunt_full['U'].notnull()) &
    (df_hunt_full['V'].notnull()) &
    (df_hunt_full['W'].notnull()) &
    (df_hunt_full['x'].notnull()) &
    (df_hunt_full['y'].notnull()) &
    (df_hunt_full['z'].notnull()) &
    (df_hunt_full['age_myr'].notnull()) &
    (df_hunt_full['x'].between(-4000.0, 4000.0)) &
    (df_hunt_full['y'].between(-4000.0, 4000.0)) &
    (df_hunt_full['n_rvs_2026'] >= 3)
]
df_hunt_60 = df_hunt_good.loc[df_hunt_good['age_myr'] < 60]
df_hunt_young = df_hunt_good.loc[df_hunt_good['age_myr'] < 15]
df_hunt_mid = df_hunt_good.loc[df_hunt_good['age_myr'].between(15, 30)]
df_hunt_old = df_hunt_good.loc[df_hunt_good['age_myr'].between(30, 60)]


chronos_cluster_columns = set(pd.read_csv(
    '/Users/swiggumc/Desktop/astro_research/chronos_fasrc/runs/current/chronos/parsec_allhunt_46w_500b_5000s_dustav_12gyr_linearage_192shards/cluster_results.csv',
    nrows=0,
).columns)
chronos_model_key = 'parsec'
chronos_model_display = 'PARSEC'
chronos_age_value_col = 'parsec_age_mode'
chronos_required_cols = [
    'name',
    'parsec_age_lo',
    chronos_age_value_col,
    'parsec_age_hi',
    'parsec_mass_cluster_imf_corrected',
    'parsec_status',
]
chronos_missing_cols = [col for col in chronos_required_cols if col not in chronos_cluster_columns]
if chronos_missing_cols:
    raise RuntimeError(
        f"Could not find {chronos_missing_cols!r} in '/Users/swiggumc/Desktop/astro_research/chronos_fasrc/runs/current/chronos/parsec_allhunt_46w_500b_5000s_dustav_12gyr_linearage_192shards/cluster_results.csv'."
    )
chronos_cluster_ages = pd.read_csv(
    '/Users/swiggumc/Desktop/astro_research/chronos_fasrc/runs/current/chronos/parsec_allhunt_46w_500b_5000s_dustav_12gyr_linearage_192shards/cluster_results.csv',
    usecols=chronos_required_cols,
)
chronos_cluster_ages = chronos_cluster_ages.rename(
    columns={
        'parsec_age_lo': 'chronos_age_lo_myr',
        chronos_age_value_col: 'chronos_age_myr',
        'parsec_age_hi': 'chronos_age_hi_myr',
        'parsec_mass_cluster_imf_corrected': 'chronos_initial_mass_msun',
        'parsec_status': 'chronos_status',
    }
)
for _chronos_col in [
    'chronos_age_lo_myr',
    'chronos_age_myr',
    'chronos_age_hi_myr',
    'chronos_initial_mass_msun',
]:
    chronos_cluster_ages[_chronos_col] = pd.to_numeric(
        chronos_cluster_ages[_chronos_col],
        errors='coerce',
    )
chronos_cluster_ages = chronos_cluster_ages.loc[
    chronos_cluster_ages['chronos_status'].eq('success')
    & chronos_cluster_ages['chronos_age_myr'].notnull()
].copy()
df_hunt_chronos_full = df_hunt_full.drop(
    columns=[
        'chronos_age_lo_myr',
        'chronos_age_myr',
        'chronos_age_hi_myr',
        'chronos_initial_mass_msun',
        'chronos_status',
    ],
    errors='ignore',
).merge(
    chronos_cluster_ages,
    on='name',
    how='left',
    validate='m:1',
)
if 'hunt_age_myr' not in df_hunt_chronos_full.columns:
    df_hunt_chronos_full['hunt_age_myr'] = pd.to_numeric(
        df_hunt_chronos_full['age_myr'] if 'age_myr' in df_hunt_chronos_full.columns else pd.Series(np.nan, index=df_hunt_chronos_full.index),
        errors='coerce',
    )
if 'hunt_initial_mass_msun' not in df_hunt_chronos_full.columns:
    df_hunt_chronos_full['hunt_initial_mass_msun'] = pd.to_numeric(
        df_hunt_chronos_full['mass_all_previous'] if 'mass_all_previous' in df_hunt_chronos_full.columns else (
            df_hunt_chronos_full['mass_all'] if 'mass_all' in df_hunt_chronos_full.columns else pd.Series(np.nan, index=df_hunt_chronos_full.index)
        ),
        errors='coerce',
    )
df_hunt_chronos_full['display_age_myr'] = df_hunt_chronos_full['chronos_age_myr'].combine_first(
    df_hunt_chronos_full['hunt_age_myr']
)
df_hunt_chronos_full['chronos_initial_mass_msun'] = df_hunt_chronos_full['chronos_initial_mass_msun'].combine_first(
    df_hunt_chronos_full['hunt_initial_mass_msun']
)
df_hunt_chronos_full['display_initial_mass_msun'] = df_hunt_chronos_full['chronos_initial_mass_msun']
df_hunt_chronos_full['mass_all'] = df_hunt_chronos_full['display_initial_mass_msun']
main_figure_cluster_xy_min_pc = -4000.0
main_figure_cluster_xy_max_pc = 4000.0
main_figure_cluster_young_age_myr = 15
main_figure_cluster_blue_max_age_myr = 120
main_figure_cluster_grey_max_age_myr = 150
df_hunt_chronos_sample = df_hunt_chronos_full.loc[
    (df_hunt_chronos_full['U_err'] < 10) &
    (df_hunt_chronos_full['V_err'] < 10) &
    (df_hunt_chronos_full['W_err'] < 10) &
    (df_hunt_chronos_full['U'].notnull()) &
    (df_hunt_chronos_full['V'].notnull()) &
    (df_hunt_chronos_full['W'].notnull()) &
    (df_hunt_chronos_full['x'].notnull()) &
    (df_hunt_chronos_full['y'].notnull()) &
    (df_hunt_chronos_full['z'].notnull()) &
    (df_hunt_chronos_full['x'].between(main_figure_cluster_xy_min_pc, main_figure_cluster_xy_max_pc)) &
    (df_hunt_chronos_full['y'].between(main_figure_cluster_xy_min_pc, main_figure_cluster_xy_max_pc)) &
    (df_hunt_chronos_full['display_age_myr'].notnull()) &
    (df_hunt_chronos_full['n_rvs_2026'] >= 3)
].copy()
df_hunt_chronos_sample['age_myr'] = df_hunt_chronos_sample['display_age_myr']
df_hunt_chronos_sample['age_source'] = np.where(
    df_hunt_chronos_sample['chronos_age_myr'].notnull(),
    chronos_model_display,
    'Hunt',
)
df_hunt_60_chronos = df_hunt_chronos_sample.loc[
    df_hunt_chronos_sample['age_myr'] < main_figure_cluster_blue_max_age_myr
].copy()
df_hunt_0_to_150_chronos = df_hunt_chronos_sample.loc[
    df_hunt_chronos_sample['age_myr'] < main_figure_cluster_grey_max_age_myr
].copy()
df_hunt_young_chronos = df_hunt_chronos_sample.loc[
    df_hunt_chronos_sample['age_myr'] < main_figure_cluster_young_age_myr
].copy()
print(
    f"Using {chronos_model_display} ages from /Users/swiggumc/Desktop/astro_research/chronos_fasrc/runs/current/chronos/parsec_allhunt_46w_500b_5000s_dustav_12gyr_linearage_192shards/cluster_results.csv for displayed cluster samples: "
    f"{len(df_hunt_60_chronos)} clusters <{main_figure_cluster_blue_max_age_myr:g} Myr, "
    f"{len(df_hunt_0_to_150_chronos)} clusters 0-{main_figure_cluster_grey_max_age_myr:g} Myr, "
    f"{len(df_hunt_young_chronos)} clusters <{main_figure_cluster_young_age_myr:g} Myr. "
    f"Displayed ages use {chronos_age_value_col} with Hunt fallbacks when Chronos is missing. "
    f"x/y within [{main_figure_cluster_xy_min_pc:g}, {main_figure_cluster_xy_max_pc:g}] pc; no z cut applied."
)
c_s24 = pd.read_csv('/Users/swiggumc/Downloads/cluster_sample_data.csv')
ap_family = df_hunt_60.loc[df_hunt_60['name'].isin(c_s24.loc[c_s24['family'] == 'alphaPer']['name'])]
cr135_family = df_hunt_60.loc[df_hunt_60['name'].isin(c_s24.loc[c_s24['family'] == 'cr135']['name'])]
m6_family = df_hunt_60.loc[df_hunt_60['name'].isin(c_s24.loc[c_s24['family'] == 'm6']['name'])]

# Lacerta family
lacerta_family_names = ['CWNU_1243', 'HSC_661', 'UPK_109', 'HSC_705', 'CWNU_96', 'Theia_420', 'UPK_166', 'UPK_168', 'Teutsch_39', 'OCSN_32', 'Theia_100', 'CWNU_311']
trumpler3_family_names = ['Trumpler_3', 'Theia_850', 'CWNU_525', 'FSR_0686', 'Theia_1722', 'UPK_325', 'UPK_307', 'FSR_0732', 'NGC_1960', 'HSC_1350', 'HSC_1341', 'HXHWL_18', 'HXWHB_8', 'NGC_1502', 'CWNU_409', 'CWNU_364', 'NGC_1444', 'UBC_51', 'Berkeley_14A', 'CWNU_518', 'CWNU_205', 'HSC_1308']
proto_orion_family_names = ['HSC_1340', 'NGC_2232', 'CWNU_1057', 'CWNU_1111', 'CWNU_1052', 'Theia_71', 'ZHBJZ_1']

lacerta_family = df_hunt_good.loc[df_hunt_good['name'].isin(lacerta_family_names)]
trumpler3_family = df_hunt_good.loc[df_hunt_good['name'].isin(trumpler3_family_names)]
proto_orion_family = df_hunt_good.loc[df_hunt_good['name'].isin(proto_orion_family_names)]

rw_region = pd.read_csv('/Users/swiggumc/Downloads/rw_clusters.csv')['name'].values
split_region = pd.read_csv('/Users/swiggumc/Downloads/the_split_selection.csv')['name'].values
cepheus_region = pd.read_csv('/Users/swiggumc/Downloads/cepheus_spur_region.csv')['name'].values
vel_sag_region = ['Pismis_5', 'Collinder_197', 'FSR_1421', 'BH_56', 'HSC_2144', 'ASCC_79', 'UPK_604', 'FSR_1694', 'NGC_6178', 'OC_0684', 'HSC_2911', 'CWNU_170', 'NGC_6383']

rw = df_hunt_young.loc[df_hunt_young['name'].isin(rw_region)]
split = df_hunt_young.loc[df_hunt_young['name'].isin(split_region)]
cepheus = df_hunt_young.loc[df_hunt_young['name'].isin(cepheus_region)]
vel_sag = df_hunt_young.loc[df_hunt_young['name'].isin(vel_sag_region)]

moca = pd.read_csv('/Users/swiggumc/Downloads/MOCA DB - Web UI.csv')
moca = moca.rename(columns={'u' : 'U', 'v' : 'V', 'w' : 'W'})

# If you previously did: moca['name'] = moca['name'] + moca['coeff_index']
# then use the original association label instead of this modified name.
group_col = "name"

median_cols = ["age_myr", "x", "y", "z", "U", "V", "W"]

# Ensure numeric types
for c in median_cols + ["nstars"]:
    moca[c] = pd.to_numeric(moca[c], errors="coerce")

moca_med = (
    moca.groupby(group_col, as_index=False)
        .agg(
            moca_aid=("moca_aid", "first"),
            age_ref=("age_ref", "first"),
            nstars=("nstars", "sum"),          # use "first" if nstars is per-association already
            n_coeff=("coeff_index", "nunique"),
            **{c: (c, "median") for c in median_cols},
        )
)
moca_med = moca_med.loc[
    (moca_med['age_myr'].notnull()) &
    (moca_med['age_myr'] < 60)
    ]  # filter to associations with at least 3 stars

sun = pd.DataFrame({'name' : ['Sun'], 'family' : ['Sun'], 'age_myr' : [8000], 'U' : [0.], 'V' : [0.], 'W' : [0.], 'x' : [0.], 'y' : [0.], 'z' : [27.], 'n_stars' : [1]})

young_trace = Trace(df_hunt_young_chronos, data_name = 'Clusters (< 15 Myr)', min_size = 1.25, max_size = 7, color = 'red', opacity = 1, marker_style = 'circle', show_tracks = False, size_by_n_stars=False)
full_sample_trace = Trace(df_hunt_60_chronos, data_name = 'Clusters (< 60 Myr)', min_size = 0, max_size = 7, color = '#2f80ff', opacity = .7, marker_style = 'circle', show_tracks = False, size_by_n_stars=False)

ap_trace = Trace(ap_family, data_name = 'Alpha Persei Family', min_size = 0., max_size = 7, color = 'violet', opacity = 1, marker_style = 'circle', show_tracks = False, size_by_n_stars=False)
cr135_trace = Trace(cr135_family, data_name = 'Cr 135 Family', min_size = 0., max_size = 7, color = 'orange', opacity = 1, marker_style = 'circle', show_tracks = False, size_by_n_stars=False)
m6_trace = Trace(m6_family, data_name = 'M6 Family', min_size = 0., max_size = 7, color = 'cyan', opacity = 1, marker_style = 'circle', show_tracks = False, size_by_n_stars=False)

lacerta_trace = Trace(lacerta_family, data_name = 'Lacerta Family', min_size = 0., max_size = 7, color = '#F0E68C', opacity = 1, marker_style = 'circle', show_tracks = False, size_by_n_stars=False)
trumpler3_trace = Trace(trumpler3_family, data_name = 'Trumpler 3 Family', min_size = 0., max_size = 7, color = '#00CED1', opacity = 1, marker_style = 'circle', show_tracks = False, size_by_n_stars=False)
proto_orion_trace = Trace(proto_orion_family, data_name = 'Proto Orion Family', min_size = 0., max_size = 7, color = '#DCDCDC', opacity = 1, marker_style = 'circle', show_tracks = False, size_by_n_stars=False)

rw_trace = Trace(rw, data_name = 'Radcliffe Wave Clusters', min_size = 1.25, max_size = 7, color = 'red', opacity = 1, marker_style = 'circle', show_tracks = True, size_by_n_stars=False)
split_trace = Trace(split, data_name = 'Split Selection Clusters', min_size = 1.25, max_size = 7, color = 'red', opacity = 1, marker_style = 'circle', show_tracks = True, size_by_n_stars=False)
cepheus_trace = Trace(cepheus, data_name = 'Cepheus Spur Clusters', min_size = 1.25, max_size = 7, color = 'red', opacity = 1, marker_style = 'circle', show_tracks = True, size_by_n_stars=False)
vel_sag_trace = Trace(vel_sag, data_name = 'Vela-Sagittarius Clusters', min_size = 1.25, max_size = 7, color = 'red', opacity = 1, marker_style = 'circle', show_tracks = True, size_by_n_stars=False)

#moca_trace = Trace(moca_med, data_name = 'MOCA Sample', min_size = 0, max_size = 12, color = 'blue', opacity = 1, marker_style = 'circle', show_tracks = False, size_by_n_stars=False)

sun_trace = Trace(sun, data_name = 'Sun', min_size = 5, max_size = 5, color = 'yellow', opacity = 1, marker_style = 'circle', show_tracks = False, size_by_n_stars=False)



ratzenboeck_scocen = pd.read_csv('/Users/swiggumc/Desktop/astro_research/radcliffe/cluster_data/files_formatted/ratzi_2022.csv')
from astropy.io import fits as _ratzenboeck_fits
with _ratzenboeck_fits.open('/Users/swiggumc/Downloads/ScoCen_SigMA2_groups_averages_34+KERR3_CenFar-to-Cham_37clsuters.fits', memmap=True) as _ratzenboeck_hdus:
    _ratzenboeck_group_data = _ratzenboeck_hdus[1].data
    ratzenboeck_scocen_metadata = pd.DataFrame({
        'sigma_id': np.asarray(_ratzenboeck_group_data['sigma_label'], dtype=int),
        'sigma_region': [
            str(_value).strip()
            for _value in _ratzenboeck_group_data['scocen_region']
        ],
        'ratzenboeck_group': [
            str(_value).strip()
            for _value in _ratzenboeck_group_data['cluster_name']
        ],
        'n_sigma_sources': np.asarray(
            _ratzenboeck_group_data['nr_sources'],
            dtype=int,
        ),
    })
ratzenboeck_required_cols = ['name', 'age_myr', 'x', 'y', 'z', 'U', 'V', 'W']
ratzenboeck_missing_cols = [
    _col for _col in ratzenboeck_required_cols
    if _col not in ratzenboeck_scocen.columns
]
if ratzenboeck_missing_cols:
    raise RuntimeError(
        f'Ratzenboeck/SigMA catalog is missing {ratzenboeck_missing_cols!r}.'
    )
ratzenboeck_metadata_required_cols = [
    'sigma_id',
    'sigma_region',
    'ratzenboeck_group',
    'n_sigma_sources',
]
ratzenboeck_metadata_missing_cols = [
    _col for _col in ratzenboeck_metadata_required_cols
    if _col not in ratzenboeck_scocen_metadata.columns
]
if ratzenboeck_metadata_missing_cols:
    raise RuntimeError(
        f'Ratzenboeck/SigMA metadata is missing {ratzenboeck_metadata_missing_cols!r}.'
    )
if ratzenboeck_scocen['name'].duplicated().any():
    raise RuntimeError('Ratzenboeck/SigMA catalog has duplicate cluster names.')
if ratzenboeck_scocen_metadata['ratzenboeck_group'].duplicated().any():
    raise RuntimeError('Ratzenboeck/SigMA metadata has duplicate cluster names.')
ratzenboeck_scocen = ratzenboeck_scocen.merge(
    ratzenboeck_scocen_metadata[ratzenboeck_metadata_required_cols],
    left_on='name',
    right_on='ratzenboeck_group',
    how='left',
    validate='1:1',
)
if len(ratzenboeck_scocen) != 37:
    raise RuntimeError(
        f'Expected 37 Ratzenboeck/SigMA clusters, found {len(ratzenboeck_scocen)}.'
    )
if ratzenboeck_scocen[ratzenboeck_metadata_required_cols].isnull().any().any():
    raise RuntimeError('Could not match every Ratzenboeck/SigMA cluster to metadata.')
ratzenboeck_scocen['n_stars'] = pd.to_numeric(
    ratzenboeck_scocen['n_sigma_sources'],
    errors='raise',
)
ratzenboeck_scocen['name_all'] = (
    'SigMA ' + ratzenboeck_scocen['sigma_id'].astype(int).astype(str)
    + '; ' + ratzenboeck_scocen['sigma_region'].astype(str)
)
ratzenboeck_scocen_trace = Trace(
    ratzenboeck_scocen,
    data_name='Sco-Cen (Ratzenböck/SigMA)',
    min_size=2.5,
    max_size=10.0,
    color='#c77dff',
    opacity=0.95,
    marker_style='circle',
    show_tracks=True,
    size_by_n_stars=True,
)
print(
    f'Loaded {len(ratzenboeck_scocen)} Ratzenboeck/SigMA Sco-Cen clusters '
    f'from /Users/swiggumc/Desktop/astro_research/radcliffe/cluster_data/files_formatted/ratzi_2022.csv.'
)
full_catalog_trace = Trace(
    df_hunt_good,
    data_name='Full Cluster Catalog',
    min_size=0.0,
    max_size=4.0,
    color='#8f959d',
    opacity=0.16,
    marker_style='circle',
    show_tracks=False,
    size_by_n_stars=False,
)
from astropy.io import fits
hdu = fits.open('/Users/swiggumc/Downloads/ONeill2024_LocalBubble_Shell_xyz.fits')
hdu.info()

from astropy.io import fits
hdu = fits.open('/Users/swiggumc/Downloads/mean_and_std_xyz-2.fits')
hdu.info()

traces = TraceCollection([
    sun_trace,

    full_catalog_trace,
    ratzenboeck_scocen_trace,
    full_sample_trace,
    young_trace, 
    ap_trace, 
    cr135_trace, 
    m6_trace, 
    lacerta_trace,
    trumpler3_trace,
    proto_orion_trace,
    rw_trace, 
    split_trace, 
    cepheus_trace, 
    vel_sag_trace
    #moca_trace
    ])

trace_groupings = {
    "Clusters": ['Sun', 'Clusters (< 60 Myr)', 'Clusters (< 15 Myr)'],
    "Cluster Families": ['Sun', 'Alpha Persei Family', 'Cr 135 Family', 'M6 Family', 'Lacerta Family', 'Trumpler 3 Family', 'Proto Orion Family'],
    "Dust Structures": ['Sun', 'Radcliffe Wave Clusters', 'Split Selection Clusters', 'Cepheus Spur Clusters', 'Vela-Sagittarius Clusters'],
    "Dust Structures and Families": ['Sun', 'Alpha Persei Family', 'Cr 135 Family', 'M6 Family', 'Lacerta Family', 'Trumpler 3 Family', 'Proto Orion Family', 'Radcliffe Wave Clusters', 'Split Selection Clusters', 'Cepheus Spur Clusters', 'Vela-Sagittarius Clusters']
}


### NEW THREEJS
xyz_widths = (2000, 2000, 400)
_all_group_names = ['Sun', 'Full Cluster Catalog']
for _group_names in trace_groupings.values():
    for _name in _group_names:
        if _name not in _all_group_names and _name != 'Full Cluster Catalog':
            _all_group_names.append(_name)
trace_groupings = {
    'All': _all_group_names,
    **{
        key: ['Sun'] + [
            name for name in value
            if name not in {'Sun', 'Full Cluster Catalog'}
        ]
        for key, value in trace_groupings.items()
    },
}

ratzenboeck_scocen_trace_name = 'Sco-Cen (Ratzenböck/SigMA)'
for _ratzenboeck_group_name in ['All', 'Clusters']:
    if ratzenboeck_scocen_trace_name not in trace_groupings.get(
        _ratzenboeck_group_name,
        [],
    ):
        trace_groupings.setdefault(_ratzenboeck_group_name, ['Sun']).append(
            ratzenboeck_scocen_trace_name
        )
trace_groupings[ratzenboeck_scocen_trace_name] = [
    'Sun',
    ratzenboeck_scocen_trace_name,
]

plot_3d = Animate3D(
    data_collection = traces, 
    xyz_widths = xyz_widths, 
    figure_theme = 'dark', 
    trace_grouping_dict=trace_groupings
    )
#save_name = '/Users/swiggumc/Desktop/astro_research/cam_website/interactive_figures/main_figure_hi_res_dust.html'
save_name = None

edenhofer_volume = {
    "name": "Edenhofer+2024 Dust",
    "path": "/Users/swiggumc/Downloads/mean_and_std_xyz-2.fits",
    "hdu": "MEAN",
    "clip_bounds": {"z": [-400.0, 400.0]},
    #"max_resolution": 512,
    "max_resolution": 512,
    "max_resolution_cap": 512,
    "opacity": 1.0,
    "samples": 200,
    "alpha_coef": 105.0,
    "vmin": 0,
    "vmax": .0385,
    "colormap": greys_cmap,   # or just "ice"
    "supports_show_all_times": True,
    "co_rotate_with_frame": True,
}
vergely_volume = {
    "name": "Vergely 3D Dust",
    "path": "/Users/swiggumc/Downloads/vergely_3D_Dust.fits",
    'hdu' : 'Primary',
    "max_resolution": 100,
    "opacity": 1.,
    "samples": 24,
    "alpha_coef": 50.0,
    "colormap": greys_cmap,   # or "inferno"
    # no hdu needed
}
oneilLB_volume = {
    "name": "O'Neill+2024 Local Bubble Shell",
    "path": "/Users/swiggumc/Downloads/ONeill2024_LocalBubble_Shell_xyz.fits",
    'hdu' : 'Shell',
    "max_resolution": 100,
    "opacity": 1.,
    "samples": 24,
    "alpha_coef": 200,
    "colormap": 'magma',   # or "inferno"
     # no hdu needed
}
mccallum_ne = {
    "name": "McCallum+ NE",
    "path": "/Users/swiggumc/Downloads/ne_grid.fits",
    'hdu' : 'PRIMARY',
    "max_resolution": 100,
    "opacity": 1.,
    "samples": 24,
    "alpha_coef": 200,
    "colormap": 'magma',   # or "inferno"
     # no hdu needed
}

import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")
os.environ.setdefault("MPLBACKEND", "Agg")

from pathlib import Path
import sys

SUPERNOVAE_MAP_ROOT = Path('/Users/swiggumc/Desktop/astro_research/supernovae_map')
SUPERNOVAE_CATALOG_PATH = Path('/Users/swiggumc/Desktop/astro_research/supernovae_map/paper/solar_encounter_catalog_current.csv.gz')
if str(SUPERNOVAE_MAP_ROOT) not in sys.path:
    sys.path.insert(0, str(SUPERNOVAE_MAP_ROOT))

from mapper import OrbitConfig
from mapper.orbits import ClusterOrbiter


def _sn_build_edges(half_width_pc, voxel_size_pc):
    n_bins = int(np.ceil((2.0 * float(half_width_pc)) / float(voxel_size_pc)))
    n_bins = max(n_bins, 1)
    half_extent = 0.5 * n_bins * float(voxel_size_pc)
    return np.linspace(-half_extent, half_extent, n_bins + 1, dtype=float)


def _sn_nearest_values(sorted_times, query_times):
    idx = np.searchsorted(sorted_times, query_times, side="left")
    idx = np.clip(idx, 1, len(sorted_times) - 1)
    left = sorted_times[idx - 1]
    right = sorted_times[idx]
    choose_left = (query_times - left) <= (right - query_times)
    nearest_idx = np.where(choose_left, idx - 1, idx)
    return sorted_times[nearest_idx]


def _sn_rotating_local_to_fixed_galactocentric(x_rot_pc, y_rot_pc, z_rot_pc, time_myr, *, orbit_cfg):
    r_sun_pc = float(orbit_cfg.ro_kpc) * 1000.0
    omega_rot_per_myr = (float(orbit_cfg.vo_kms) / float(orbit_cfg.ro_kpc)) / 10.0
    time_code_units = np.asarray(time_myr, dtype=float) * 0.01022

    local_x = r_sun_pc - np.asarray(x_rot_pc, dtype=float)
    local_y = np.asarray(y_rot_pc, dtype=float)
    r_pc = np.hypot(local_x, local_y)
    phi = np.arctan2(local_y, local_x)
    theta = phi + omega_rot_per_myr * time_code_units - (0.5 * np.pi)

    x_gc_pc = r_pc * np.sin(theta)
    y_gc_pc = r_pc * np.cos(theta)
    z_gc_pc = np.asarray(z_rot_pc, dtype=float)
    return x_gc_pc, y_gc_pc, z_gc_pc


def _sn_build_reference_orbit_df(time_grid, *, orbit_cfg):
    rf = pd.DataFrame(
        {
            "name": ["rf"],
            "x": [0.0],
            "y": [0.0],
            "z": [0.0],
            "U": [-11.1],
            "V": [-12.24],
            "W": [-7.25],
            "x_err": [0.0],
            "y_err": [0.0],
            "z_err": [0.0],
            "U_err": [0.0],
            "V_err": [0.0],
            "W_err": [0.0],
        }
    )
    rf_df = ClusterOrbiter(rf, config=orbit_cfg).run_mc_integration(time_grid, n_samples=1).copy()
    if rf_df.empty:
        raise ValueError("Reference orbit integration returned no rows.")
    rf_df = rf_df.loc[rf_df["sample_id"] == 0, ["time_myr", "x_gc_pc", "y_gc_pc", "z_gc_pc"]].copy()
    rf_df = rf_df.rename(
        columns={
            "x_gc_pc": "rf_x_gc_pc",
            "y_gc_pc": "rf_y_gc_pc",
            "z_gc_pc": "rf_z_gc_pc",
        }
    )
    if not np.isclose(rf_df["time_myr"].to_numpy(dtype=float), 0.0, atol=1e-9).any():
        x0_gc_pc, y0_gc_pc, z0_gc_pc = _sn_rotating_local_to_fixed_galactocentric(
            np.array([0.0]),
            np.array([0.0]),
            np.array([0.0]),
            np.array([0.0]),
            orbit_cfg=orbit_cfg,
        )
        rf_df = pd.concat(
            [
                rf_df,
                pd.DataFrame(
                    {
                        "time_myr": [0.0],
                        "rf_x_gc_pc": [float(x0_gc_pc[0])],
                        "rf_y_gc_pc": [float(y0_gc_pc[0])],
                        "rf_z_gc_pc": [float(z0_gc_pc[0])],
                    }
                ),
            ],
            ignore_index=True,
        )
    return rf_df.sort_values("time_myr").reset_index(drop=True)


def _sn_assign_supernova_positions_in_trace_frame(sne_df, *, rf_df, orbit_cfg):
    df = sne_df.copy()
    time_grid = np.sort(rf_df["time_myr"].unique())
    if len(time_grid) == 0:
        raise ValueError("Reference orbit lookup is empty.")

    df["trace_time_myr"] = _sn_nearest_values(time_grid, df["time_of_death_myr"].to_numpy(dtype=float))
    joined = df[["trace_time_myr"]].merge(
        rf_df.rename(columns={"time_myr": "trace_time_myr"}),
        on="trace_time_myr",
        how="left",
        sort=False,
    )
    x_gc_pc, y_gc_pc, z_gc_pc = _sn_rotating_local_to_fixed_galactocentric(
        df["x_pc"].to_numpy(dtype=float),
        df["y_pc"].to_numpy(dtype=float),
        df["z_pc"].to_numpy(dtype=float),
        df["trace_time_myr"].to_numpy(dtype=float),
        orbit_cfg=orbit_cfg,
    )
    df["x_trace_frame_pc"] = x_gc_pc - joined["rf_x_gc_pc"].to_numpy(dtype=float)
    df["y_trace_frame_pc"] = y_gc_pc - joined["rf_y_gc_pc"].to_numpy(dtype=float)
    df["z_trace_frame_pc"] = z_gc_pc - joined["rf_z_gc_pc"].to_numpy(dtype=float)
    return df


def _sn_build_supernova_volume_layers(
    sne_df,
    *,
    time_grid,
    x_edges,
    y_edges,
    z_edges,
    gaussian_sigma_vox=0.8,
    time_window_half_width_myr=5.0,
):
    work = sne_df.copy()
    work = work[np.isfinite(work["time_of_death_myr"])]
    work = work[
        np.isfinite(work["x_trace_frame_pc"])
        & np.isfinite(work["y_trace_frame_pc"])
        & np.isfinite(work["z_trace_frame_pc"])
    ]
    work = work.loc[
        (work["time_of_death_myr"] <= float(np.max(time_grid)) + float(time_window_half_width_myr))
        & (work["time_of_death_myr"] >= float(np.min(time_grid)) - float(time_window_half_width_myr))
    ].copy()

    try:
        from scipy.ndimage import gaussian_filter
    except ImportError:
        gaussian_filter = None

    smoothed_by_time = {}
    positive_values = []

    for time_value in time_grid:
        frame_events = work.loc[
            np.abs(work["time_of_death_myr"].to_numpy(dtype=float) - float(time_value))
            <= float(time_window_half_width_myr) + 1e-9
        ]
        if frame_events.empty:
            cube_zyx = np.zeros(
                (len(z_edges) - 1, len(y_edges) - 1, len(x_edges) - 1),
                dtype=np.float32,
            )
        else:
            sample = frame_events[["z_trace_frame_pc", "y_trace_frame_pc", "x_trace_frame_pc"]].to_numpy(dtype=float)
            cube_zyx, _ = np.histogramdd(sample, bins=(z_edges, y_edges, x_edges))
            cube_zyx = cube_zyx.astype(np.float32, copy=False)

        if gaussian_filter is not None and float(gaussian_sigma_vox) > 0:
            cube_zyx = gaussian_filter(cube_zyx, sigma=float(gaussian_sigma_vox), mode="constant").astype(
                np.float32,
                copy=False,
            )

        smoothed_by_time[float(time_value)] = cube_zyx
        positive = cube_zyx[cube_zyx > 0]
        if positive.size:
            positive_values.append(positive)

    if positive_values:
        all_positive = np.concatenate(positive_values)
        data_max = float(np.nanmax(all_positive))
        default_vmin = float(np.nanquantile(all_positive, 0.70))
        default_vmax = float(np.nanquantile(all_positive, 0.995))
    else:
        data_max = 1.0
        default_vmin = 0.05
        default_vmax = 1.0

    if not data_max > 0:
        data_max = 1.0
    if not default_vmax > default_vmin:
        default_vmin = 0.05 * data_max
        default_vmax = data_max

    bounds = {
        "x": [float(x_edges[0]), float(x_edges[-1])],
        "y": [float(y_edges[0]), float(y_edges[-1])],
        "z": [float(z_edges[0]), float(z_edges[-1])],
    }

    volumes = []
    for time_index, time_value in enumerate(time_grid):
        volumes.append(
            {
                "key": f"supernova-density-{time_index:03d}",
                "state_key": "supernova-density",
                "state_name": "Supernova Density",
                "name": f"Supernova Density | {abs(float(time_value)):.0f} Myr",
                "time_myr": float(time_value),
                "data": smoothed_by_time[float(time_value)],
                "bounds": bounds,
                "data_range": [0.0, float(data_max)],
                "vmin": float(default_vmin),
                "vmax": float(default_vmax),
                "opacity": 0.24,
                "alpha_coef": 95.0,
                "gradient_step": 0.01,
                "samples": 240,
                "colormap": "ice",
                "interpolation": True,
                "visible": False,
            }
        )
    return volumes


def _build_supernova_volumes_for_main_figure(time_values):
    if not SUPERNOVAE_CATALOG_PATH.exists():
        print(f"Skipping supernova volumes; missing catalog: {SUPERNOVAE_CATALOG_PATH}")
        return []

    time_grid = np.asarray(time_values, dtype=float)
    orbit_cfg = OrbitConfig(backend="galpy", spiral_model="none")
    sne_df = pd.read_csv(
        SUPERNOVAE_CATALOG_PATH,
        usecols=["time_of_death_myr", "x_pc", "y_pc", "z_pc"],
    )
    rf_df = _sn_build_reference_orbit_df(time_grid, orbit_cfg=orbit_cfg)
    sne_trace_df = _sn_assign_supernova_positions_in_trace_frame(
        sne_df,
        rf_df=rf_df,
        orbit_cfg=orbit_cfg,
    )
    x_edges = _sn_build_edges(2000.0, 50.0)
    y_edges = _sn_build_edges(2000.0, 50.0)
    z_edges = _sn_build_edges(400.0, 50.0)
    return _sn_build_supernova_volume_layers(
        sne_trace_df,
        time_grid=time_grid,
        x_edges=x_edges,
        y_edges=y_edges,
        z_edges=z_edges,
        gaussian_sigma_vox=0.8,
        time_window_half_width_myr=5.0,
    )


supernova_volumes = _build_supernova_volumes_for_main_figure(time_int)


import os
from pathlib import Path as _OvizPath

MCCALLUM_NE_GRID_PATH = _OvizPath(
    os.environ.get("OVIZ_MCCALLUM_NE_FITS", '/Users/swiggumc/Downloads/ne_grid.fits')
).expanduser()


def _build_mccallum_ne_volumes_for_main_figure():
    if not MCCALLUM_NE_GRID_PATH.exists():
        print(f"Skipping McCallum+2025 electron-density volume; missing FITS cube: {MCCALLUM_NE_GRID_PATH}")
        return []

    return [
        {
            "key": "mccallum-ne",
            "state_key": "mccallum-ne",
            "state_name": "McCallum Hα",
            "name": "McCallum Hα",
            "path": str(MCCALLUM_NE_GRID_PATH),
            "hdu": "PRIMARY",
            "clip_bounds": {"z": [-400.0, 400.0]},
            "max_resolution": 128,
            "max_resolution_cap": 128,
            "opacity": 0.12,
            "samples": 160,
            "alpha_coef": 200,
            "gradient_step": 0.006,
            "stretch": "log10",
            "default_vmin_quantile": 0.90,
            "default_vmax_quantile": 0.9995,
            "colormap": "inferno",
            "unit_label": "cm^-3",
            "visible": False,
            "only_at_t0": True,
            "supports_show_all_times": False,
            "show_all_times": False,
        }
    ]


mccallum_ne_volumes = _build_mccallum_ne_volumes_for_main_figure()

import os
from pathlib import Path as _OvizPath

VERGELY_DUST_PATH = _OvizPath(
    os.environ.get("OVIZ_VERGELY_DUST_FITS", '/Users/swiggumc/Downloads/vergely_3D_Dust.fits')
).expanduser()


def _build_vergely_dust_volumes_for_main_figure():
    if not VERGELY_DUST_PATH.exists():
        print(f"Skipping Vergely 3D dust volume; missing FITS cube: {VERGELY_DUST_PATH}")
        return []

    return [
        {
            "key": "vergely-dust",
            "state_key": "vergely-dust",
            "state_name": "Vergely 3D Dust",
            "name": "Vergely 3D Dust",
            "path": str(VERGELY_DUST_PATH),
            "hdu": "PRIMARY",
            "clip_bounds": {"z": [-400.0, 400.0]},
            "max_resolution": 128,
            "max_resolution_cap": 128,
            "sky_overlay_max_resolution": 128,
            "data_encoding": "png_atlas_uint8",
            "opacity": 0.46,
            "samples": 160,
            "alpha_coef": 200.0,
            "gradient_step": 0.006,
            "stretch": "asinh",
            "vmin": 0.0,
            "vmax": 0.02,
            "default_vmin_quantile": 0.7,
            "default_vmax_quantile": 0.995,
            "colormap": "Greys",
            "unit_label": "mag pc^-1",
            "visible": False,
            "only_at_t0": True,
            "supports_show_all_times": False,
            "show_all_times": False,
        }
    ]


vergely_dust_volumes = _build_vergely_dust_volumes_for_main_figure()

optional_static_volume_state = {
    str(_volume.get("state_key") or _volume.get("key")): {
        "visible": False,
        "showAllTimes": False,
    }
    for _volume in [*mccallum_ne_volumes, *vergely_dust_volumes]
}
optional_static_legend_state = {
    _state_key: False
    for _state_key in optional_static_volume_state
}
fig3d = plot_3d.make_plot(
    time=time_int,
    show=False,
    save_name=save_name,
    static_traces=None,
    static_traces_times=None,
    static_traces_legendonly=True,
    focus_group=None,
    fade_in_time=8,
    fade_in_and_out=False,
    show_age_kde_inset=True,
    age_kde_bandwidth_myr=2,
    include_spiral_arms=False,
    galactic_mode=True,
    show_galactic_guides=True,
    camera_zoom_factor=5,
    show_gc_line=True,
    #cluster_members_file="/Users/swiggumc/Downloads/members-2.csv",
    show_milky_way_model=False,
    renderer="threejs",
    threejs_initial_state={'current_group': 'Clusters', 'click_selection_enabled': False, 'compact_payload_enabled': True, 'compact_widget_payload_enabled': False, 'scene_float_precision': 2, 'active_volume_key': ('supernova-density' if supernova_volumes else 'volume-0'), 'mobile_defer_volumes': False, 'legend_state': ({'volume-0': False, **optional_static_legend_state, 'supernova-density': True} if supernova_volumes else {'volume-0': True, **optional_static_legend_state}), 'volume_state_by_key': ({'volume-0': {'visible': False, 'opacity': 1.0, 'stretch': 'asinh', 'vmax': 0.07}, **optional_static_volume_state, 'supernova-density': {'visible': True}} if supernova_volumes else {'volume-0': {'visible': True, 'opacity': 1.0, 'stretch': 'asinh', 'vmax': 0.07}, **optional_static_volume_state}), 'galaxy_image': True, 'galaxy_image_path': '/Users/swiggumc/Downloads/Top-down_view_of_the_Milky_Way.jpg', 'galaxy_image_size_pc': 40000.0, 'galaxy_image_opacity': 0.35, 'galaxy_image_hide_below_scale_bar_pc': 420.0, 'galaxy_image_fade_start_scale_bar_pc': 700.0, 'sky_dome_enabled': True, 'sky_dome_background_mode': 'live_aladin', 'sky_dome_source': 'aladin', 'sky_dome_projection': 'TAN', 'sky_dome_capture_width_px': 4096, 'sky_dome_capture_height_px': 2048, 'sky_dome_capture_format': 'image/jpeg', 'sky_dome_capture_quality': 0.94, 'sky_dome_radius_pc': 40000.0, 'sky_dome_opacity': 1.0, 'sky_dome_force_visible': False, 'sky_dome_full_opacity_scale_bar_pc': 120.0, 'sky_dome_fade_out_scale_bar_pc': 360.0, 'sky_layers': [{'key': 'P/Mellinger/color', 'label': 'Mellinger Color', 'survey': 'P/Mellinger/color', 'opacity': 1.0, 'visible': True}, {'key': 'P/PLANCK/R2/HFI/color', 'label': 'Planck Dust Emission Color', 'survey': 'P/PLANCK/R2/HFI/color', 'opacity': 1.0, 'visible': True, 'stretch': 'asinh', 'cut_max': 0.07}], 'active_sky_layer_key': 'P/Mellinger/color'},
    enable_sky_panel=True,
    cluster_members_file='/Users/swiggumc/Downloads/members-2.csv',
    show_cluster_members_in_sky=True,
    volumes=[edenhofer_volume, *mccallum_ne_volumes, *vergely_dust_volumes, *supernova_volumes]
)


### OLD PLOTLY

### Make figure
# xyz_widths = (2000, 2000, 400)
# plot_3d = Animate3D(
#     data_collection = traces, 
#     xyz_widths = xyz_widths, 
#     figure_theme = 'dark', 
#     trace_grouping_dict=trace_groupings
#     )
# save_name = '/Users/swiggumc/Desktop/astro_research/radcliffe/cam_website_clone/cam_website/interactive_figures/main_figure.html'


# fig3d = plot_3d.make_plot(
#     time = time_int,
#     show = False, 
#     save_name =save_name, 
#     static_traces = static_traces, 
#     static_traces_times = static_traces_times, 
#     static_traces_legendonly=True,
#     #focus_group =  None,
#     focus_group = None,
#     fade_in_time = 8, # Myr,
#     galactic_mode=True,
#     camera_zoom_factor = 5,
#     include_spiral_arms = False,
#     show_age_kde_inset=True,
#     age_dakde_bandwidth_myr = 1,
#     show_galactic_guides = False, 
#     #coord_system = 'rot',
#     fade_in_and_out = False,
#     show_gc_line = True
# )

# Launch Dash with the new sky panel enabled.
# Click a cluster member to draw the 3D cone footprint and load Aladin imagery.
# app = run_dash_app_in_notebook(
#     figure=fig3d,
#     mode='external',
#     host='127.0.0.1',
#     port=8061,
#     debug=False,
#     title='Main Figure App',
#     height=850,
#     enable_age_filter=False,
#     enable_sky_panel=True,
#     sky_radius_deg=7.0,
#     sky_frame='galactic',
#     sky_survey='P/DSS2/color',
#     cluster_members_file='/Users/swiggumc/Downloads/members-2.csv',
# )

# taurus_names = ['HSC_1318', 'Theia_65', 'Theia_66', 'CWNU_1129', 'Theia_54', 'Theia_93']
# taurus = df_hunt_full.loc[df_hunt_full['name'].isin(taurus_names)]

# # Traces
# taurus_trace = Trace(taurus, data_name = 'Taurus Region Clusters', min_size = 2., max_size = 12, color = 'violet', opacity = 1, marker_style = 'circle', show_tracks = True, size_by_n_stars=False)

# # Redefine some traces
# full_sample_trace = Trace(df_hunt_60, data_name = 'Clusters (< 60 Myr)', min_size = 0, max_size = 12, color = 'grey', opacity = .3, marker_style = 'circle', show_tracks = False, size_by_n_stars=False)
# ap_trace = Trace(ap_family, data_name = 'Alpha Persei Family', min_size = 0., max_size = 12, color = 'violet', opacity = .3, marker_style = 'circle', show_tracks = False, size_by_n_stars=False)
# rw_trace = Trace(rw, data_name = 'Radcliffe Wave Clusters', min_size = 0., max_size = 12, color = 'red', opacity = .3, marker_style = 'circle', show_tracks = False, size_by_n_stars=False)

# traces = TraceCollection([
#     sun_trace,
#     full_sample_trace,
#     taurus_trace,
#     rw_trace,
#     ap_trace
#     ])

# trace_groupings = {
#     "Taurus, Alpha Persei family, and Radcliffe Wave clusters": ['Sun', 'Taurus Region Clusters', 'Radcliffe Wave Clusters', 'Alpha Persei Family'],
# }




# ### Make figure
# xyz_widths = (1001, 1001, 1001)
# plot_3d = Animate3D(
#     data_collection = traces, 
#     xyz_widths = xyz_widths, 
#     figure_theme = 'dark', 
#     trace_grouping_dict=trace_groupings
#     )
# save_name = '/Users/swiggumc/Desktop/astro_research/radcliffe/cam_website_clone/cam_website/interactive_figures/taurus_figure.html'


# fig3d = plot_3d.make_plot(
#     time = time_int,
#     show = False, 
#     save_name =save_name, 
#     static_traces = static_traces, 
#     static_traces_times = static_traces_times, 
#     static_traces_legendonly=False,
#     #focus_group =  None,
#     focus_group = None,
#     fade_in_time = 8, # Myr,
#     galactic_mode=False,
#     #coord_system = 'rot',
#     fade_in_and_out = False,
#     show_gc_line = True
# )
