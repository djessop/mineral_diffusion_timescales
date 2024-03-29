# mineral diffusion timescales

Code and sample data for calculating Fe-Mg diffusion timescales in minerals, as described in Metcalf et al. (2021), "Magmatic Processes at La Soufrière de Guadeloupe: Insights from Crystal Studies and Diffusion Timescales for eruption onset", Frontiers in Earth Science, DOI:10.3389/feart.2021.617294.  This software is distributed under a GPL-v3.0 licence.

If you use this package please cite the above article in your study.  Citation files (.bib and .ris) are included in the repository.

## Installation 
The code is written in python and is compatible with all versions >3.6.  It can be run from the command line, within an interactive python environment (e.g. ipython), or using the supplied jupyter notebook.

To install, first clone the git repository:

```bash
git clone https://github.com/djessop/mineral_diffusion_timescales.git
```

or download the zip file.

Second, check that the required packages are installed

```bash
pip install -r diffusion_timescale_modelling/requirements.txt
```

Note that the "--user" flag may be required depending on the user's level of access.  You may wish to do this within a virtual environment, so as not to affect the function of other projects (in which case the "--user" flag is unnecessary),

```bash
python3 -m venv [your venv name]
```


## Development

Further development of this code is planned in the near future.  This will include modules for diffusion models in various mineral systems (olivine, opx...)  Please contact the author of the code (d.jessop@opgc.fr) to discuss adding additional systems or making modifications to existing ones.


## Usage

The core of the diffusion_timescale_modelling code is the fitting of an analytical solution for concentration profile to SEM traverse data.  In doing this, we can determine a timescale for the diffusion processes that created the profile.  Hence, to run the code one needs data, and we have provided our SEM data for several eruptions of la Soufrière volcano, Guadeloupe as reported in Metcalfe et al (2021), hereto refered to as Met2021.

The code is situated in the diffusion_timescale_modelling/diffusion_timescale_modelling.py file.  To use it, do the usual:

```python
from diffusion_timescale_modelling.diffusion_timescale_modelling import <function_name>
```

The module provides the following functions (see the individual function documentation for help):
- opx_model
- diffusion_coeff
- temp_from_diff_coeff
- analytical_soln
- fit_wrapper
- coefficient_determination
- eval_ts_data
- unique_sample_names
- sorted_data_to_df
- run_model_fitting
- plot_data_model
- calc_fit_plot
- read_raw_data

For example, to calculate the diffusion coefficient for a given absolute temperature, T, (in K) and O2 partial pressure, fO2, (in Pa) we would use the "diffusion_coeff" function, which itself makes a call to "opx_model":
```python
print(diffusion_model(T, fO2))
```

Note that fO2 takes a default value of 6.31e-7 Pa, which is the case for the Mg-Fe analyses in Met2021.

Going further, to estimate the diffusion timescale and other parameters via a fit of the empirical data, we would use the "fit_wrapper" function, which itself requires the data, an initial guess at the solution and bounds for the solution domain as inputs.  Here is an example of how to do this, using the 

```python
filename = 'SEM_traverse_data/1010CE/0111A_A_10-3.xls'
data = pd.read_excel(filename, sheet_name='raw')
## the 'raw' sheet contains 'distance' and 'greyscale' columns
x, y = data[['distance', 'greyscale']].values.T
## The initial guess requires a fairly precise value of the sample
## temperature, obtained using the opx-melt geothermometer model of
## Putirka (2008) and contained in the same excel workbook
wb   = xlrd.open_workbook(filename) 
ws   = wb.sheet_by_name('Dan, WH37 processing, usabl (2)')
T    = ws.cell(3, 13).value
D    = diffusion_coeff(T)
Dlow = diffusion_coeff(T - 30)
Dupp = diffusion_coeff(T + 30)
## Initial guess containing Cmax, Cmin, mu, tau, D.  See Met2021 for details
p0 = [y.max(), y.min(), x.mean(), 2e6, D]
## bounds for solution domain: unconstrained except for timescale (+ve
## values only!), and diffusion coefficient:
lower_bounds = [-np.inf, -np.inf, -np.inf, 0, Dlow]
upper_bounds = [ np.inf,  np.inf,  np.inf, np.inf, Dupp]
popt, pcov, (timescale, ts_sigma, Test, diffusion) = fit_wrapper(
    x, y, p0, bounds=[lower_bounds, upper_bounds])
```

In future versions, where other diffusion models and mineral systems are planned, the "function_to_fit" argument may be set accordingly.  Currently, it defaults to the opx/Fe-Mg system ("analytical_solution").

With the above done, it may be interesting to look at a plot of the raw data and model solution.  This can be accomplished using "plot_data_model":

```python
data, R2, (fig, ax) = plot_data_model(filename, popt, pcov, savefig=False)
```

### Bulk runs and quality-of-fit assessment
It might be interesting to batch run the entire database in one go.  Assuming that the structure is as per this archive ('SEM_traverse_data/<eruption_name>/<sample_file>.xls', see below), this can be achieved via

```python
## use the glob module to find all the xls file under the archive
filenames = sorted(glob('SEM_traverse_data/**/*.xls', recursive=True))
df = run_model_fitting(filenames)
## apply the data quality and 'goodness-of-fit' assessments
sorted_df = sorted_data_to_df(df)

## Write each eruption as a tab in an excel file
writer = pd.ExcelWriter('sorted_timescales.xlsx',
    engine='xlsxwriter')

eruptions = sorted(list(set(sorted_df['eruption'])))
for eruption in eruptions: 
    subdf = sorted_df[sorted_df['eruption'] == eruption]
    subdf.to_excel(writer, sheet_name=eruption, index=False) 
    
writer.save()
```

The output is a couple of pandas dataframes ('df' and 'sorted_df') containing the measured and estimated (fitted) melt temperatures, timescale and its standard error, and the R2 (Pearson's correlation coefficient) value.  Currently, this function is not configured for running analyses in parallel, though this is planned for future versions.  For the ~400 samples in the Met2021 study, this takes only a couple of minutes on a computer with a reasonably fast processor.  You will also find an excel file named 'sorted_timescales.xlsx' in the current working directory containing the fitted parameters for each sample, with separate worksheets for each eruption.  


### Population timescales

TO DO!!!


## File structure
```
.
├── citaions
|   ├── Metcalfe_etal_2021_Frontiers.bib
|   ├── Metcalfe_etal_2021_Frontiers.ris
├── diffusion_timescale_modelling
│   ├── diffusion_timescale_modelling.py
│   ├── __init__.py
│   └── requirements.txt
├── ml_stats
|   ├── mle_gamma.py
|   ├── __init__.py
├── model_fits   
├── mineral_diffusion_timescales.ipynb
├── mineral_diffusion_timescales.pdf
├── README.md
└── SEM_traverse_data
    ├── 1010CE
    │   ├── 0111A_A_02-1.xls
    │   ├── 0111A_A_02-2.xls
    │   ├── 0111A_A_02-3.xls
    .   .
    .   .
    .   .
    │   ├── 0111A_B_34-1.xls
    │   ├── 0111A_B_34-2.xls
    │   └── 0111A_B_34-3.xls
    ├── 1657CE
    │   ├── 0913D_A_01-1.xls
    │   ├── 0913D_A_01-2.xls
    │   ├── 0913D_A_01-3.xls
    .   .
    .   .
    .   .
    │   ├── 0913D_C_40-1.xls
    │   ├── 0913D_C_40-2.xls
    │   └── 0913D_C_40-3.xls
    ├── 5680BCE
    │   ├── 1101via_A_07-1.xls
    │   ├── 1101via_A_07-2.xls
    │   ├── 1101via_A_07-3.xls
    .   .
    .   .
    .   .
    │   ├── 1101via_B_20_1.xls
    │   ├── 1101via_B_20_2.xls
    │   └── 1101via_B_20_3.xls
    └── 341CE
        ├── 1904A-01-1.xls
        ├── 1904A-01-2.xls
        ├── 1904A-01-3.xls
        .
	.
	.
        ├── 1904A-83-1.xls
        ├── 1904A-83-2.xls
        └── 1904A-83-3.xls
```
