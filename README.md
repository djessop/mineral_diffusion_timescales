# mineral diffusion timescales

Code and sample data for calculating Fe-Mg diffusion timescales in minerals, as described in Metcalf et al. (2021), "Magmatic Processes at La Soufrière de Guadeloupe: Insights from Crystal Studies and Diffusion Timescales for eruption onset", Frontiers in Earth Science, DOI:xxx.  This software is distributed under a GPL-v3.0 licence.

If you use this package please cite the above article in your study.  Citation files (.bib and .ris) are included in the repository.

## Installation 
The code is written in python and is compatible with all versions >3.6.  It can be run from the command line, within an interactive python environment (e.g. ipython), or using the supplied jupyter notebook.

To install, first clone the git repository:
> git clone https://github.com/djessop/mineral_diffusion_timescales.git

or download the zip file.

Second, check that the required packages are installed
> pip install -r diffusion_timescale_modelling/requirements.txt

You may wish to do this within a virtual environment, so as not to affect the function of other projects.
> python3 -m venv [your venv name]



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
├── mineral_diffusion_timescales.ipynb
├── mineral_diffusion_timescales.pdf
├── README.md
└── SEM_traverse_data
    ├── 1010CE
    │   ├── 0111A_A_02-1.xls
    │   ├── 0111A_A_02-2.xls
    │   ├── 0111A_A_02-3.xls
    .
    .
    .
    │   ├── 0111A_B_34-1.xls
    │   ├── 0111A_B_34-2.xls
    │   └── 0111A_B_34-3.xls
    ├── 1657CE
    │   ├── 0913D_A_01-1.xls
    │   ├── 0913D_A_01-2.xls
    │   ├── 0913D_A_01-3.xls
    .
    .
    .
    │   ├── 0913D_C_40-1.xls
    │   ├── 0913D_C_40-2.xls
    │   └── 0913D_C_40-3.xls
    ├── 5680BCE
    │   ├── 1101via_A_07-1.xls
    │   ├── 1101via_A_07-2.xls
    │   ├── 1101via_A_07-3.xls
    .
    .
    .
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
