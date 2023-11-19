# Interactive-Distance-Estimation
The Jupyter Notebook "Interactive Distance Estimation.ipynb" can be used to plot different prior models as well as their respective posteriors and estimated distances and distance errors. 

First you choose your model. There currently are three options:

1. GGD: uses generalized gamma distribution as prior
2. EDSD: uses exponentially decreasing density distribution as prior
3. Photogeometric: uses combination of GGD prior and information on colour and apparent magnitude 

Next, you coose your type of input.
1. Custom
3. Single
5. CSVfile (source_id only)
6. CSVfile (name, parallax, parallax_error, ra, dec)

With "Custom", you can use the sliders and respective fields to manually type in the inputs for your model.The sliders are limited but the boxes next to them are not, so if you want to exceed the maximum or minimum value of the sliders, type your input in the textfield. The fields "source_id" and "Input .csv-file" are only used for the other input types and are ignored here. The custom input is only available for the  geometric models (GGD and EDSD), not the photogeometric one. The parameters "alpha" and "beta" are only used if you chose the GGD prior, since they are fixed for the EDSD prior. Then you can choose a name for your output files with "Output .pdf/.csv-file". There are three outputs with this type which are being saved in the "results"-folder: 
- 'your_name'_MCMCsamples.csv_: csv-file with MCMC-samples 
- 'your_name'_summary.csv_: csv-file with summary statistics
- 'your_name'.pdf: pdf-file with plot
If you have everything adjusted, click "start" to get your outputs.

With "Single", you can choose a single source_id or the name of the star in the field "source_id/name". You can again choose a name for your output. All the other fields are being ignored. The parameters needed are automatically queried from Gaia. The distance inference includes a automatical zeropoint correction of the parallax. On the resulting plot, you can see a second distace esimate, which is the one from the published distance catalogue (blue).

With "CSVfile (source_id only)" you can enter a csv-file containing a table with only source_id's and header "source_id". After selecting your name and klicking "start", it produces the same output as with "single", only for all the sources in the file. In addition to that, in saves comparison plots in 'your_name'_comparison-plots.pdf_, in which the inferred distances are compared to the ones from the catalogue. 

With "CSVfile (name, parallax, parallax_error, ra, dec)", you can enter a csv-file containing the columns listed (the name does not automatically have to be the source_id). With this, no data is queried by gaia and no zeropoint correction is being done (you can use your own correction). There is no comparison to the catalogue distances, as the output does not depend on the source_id. This only works for the geometric priors (GGD and EDSD) because for the Photogeometric, you would need more information. 

In the "results" folder, there is already some data from the test-table for the input "CSVfile (source_id only)" provided in the "data"-folder. You should delete this data before entering your own.

Requirements: 

- R version 4.0.4 (2021-02-15)
- Python 3.9.13 (main, Aug 25 2022, 23:26:10)
- Jupyter Notebook or JupyterLab

Python Packages: 

numpy, matplotlib, ipywidgets, IPython, scipy, math, astroquery, astropy, astropy_healpix, gaiadr3-zeropoint, rpy2

R libraries: 

data.table, bit64, fields, mvtnorm, PolynomF
