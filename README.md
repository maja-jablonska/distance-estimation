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

With "Custom", you can use the sliders and respective fields to manually type in the inputs for your model. The fields "source_id" and "Input .csv-file" are only used for the other input types and are ignored here. The custom input is only available for the  geometric models (GGD and EDSD), not the photogeometric one. The parameters "alpha" and "beta" are only used if you chose the GGD prior, since they are fixed for the EDSD prior. Then you can choose a name for your output files with "Output .pdf/.csv-file". There are three outputs with this type which are being saved in the "results"-folder: 
- 'your_name'_MCMCsamples.csv_: csv-file with MCMC-samples 
- 'your_name'_summary.csv_: csv-file with summary statistics
- 'your_name'.pdf: pdf-file with plot
If you have everything adjusted, click "start" to get your outputs.

With "Single", you can choose a single source_id in the field "source_id". You can again choose a name for your output. All the other fields are being ignored. The parameters needed are automatically queried from Gaia. The distance inference includes a automatical zeropoint correction of the parallax. On the resulting plot, you can see a second disnace esimate, which is the one from the published distance catalogue.
