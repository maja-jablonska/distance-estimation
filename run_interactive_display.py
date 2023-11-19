# In this file, all interactive buttons are displayed and connected to the interactive function

import packages
import functions
from functions import *
import interactive_function
from interactive_function import *

#button to select the model

model = widgets.RadioButtons(
    options=['GGD', 'EDSD','Photogeometric'],
    description='Model:',
    disabled=False
)

#field to insert source_id

source_id = widgets.Text(
    value='28147497671505664',
    description='source_id/name',
    disabled=False,
    style = {'description_width': 'initial'}
)

#field to insert name of file containing source_ids or (name, parallax, parallax_error, ra, dec)

filename_in = widgets.Text(
    value='test',
    description='Input .csv-file',
    disabled=False
)

#field to name the output file

filename_out = widgets.Text(
    value='results',
    description='Output .pdf/.csv-file',
    disabled=False,
    style = {'description_width': 'initial'}
)

# Buttons to select the mode

custom=widgets.Select(
    options=['Custom', 'Single', 'CSVfile (source_id only)','CSVfile (name, parallax, parallax_error, ra, dec)'],
    value='Custom',
    # rows=10,
    description='Mode:',
    layout={'width': 'max-content'},
    disabled=False
)

#If custom mode is selected:

#slider for custom parallax

w0 = widgets.FloatSlider(
    value=1.00,
    min=-2.00,
    max=2.00,
    step=0.01,
    description='Parallax [mas]:',
    disabled=False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.2f',
    style = {'description_width': 'initial'},
    
    
)

#textfield connected to parallax-slider

w = widgets.FloatText()
widgets.jslink((w0, 'value'), (w, 'value'))

#slider for custom parallax standard deviation

wsd0 = widgets.FloatSlider(
    value=0.2,
    min=0,
    max=2.00,
    step=0.01,
    description='Parallax error [mas]:',
    disabled=False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.2f',
    style = {'description_width': 'initial'}
    
)

#textfield connected to parallax std-slider

wsd = widgets.FloatText()
widgets.jslink((wsd0, 'value'), (wsd, 'value'))

#rlen slider

rlen0 = widgets.FloatSlider(
    value=200,
    min=0,
    max=1000,
    step=0.01,
    description='Length scale [pc]:',
    disabled=False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.2f',
    style = {'description_width': 'initial'}
)

#textfield connected to parallax rlen-slider

rlen = widgets.FloatText()
widgets.jslink((rlen0, 'value'), (rlen, 'value'))

# alpha-slider

alpha0 = widgets.FloatSlider(
    value=1.00,
    min=0,
    max=3,
    step=0.01,
    description='alpha for GGD prior',
    disabled = False ,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.2f',
    style = {'description_width': 'initial'}
)

#textfield connected to parallax alpha-slider

alpha = widgets.FloatText()
widgets.jslink((alpha0, 'value'), (alpha, 'value'))

#beta-slider

beta0 = widgets.FloatSlider(
    value=2.00,
    min=0,
    max=3.00,
    step=0.01,
    description='beta for GGD prior',
    disabled = False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.2f',
    style = {'description_width': 'initial'}
)

#textfield connected to parallax beta-slider

beta = widgets.FloatText()
widgets.jslink((beta0, 'value'), (beta, 'value'))

#metrop_start slider

metrop_start0 = widgets.FloatSlider(
    value=1000,
    min=0,
    max=2000,
    step=1,
    description='starting value of Metropolis algorithm [pc]:',
    disabled = False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.2f',
    style = {'description_width': 'initial'}
)

#textfield connected to parallax metrop-start_slider

metrop_start = widgets.FloatText()
widgets.jslink((metrop_start0, 'value'), (metrop_start, 'value'))

#Metrop_step slider

metrop_step0 = widgets.FloatSlider(
    value=250,
    min=0,
    max=500,
    step=1,
    description='stepsize of Metropolis algorithm [pc]:',
    disabled = False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.2f',
    style = {'description_width': 'initial'}
)

#textfield connected to parallax metrop-step_slider

metrop_step = widgets.FloatText()
widgets.jslink((metrop_step0, 'value'), (metrop_step, 'value'))

#metrop_Nsamp slider

metrop_Nsamp0 = widgets.IntSlider(
    value = 5000,
    min=0,
    max=10000,
    step=1,
    description='number of posterior samples:',
    disabled = False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.2f',
    style = {'description_width': 'initial'}
)

#textfield connected to parallax metrop_Nsamp_slider

metrop_Nsamp = widgets.IntText()
widgets.jslink((metrop_Nsamp0, 'value'), (metrop_Nsamp, 'value'))

#metrop_Nburnin slider

metrop_Nburnin0 = widgets.IntSlider(
    value = 500,
    min=0,
    max=1000,
    step=1,
    description='number of burn-in samples:',
    disabled = False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.2f',
    style = {'description_width': 'initial'}
)

#textfield connected to metrop_Nburnin slider

metrop_Nburnin = widgets.IntText()
widgets.jslink((metrop_Nburnin0, 'value'), (metrop_Nburnin, 'value'))

# display all widgets

display(model)
display(custom)

display(source_id)
display(filename_in)
display(filename_out)


display(widgets.HBox([w0,w]))
display(widgets.HBox([wsd0,wsd]))
display(widgets.HBox([rlen0,rlen]))
display(widgets.HBox([alpha0,alpha]))
display(widgets.HBox([beta0,beta]))
display(widgets.HBox([metrop_start0,metrop_start]))
display(widgets.HBox([metrop_step0,metrop_step]))
display(widgets.HBox([metrop_Nburnin0,metrop_Nburnin]))
display(widgets.HBox([metrop_Nsamp0,metrop_Nsamp]))




out = widgets.Output()

# submit button

submit_button = widgets.Button(description='start')
display(submit_button)

# submit function: when clicking submit, the simulation is run for all settings that have been made

def submit(button):
    out.clear_output()
    with out:
        interactive_function(custom.value,source_id.value,filename_in.value,filename_out.value,model.value,w.value,wsd.value,rlen.value,alpha.value,beta.value,metrop_start.value,metrop_step.value,metrop_Nsamp.value,metrop_Nburnin.value)
        
# tie submit button to a function
submit_button.on_click(submit)

out 
