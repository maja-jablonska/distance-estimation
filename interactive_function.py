# In this file, all relevant functions for the interactive function as well as the interactive function itself which is connected to the widget buttons in the interactive notebook and performs the plotting etc. is defined.

import packages
import functions
from functions import *

# plots the comparison plots to the distance estimation from the published distance catalogue; saves the plot as pdf file
# input: quantiles rmeds, rhis, rlos from interactive distance estimation
#        quantiles rmeds_cat, rhis_cat, rlos_cat from catalogue

def plot_comparison_plots(rmeds,rhis,rlos,rmeds_cat,rhis_cat,rlos_cat,filename_out):

    fig2,ax2 = plt.subplots(3,2,figsize=(7,7))
    fig2.suptitle('Comparison to published distances in gaiaedr3_distance')
    ax2[0,0].scatter(rmeds_cat,rmeds)
    ax2[0,0].set_ylabel('estimated distance $r_{est}$ [pc]')
    ax2[0,0].set_xlabel('distance from catalogue $r_{cat}$ [pc]')
    ax2[0,1].scatter(rmeds_cat,rmeds-rmeds_cat)
    ax2[0,1].set_ylabel('$r_{est}-r_{cat}$ [pc]')
    ax2[0,1].set_xlabel('$r_{cat}$ [pc]')
    sigma_x = ((rhis-rmeds) + (rmeds-rlos))/2
    sigma_y = ((rhis_cat-rmeds_cat) + (rmeds_cat-rlos_cat))/2             
    comp = (rmeds-rmeds_cat)/np.sqrt(sigma_x**2+sigma_y**2)
    
    ax2[1,0].scatter(rmeds_cat,comp)
    ax2[1,0].set_xlabel('$r_{cat}$ [pc]')
    ax2[1,0].set_ylabel('$\\frac{r_{est}-r_{cat}}{\sqrt{\sigma_{est}^2+\sigma_{cat}^2}}$')
    
    ax2[1,1].hist(comp)
    ax2[1,1].set_xlabel('$\\frac{r_{est}-r_{cat}}{\sqrt{\sigma_{est}^2+\sigma_{cat}^2}}$')
    
    ax2[2,0].scatter(rmeds_cat,(rmeds-rmeds_cat)/rmeds_cat)
    ax2[2,0].set_ylabel('$\\frac{r_{est}-r_{cat}}{r_{cat}}$')
    ax2[2,0].set_xlabel('$r_{cat}$ [pc]')
    
    ax2[2,1].remove()
    
    fig2.tight_layout()
    plt.savefig("./results/{}_comparison-plots.pdf".format(filename_out), format="pdf", bbox_inches="tight")
    plt.show()
    
# function that gets called in the interactive function if the input is a csv-file containing (name, parallax, parallax_error, ra, dec). The function returns a pdf with distance estimation plots for each source in the file depending on the parallax, parallax_error, ra, dec as well as a csv-file containing the distance estimation and a csv-file containing the MCMC-chains to the directory /results. 
    
# input:
# model: EDSD or GGD
# filname_in: name of input-csv-file
# filename_out: name of output files (pdf with plots, summary csv and MCMC csv)
# rows_prior_summary: list with information on geometric priors where index hp+1 corresponds to healpixel hp
# HDIprob, rlo, rhi, rplotlo: plotting parameters

def input_csv_parallax(model,filename_out,filename_in,rows_prior_summary,HDIprob,rlo,rhi,rplotlo):
    
    if model == 'Photogeometric':
        print("\033[1;31m Photogeometric model only available for modes 'single' and 'CSVfile(source_id only)'\033[0m")
                
    else:
        # define output pdf file, later contains plots of all stars
        
        p = PdfPages('./results/{}.pdf'.format(filename_out))
        
        # read in the data
        
        with open('./data/{}.csv'.format(filename_in), newline='') as f1:
            reader = csv.reader(f1)
            rows = []
            for row in reader:
                rows.append(row)
            
            # check if header of input file is ['name', 'parallax', 'parallax_error', 'ra', 'dec']
            
            if all(x in rows[0] for x in ['name', 'parallax', 'parallax_error', 'ra', 'dec']):
                
                # create new header for output csv files; one containing all estimated distances and data of the stars and the other containing the 
                # MC chains
                
                header = ['name','HEALpixel level 5','w [mas]', 'wsd [mas]', 'wsd/w','glon for HEALpixel [deg]','glat for HEALpixel [deg]','rest [pc]','rlo [pc]','rhi [pc]','rlen [pc]']
                header_MC = ['name','samples']
                
                # open output csv files and write in the headers
                
                with open('./results/{}_summary.csv'.format(filename_out), 'w', encoding='UTF8') as f2:
                    writer = csv.writer(f2)
                    
                    # write the header
                    
                    writer.writerow(header)
                    
                    with open('./results/{}_MCMCsamples.csv'.format(filename_out), 'w', encoding='UTF8') as f3:
                        
                        writer3 = csv.writer(f3)
                        writer3.writerow(header_MC)
                    
                    # Now a loop is created in which the distance estimation is performed for each source in the csv-file
                    
                        #extract the data from the input file
                        
                        for  i in rows[1:]:
                            
                            name = int(i[rows[0].index('name')])
                            w = float(i[rows[0].index('parallax')])
                            wsd = float(i[rows[0].index('parallax_error')])
                            ra = float(i[rows[0].index('ra')])
                            dec = float(i[rows[0].index('dec')])
                            HP = HEALPix(nside=2**5, order='nested') #level 5
                            hp = HP.lonlat_to_healpix(ra*u.deg,dec*u.deg) 
                            
                            print('name',name)
                            
                            # extract data on prior from the prior summary file
                            
                            if model == 'EDSD':
                                rlen = float(rows_prior_summary[hp+1][10])  #EDSDrlen in csv-file
                            if model == 'GGD':  
                                rlen = float(rows_prior_summary[hp+1][5]) # GGDrlen in csv-file
                                
                            rlen_EDSD = float(rows_prior_summary[hp+1][10])      
                            alpha = float(rows_prior_summary[hp+1][6])
                            beta = float(rows_prior_summary[hp+1][7])
                            glon = float(rows_prior_summary[hp+1][1])
                            glat = float(rows_prior_summary[hp+1][2])
                            
                            
                            # initial guess of distance; use the mode of the EDSD model
                            
                            rInit = float(mode_post3(w=w*1e-3,wsd=wsd*1e-3,rlen = rlen_EDSD,retall = False)) 
                            rStep = 0.75*rInit*min(1/3, abs(wsd/w))
                            Nsamp = int(5e3)
                            Nburnin = int(5e2)
                            
                            # Convert inputs
                            
                            
                            w = 1e-3*w
                            wsd = 1e-3*wsd
                            
                            # get the quantiles of the distance posteriors 
                            #rQuant0: (output dictionary, MCMC samples)
                            #rQuant: only output dictionary
        
                            if model == 'EDSD':
        
                                failMessage = np.nan
                                rQuant0 = quantile_distpost3(w=w,wsd=wsd,rlen=rlen,rInit=rInit,rStep=rStep,Nburnin=Nburnin,Nsamp=Nsamp,probs=np.array([0.5,(1-HDIprob)/2,(1+HDIprob)/2]))
                                rQuant = rQuant0[0]
        
                            if model == 'GGD':
        
                                failMessage = np.nan
                                rQuant0 = quantile_distpost4(w=w,wsd=wsd,rlen=rlen,alpha=alpha,beta=beta,rStep=rStep,rInit=rInit,Nburnin=Nburnin,Nsamp=Nsamp,probs=np.array([0.5,(1-HDIprob)/2,(1+HDIprob)/2]))
                                rQuant = rQuant0[0]

                            # if code = 1: rRes output with estimated quantiles, else change failMessage from nan to string with the message 
                            # if the estimation has worked, rRes containes quantiles and failMessage = nan , otherwise rRes containes nan's and
                            # failMessage is a string 
    
                            if rQuant['code'] == 1:
                                rRes = np.array([rQuant['val'][0],rQuant['val'][1],rQuant['val'][2],2])
                                
                            else:
                                failMessage = rQuant['message']
    
                            if failMessage == rQuant['message']:
                                rRes = np.empty((6,))
                                rRes[:] = np.nan
                                rRes = np.array([rRes[0],rRes[1],rRes[2],rRes[3],rRes[4],0])    
        
                
                            if type(failMessage) == str :
                                print(f"\033[1;31m {failMessage} \033[0m")
                            
                            if rRes[3] == 0:
                                 print('result_flag = 0.Attempting to plot ...')
                            
                            # set upper plotting limits on x-axis rplothi. If at least one of the estimated quantiles is nan: 5*rlen, otherwise 2*r_hi
                            
                            if np.isnan(rRes[0]) or np.isnan(rRes[1]) or np.isnan(rRes[2]) or np.isnan(rRes[3]) :
                                rplothi = 5*rlen
                            else:
                                rplothi = 2*rRes[2]
                            
                            data = [name,hp,1e3*w,1e3*wsd,wsd/w,glon,glat,rRes[0],rRes[1],rRes[2],rlen]
                            
                            writer.writerow(data)
                            
                            if type(rQuant0[1]) == float:   
                                writer3.writerow([name,rQuant0[1]])
                            else:
                                writer3.writerow([name]+rQuant0[1].tolist())
                            
                            # Plots 
                            
                            samp = rQuant0[1]
                            
                            if not np.isnan(samp).all():
                                
                                fig,ax = plt.subplots(2,1,figsize=(5,5),gridspec_kw={'height_ratios': [1,2]})
                                r1 = np.arange(0,len(samp))
                                #plot chains
                                ax[0].plot(r1[::50],samp[::50])
                                ax[0].set_xlabel('samplenumber')
                                Nplot = 1e3
                                s = np.arange(1/(2*Nplot),1/Nplot*(Nplot+1),1/Nplot)
                                rplot = s*(rplothi-rplotlo) + rplotlo
                                
                                # normalize posterior
                                
                                if model == 'GGD':
                                    
                                    Z = scipy.integrate.quad(ud_distpost4,rlo,rhi,args=(w,wsd,rlen,alpha,beta))[0] #normalization factor: 1/Z
                                    dprior = 1/gamma((beta+1)/alpha)*alpha/(rlen**(beta+1))*rplot**beta*np.exp(-(rplot/rlen)**alpha)
                                    if not (Z == 0 or 1/Z == np.inf): #avoid divide by 0 error
                                        dpost  = ud_distpost4(r=rplot, w=w, wsd=wsd, rlen=rlen,alpha=alpha,beta=beta)/Z
                                    
                                if model == 'EDSD':
        
                                    Z = scipy.integrate.quad(ud_distpost3,rlo,rhi,args=(w,wsd,rlen))[0]
                                    dprior = (1/(2*rlen**3))*np.exp(-rplot/rlen)*rplot**2
                                    if not (Z == 0 or 1/Z == np.inf):
                                        dpost  = ud_distpost3(r=rplot, w=w, wsd=wsd, rlen=rlen)/Z
                                
                                ax[1].plot(1e-3*rplot,dprior,color='green',label = f'{model} Prior') #plot prior 
                                ax[1].set_xlim([1e-3*rplotlo,1e-3*rplothi])
                                #ax[1].set_ylim([0,1.05*max(np.array([dprior,dpost]).flatten())])
                                ax[0].set_title('{}'.format(name))
                                ax[1].set_xlabel('distance [kpc]')
                                
                                if not (Z == 0 or 1/Z == np.inf):
                                    ax[1].plot(1e-3*rplot,dpost,color = 'black',label=f'{model} Posterior') #plot posterior, if Z not 0
                                    
                                #plot distance estimates + quantiles   
                                
                                ax[1].axvline(1e-3*rRes[0],label ='$r_{med}$(quantile 0.5)')
                                ax[1].axvline(1e-3*rRes[1], linestyle ='--',label='$r_{lo}$ (quantile 0.159)')
                                ax[1].axvline(1e-3*rRes[2], linestyle ='--',label='$r_{hi}$(quantile 0.841)')
                                
                                
                                ax[1].plot([], [], ' ', label=f'w [mas] = {1e3*w:.3f}')
                                ax[1].plot([], [], ' ', label=f'wsd [mas] = {1e3*wsd:.3f}')       
                                ax[1].plot([], [], ' ', label=f'wsd/w = {wsd/w:.3f}')
                                
                    
                                ax[1].legend(fontsize='6')
                                fig.tight_layout()
                                p.savefig()
                                plt.close()
                        
                    p.close()
            else: 
                print("\033[1;31m csv-input is missing data \033[0m")

# Same as input_csv_parallax, only that the input csv-file only contains a list of source_ids. The parallax w, etc. will then be sourced from gaia within the function. Additionally, a zeropoint correction of the parallaxes will be performed. The output is the same as with input_csv_parallax, only that there additionally will be a comparison plot saved as pdf which compares the estimated distance to the distance from the distance catalogue gaiaedr3_distance. Here, we allow for the photogeometric model as well.
    
def input_csv_sourceid(model,filename_out,filename_in,rows_prior_summary,HDIprob,rlo,rhi,rplotlo):
    
    # define lists of estimated distance for summary statistics     
    
    rmeds = []
    rlos = []        
    rhis = []
    
    # define lists of catalogue distances for summary statistics
    
    rmeds_cat = []
    rlos_cat = []
    rhis_cat = []
    
    # output pdf containing the distance estimation plots
    
    p = PdfPages('./results/{}.pdf'.format(filename_out))
    
    # open input csv-file containing the source_ids
    
    with open('./data/{}.csv'.format(filename_in), newline='') as f1:
        reader = csv.reader(f1)
        rows = []
        for row in reader:
            rows.append(row)
        
        #only proceed if the first entry of the csv-file is the string 'source_id'
        
        if rows[0] == ['source_id']:
            
            #define headers for the output csv-files containing the summary statistics and the MCMC chains
            
            header = ['source_id','HEALpixel level 5','w [mas]', 'wsd [mas]', 'wsd/w','glon for HEALpixel [deg]','glat for HEALpixel [deg]','rest [pc]','rlo [pc]','rhi [pc]','rlen [pc]']
            header_MC = ['source_id','samples']
            
            # open the outputfiles and write the header into them
            
            with open('./results/{}_summary.csv'.format(filename_out), 'w', encoding='UTF8') as f2:
                writer = csv.writer(f2)
                writer.writerow(header)
                with open('./results/{}_MCMCsamples.csv'.format(filename_out), 'w', encoding='UTF8') as f3:
                    writer3 = csv.writer(f3)
                    writer3.writerow(header_MC) 
                   
                # for all source_ids in the csv-file, first source information on the star from gaiadr3
                
                    for i in rows[1:]:
                        job = Gaia.launch_job("select top 10 "
                                      "source_id, parallax, parallax_error,phot_g_mean_mag, nu_eff_used_in_astrometry, pseudocolour, ecl_lat, astrometric_params_solved,bp_rp "
                                      "from gaiadr3.gaia_source "
                                      "where source_id={}".format(int(i[0])))
    
                        r = job.get_results()
                        
                        # len(r)=0 if source_id does not exist within gaia --> print error message!
                
                        if len(r) == 0:
                
                            print("\033[1;31m source_id does not exist! Not added to output file! \033[0m")
                    
                        # if source exists, continue acquiring the data
                        
                        else:
                    
        
                            w= float(r['parallax'])
                            wsd = float(r['parallax_error'])   
                            phot_g_mean_mag = float(r['phot_g_mean_mag'])
                            nu_eff_used_in_astrometry = float(r['nu_eff_used_in_astrometry'])
                            pseudocolour = float(r['pseudocolour'])
                            ecl_lat = float(r['ecl_lat'])
                            astrometric_params_solved = float(r['astrometric_params_solved'])
                            bp_rp = float(r['bp_rp'])
                            
                            
                            #Perform a Zeropoint-correction: 
                            
                            if astrometric_params_solved == 31 or astrometric_params_solved == 95:
                            
                                if  phot_g_mean_mag == np.nan:
                                    wzp = -0.017
                                else:
                                    wzp = zpt.get_zpt(phot_g_mean_mag, nu_eff_used_in_astrometry, pseudocolour, ecl_lat, astrometric_params_solved)
                    
                                #Correct for Zeropoint:
                        
                                w = w - wzp
                    
                            source_id = int(i[0])
                            hp = math.floor(source_id / (2**(35)*4**(12-5)))
                            
                            print('source_id',source_id)
                                    
                            if model == 'EDSD':
                                rlen = float(rows_prior_summary[hp+1][10])  #EDSDrlen in csv-file
                            if model == 'GGD':  
                                rlen = float(rows_prior_summary[hp+1][5]) # GGDrlen in csv-file
                            if model == 'Photogeometric':  
                                rlen = np.nan
                                
                            rlen_EDSD = float(rows_prior_summary[hp+1][10])      
                            alpha = float(rows_prior_summary[hp+1][6])
                            beta = float(rows_prior_summary[hp+1][7])
                            glon = float(rows_prior_summary[hp+1][1])
                            glat = float(rows_prior_summary[hp+1][2])
                            
                            rInit = float(mode_post3(w=w*1e-3,wsd=wsd*1e-3,rlen = rlen_EDSD,retall = False)) 
                            rStep = 0.75*rInit*min(1/3, abs(wsd/w))
                            Nsamp = int(5e3)
                            Nburnin = int(5e2)
                            
                            # Convert inputs
                            
                            
                            w = 1e-3*w
                            wsd = 1e-3*wsd
        
                            if model == 'EDSD':
        
                                failMessage = np.nan
                                rQuant0 = quantile_distpost3(w=w,wsd=wsd,rlen=rlen,rInit=rInit,rStep=rStep,Nburnin=Nburnin,Nsamp=Nsamp,probs=np.array([0.5,(1-HDIprob)/2,(1+HDIprob)/2]))
                                rQuant = rQuant0[0]
        
                            if model == 'GGD':
        
                                failMessage = np.nan
                                rQuant0 = quantile_distpost4(w=w,wsd=wsd,rlen=rlen,alpha=alpha,beta=beta,rStep=rStep,rInit=rInit,Nburnin=Nburnin,Nsamp=Nsamp,probs=np.array([0.5,(1-HDIprob)/2,(1+HDIprob)/2]))
                                rQuant = rQuant0[0]
                            
                            if model == 'Photogeometric':
                            
                                failMessage = np.nan
                                rQuant0 = quantile_distpost5(w=w,wsd=wsd,hp=hp,phot_g_mean_mag=phot_g_mean_mag,bp_rp=bp_rp,pseudocolour=pseudocolour,rInit=rInit,Nsamp=Nsamp,Nburnin=Nburnin,probs=np.array([0.5,0.159,0.841]))
                                rQuant = rQuant0[0]
                        
                            if rQuant['code'] == 1:
                                rRes = np.array([rQuant['val'][0],rQuant['val'][1],rQuant['val'][2],2])
                                
                            else:
                                failMessage = rQuant['message']
    
                            if failMessage == rQuant['message']:
                                rRes = np.empty((6,))
                                rRes[:] = np.nan
                                rRes = np.array([rRes[0],rRes[1],rRes[2],rRes[3],rRes[4],0])    
        
                            #if rRes[3] == 1:
                             #   print("HDI: probability contained, #steps to find:",HDIinfo)
                
                            if type(failMessage) == str :
                                print(f"\033[1;31m {failMessage}\033[0m")
                            
                            if rRes[3] == 0:
                                 print('result_flag = 0.Attempting to plot ...')
        
                            if np.isnan(rRes[0]) or np.isnan(rRes[1]) or np.isnan(rRes[2]) or np.isnan(rRes[3]) :
                                rplothi = 5*rlen
                            else:
                                rplothi = 2*rRes[2]
                            
                            data = [source_id,hp,1e3*w,1e3*wsd,wsd/w,glon,glat,rRes[0],rRes[1],rRes[2],rlen]
                            
                            writer.writerow(data)
                            
                            if type(rQuant0[1]) == float:   
                                writer3.writerow([source_id,rQuant0[1]])
                            else:
                                writer3.writerow([source_id]+rQuant0[1].tolist())
                            
                            
                            #get distances from gaiaedr3_distance catalogue
                            
                            job2 = Gaia.launch_job("select top 10 *"
                                                "from external.gaiaedr3_distance "
                                                "where source_id={}".format(int(i[0])))
                            r2 = job2.get_results()
                            
                            if len(r2['r_med_geo'])==1 :
                                
                                rmeds.append(rRes[0])
                                rlos.append(rRes[1])
                                rhis.append(rRes[2])
                                
                                r_med_geo = float(r2['r_med_geo'])
                                r_lo_geo = float(r2['r_lo_geo'])
                                r_hi_geo = float(r2['r_hi_geo'])
                                
                                r_med_photogeo = float(r2['r_med_photogeo'])
                                r_lo_photogeo = float(r2['r_lo_photogeo'])
                                r_hi_photogeo = float(r2['r_hi_photogeo'])
                                
                                if model == 'Photogeometric':
                                    
                                    rmeds_cat.append(r_med_photogeo)
                                    rlos_cat.append(r_lo_photogeo)
                                    rhis_cat.append(r_hi_photogeo)
                                    
                                else:
                                    
                                    rmeds_cat.append(r_med_geo)
                                    rlos_cat.append(r_lo_geo)
                                    rhis_cat.append(r_hi_geo)
                            
                            # Plots 
                            
                            samp = rQuant0[1]
                            
                            if not np.isnan(samp).all():
                                
                                fig,ax = plt.subplots(2,1,figsize=(5,5),gridspec_kw={'height_ratios': [1,2]})
                                r1 = np.arange(0,len(samp))
                                ax[0].plot(r1[::50],samp[::50])
                                ax[0].set_xlabel('samplenumber')
                                Nplot = 1e3
                                s = np.arange(1/(2*Nplot),1/Nplot*(Nplot+1),1/Nplot)
                                rplot = s*(rplothi-rplotlo) + rplotlo
                                
                                if model == 'GGD':
                                    
                                    Z = scipy.integrate.quad(ud_distpost4,rlo,rhi,args=(w,wsd,rlen,alpha,beta))[0]
                                    dprior = 1/gamma((beta+1)/alpha)*alpha/(rlen**(beta+1))*rplot**beta*np.exp(-(rplot/rlen)**alpha)
                                    if not (Z == 0 or 1/Z == np.inf):
                                        dpost  = ud_distpost4(r=rplot, w=w, wsd=wsd, rlen=rlen,alpha=alpha,beta=beta)/Z
                                    
                                if model == 'EDSD':
        
                                    Z = scipy.integrate.quad(ud_distpost3,rlo,rhi,args=(w,wsd,rlen))[0]
                                    dprior = (1/(2*rlen**3))*np.exp(-rplot/rlen)*rplot**2
                                    if not (Z == 0 or 1/Z == np.inf):
                                        dpost  = ud_distpost3(r=rplot, w=w, wsd=wsd, rlen=rlen)/Z
                            
                                if model == 'Photogeometric':
                                    
                                    dpost = r_ud_distpost5_photogeo(r=rplot,parallax=w,parallax_error=wsd,p=hp,phot_g_mean_mag=phot_g_mean_mag,bp_rp=bp_rp,pseudocolour = pseudocolour)
                                    dpost = dpost/max(dpost)
                                    dprior = r_photogeo_dist_prior(r=rplot,rlen=rlen,beta=beta,alpha=alpha,p=hp,phot_g_mean_mag=phot_g_mean_mag,bp_rp=bp_rp,pseudocolour = pseudocolour)
                                    dprior = dprior/max(dprior)
                                    Z = max(dpost)
                    
                                
                                ax[1].plot(1e-3*rplot,dprior,color='green',label = f'{model} Prior')
                                ax[1].set_xlim([1e-3*rplotlo,1e-3*rplothi])
                                #ax[1].set_ylim([0,1.05*max(np.array([dprior,dpost]).flatten())])
                                ax[0].set_title('source_id = {}'.format(i[0]))
                                ax[1].set_xlabel('distance [kpc]')
                                if not (Z == 0 or 1/Z == np.inf):
                                    ax[1].plot(1e-3*rplot,dpost,color = 'black',label=f'{model} Posterior')
                                    
                                ax[1].axvline(1e-3*rRes[0],label ='$r_{med}$(quantile 0.5)')
                                ax[1].axvline(1e-3*rRes[1], linestyle ='--',label='$r_{lo}$ (quantile 0.159)')
                                ax[1].axvline(1e-3*rRes[2], linestyle ='--',label='$r_{hi}$(quantile 0.841)')
                                
                                if model == 'Photogeometric':
                                    
                                    if not np.isnan(np.array([r_med_photogeo,r_lo_photogeo,r_hi_photogeo])).any():
                                        ax[1].plot([], [], ' ', label= 'Distances from catalogue:')
                                        ax[1].axvline(1e-3*r_med_photogeo,label ='$r_{med_{cat}}$ (quantile 0.5)',color='cornflowerblue')
                                        ax[1].axvline(1e-3*r_lo_photogeo, linestyle ='--',label='$r_{lo_{cat}}$ (quantile 0.159)',color='cornflowerblue')
                                        ax[1].axvline(1e-3*r_hi_photogeo, linestyle ='--',label='$r_{hi_{cat}}$ (quantile 0.841)',color='cornflowerblue')
                                    
                                else:
                                    
                                    if not np.isnan(np.array([r_med_geo,r_lo_geo,r_hi_geo])).any():
                                        
                                        ax[1].plot([], [], ' ', label= 'Distances from catalogue:')
                                        ax[1].axvline(1e-3*r_med_geo,label ='$r_{med_{cat}}$ (quantile 0.5)',color='cornflowerblue')
                                        ax[1].axvline(1e-3*r_lo_geo, linestyle ='--',label='$r_{lo_{cat}}$ (quantile 0.159)',color='cornflowerblue')
                                        ax[1].axvline(1e-3*r_hi_geo, linestyle ='--',label='$r_{hi_{cat}}$ (quantile 0.841)',color='cornflowerblue')
                                
                                ax[1].plot([], [], ' ', label=f'w [mas] = {1e3*w:.3f}')
                                ax[1].plot([], [], ' ', label=f'wsd [mas] = {1e3*wsd:.3f}')       
                                ax[1].plot([], [], ' ', label=f'wsd/w = {wsd/w:.3f}')
                                
                    
                                ax[1].legend(fontsize='6')
                                fig.tight_layout()
                                p.savefig()
                                plt.close()
                    
                p.close()  
        else:
            print("\033[1;31m csv-file has to contain only one table with header 'source_id'; either it contains more than one table or header is wrong \033[0m")
    
    #Comparison plots of estimated distances and catalogue distances for entire csv-input  
    
    if rows[0] == ['source_id']:
        
        rmeds = np.array(rmeds)
        rhis = np.array(rhis)
        rlos = np.array(rlos)
        rmeds_cat = np.array(rmeds_cat)
        rhis_cat = np.array(rhis_cat)
        rlos_cats = np.array(rlos_cat)
        
        plot_comparison_plots(rmeds=rmeds,rhis=rhis,rlos=rlos,rmeds_cat=rmeds_cat,rhis_cat=rhis_cat,rlos_cat=rlos_cat,filename_out=filename_out)

# function that gets called if the input is a single source_id in the respective widget. The input here is a single source_id, for which the required data is then sourced from Gaia. All three types of priors are possible with this setting. There is no input file, but the results are still saved to the output csv- and pdf-files with the custom name. The plot will appear in the notebook as well as the summary statistics.         
        
def input_single(model,source_id,filename_out,rows_prior_summary,HDIprob,rlo,rhi,rplotlo):
    
    if not source_id.isdigit():
        
        source_id = f'{resolve_simbad_to_gaia(source_id)}'
        
    if not source_id.isdigit(): 
        
        print(f"\033[1;31m {source_id} \033[0m") 
    
    else:        
     
    
        job = Gaia.launch_job("select top 10 "
                                        "source_id, parallax, parallax_error,phot_g_mean_mag, nu_eff_used_in_astrometry, pseudocolour, ecl_lat, astrometric_params_solved,bp_rp "
                                        "from gaiadr3.gaia_source "
                                        "where source_id={}".format(int(source_id)))
        
        r = job.get_results()
        
        if len(r) == 0:
        
            print("\033[1;31m source_id does not exist! \033[0m") 
            
        else:
            
            w= float(r['parallax'])
            wsd = float(r['parallax_error'])   
            phot_g_mean_mag = float(r['phot_g_mean_mag'])
            nu_eff_used_in_astrometry = float(r['nu_eff_used_in_astrometry'])
            pseudocolour = float(r['pseudocolour'])
            ecl_lat = float(r['ecl_lat'])
            astrometric_params_solved = float(r['astrometric_params_solved'])
            bp_rp = float(r['bp_rp'])
            
            #Zeropoint:
            
            if astrometric_params_solved == 31 or astrometric_params_solved == 95:
                
                if  phot_g_mean_mag == np.nan:
                    wzp = -0.017
                                
                else:
                    wzp = zpt.get_zpt(phot_g_mean_mag, nu_eff_used_in_astrometry, pseudocolour, ecl_lat, astrometric_params_solved)
                            
                #Correct for Zeropoint:
                w = w - wzp
            
            
            
            source_id = int(r['source_id'])
            hp = math.floor(source_id / (2**(35)*4**(12-5)) )
            print('Gaia DR3 source_id',source_id)
            print('HEALpixel level 5:',hp)
            
            if model == 'EDSD' or model == 'Photogeometric':
                rlen = float(rows_prior_summary[hp+1][10])  #EDSDrlen in csv-file
            if model == 'GGD':  
                rlen = float(rows_prior_summary[hp+1][5]) # GGDrlen in csv-file
            
            rlen_EDSD = float(rows_prior_summary[hp+1][10])      
            alpha = float(rows_prior_summary[hp+1][6])
            beta = float(rows_prior_summary[hp+1][7])
            glon = float(rows_prior_summary[hp+1][1])
            glat = float(rows_prior_summary[hp+1][2])
            
            rInit = float(mode_post3(w=1e-3*w,wsd=1e-3*wsd,rlen = rlen_EDSD,retall = False)) 
            rStep = 0.75*rInit*min(1/3, abs(wsd/w))
            Nsamp = int(5e3)
            Nburnin = int(5e2)
            
            # Convert inputs        
            
            w = 1e-3*w
            wsd = 1e-3*wsd
            if astrometric_params_solved == 31 or astrometric_params_solved == 95:
                wzp = 1e-3*wzp
                
            if model == 'EDSD':
                failMessage = np.nan
                rQuant0 = quantile_distpost3(w=w,wsd=wsd,rlen=rlen,rInit=rInit,rStep=rStep,Nburnin=Nburnin,Nsamp=Nsamp,probs=np.array([0.5,0.159,0.841])) #(1-HDIprob)/2,(1+HDIprob)/2
                rQuant = rQuant0[0]
                
            if model == 'GGD':
                failMessage = np.nan
                rQuant0 = quantile_distpost4(w=w,wsd=wsd,rlen=rlen,alpha=alpha,beta=beta,rStep=rStep,rInit=rInit,Nburnin=Nburnin,Nsamp=Nsamp,probs=np.array([0.5,0.159,0.841]))
                rQuant = rQuant0[0]
                
            if model == 'Photogeometric':
                
                failMessage = np.nan
                rQuant0 = quantile_distpost5(w=w,wsd=wsd,hp=hp,phot_g_mean_mag=phot_g_mean_mag,bp_rp=bp_rp,pseudocolour=pseudocolour,rInit=rInit,Nsamp=Nsamp,Nburnin=Nburnin,probs=np.array([0.5,0.159,0.841]))
                rQuant = rQuant0[0]
                if rQuant['code'] == 1:
                    rFlag = rQuant0[2]
                    
            if rQuant['code'] == 1:
                rRes = np.array([rQuant['val'][0],rQuant['val'][1],rQuant['val'][2],2])
                
            else:
                
                failMessage = rQuant['message']
                
            if failMessage == rQuant['message']:
                rRes = np.empty((6,))
                rRes[:] = np.nan
                rRes = np.array([rRes[0],rRes[1],rRes[2],rRes[3],rRes[4],0])    
                
            #if rRes[3] == 1:
             #   print("HDI: probability contained, #steps to find:",HDIinfo)
    
            if type(failMessage) == str :
                
                print('')
                print(f"\033[1;31m {failMessage} \033[0m")
                print ('')
                
            if rRes[3] == 0:
                 print('result_flag = 0.Attempting to plot ...')
                    
            if np.isnan(rRes[0]) or np.isnan(rRes[1]) or np.isnan(rRes[2]) or np.isnan(rRes[3]) :
                rplothi = 5*rlen
            else:
                rplothi = 2*rRes[2]
                
            # Print summary statistics
            
            if astrometric_params_solved == 31 or astrometric_params_solved == 95:
                print('w[mas]',1e3*(w+wzp))
                print('Zeropoint wzp [mas]', 1e3*wzp)
                print('Zeropointcorrected parallax w[mas]',1e3*w)
            else:
                print('w[mas]',1e3*w)
                print('no zeropoint correction, because astrometric_params_solved is not 31 or 95 ')
            print('wsd/w', wsd/w)
            
            if model == 'GGD' or model == 'Photogeometric':
                print('alpha',alpha)
                print('beta',beta)
            
            if model == 'Photogeometric':
                print('G',phot_g_mean_mag)
                print('colour bp-rp',bp_rp)
                    
            print('')
            print(f'glon [deg] for HEALpixel {hp}:',glon)
            print(f'glat [deg] for HEALpixel {hp}:',glat)
            print('')
            print('rest [pc]:',rRes[0])
            print('rlo [pc]:', rRes[1])
            print('rhi [pc]:', rRes[2])
            print('rlen [pc]:', rlen)
            #print('result_flag', rRes[3]) 
            print('')
            print('MCMC initialization [pc]:', rInit)
            print('MCMC stepsize [pc]:',rStep)
            #print('MCMC Number of iterations[pc]:',Nsamp+Nburnin)
            print('MCMC number of burn-in samples:',Nburnin)
            print('MCMC number of retained iterations:',Nsamp)
            print('')
            
            if model == 'Photogeometric' and rQuant['code'] == 1:
                print('QGmodel:',rFlag[3])
                
            
            header = ['source_id','HEALpixel level 5','w [mas]', 'wsd [mas]', 'wsd/w','glon [deg]','glat [deg]','rest [pc]','rlo [pc]','rhi [pc]','rlen [pc]','result_flag']
            header_MC = ['source_id','samples'] 
            
            # save summary + MCMC chains in file
            
            with open('./results/{}_summary.csv'.format(filename_out), 'w', encoding='UTF8') as f2:
                
                header = ['source_id','HEALpixel level 5','w [mas]', 'wsd [mas]', 'wsd/w','glon [deg]','glat [deg]','rest [pc]','rlo [pc]','rhi [pc]','rlen [pc]','result_flag']
                writer = csv.writer(f2)
                # write the header
                writer.writerow(header)
                data = [source_id,hp,1e3*w,1e3*wsd,wsd/w,glon,glat,rRes[0],rRes[1],rRes[2],rlen,rRes[3]]               
                writer.writerow(data)      
                
            with open('./results/{}_MCMCsamples.csv'.format(filename_out), 'w', encoding='UTF8') as f3:
                writer3 = csv.writer(f3)
                writer3.writerow(header_MC)
                if type(rQuant0[1]) == float:   
                    writer3.writerow([source_id,rQuant0[1]])
                else:
                    writer3.writerow([source_id]+rQuant0[1].tolist())
            
            # Plots 
            
            samp = rQuant0[1]
            
            if not np.isnan(samp).all():
            
                fig,ax = plt.subplots(3,1,figsize=(7,7), gridspec_kw={'height_ratios': [1,1,2]})
                
                r1 = np.arange(0,len(samp))
                
                ax[1].plot(r1[::50],samp[::50])
                ax[1].set_xlabel('samplenumber')
                ax[0].hist(samp,bins=100,density=True,label='MCMC samples')
                r0 = np.linspace(min(samp),max(samp),num=100)
                
                Nplot = 1e3
                s = np.arange(1/(2*Nplot),1/Nplot*(Nplot+1),1/Nplot)
                rplot = s*(rplothi-rplotlo) + rplotlo
                
                #Normalization of the posterior and prior
                
                if model == 'EDSD':
                    Z = scipy.integrate.quad(ud_distpost3,rlo,rhi,args=(w,wsd,rlen))[0]
                    if not (Z == 0 or 1/Z == np.inf):
                        ax[0].plot(r0,ud_distpost3(r=r0,w=w,wsd=wsd,rlen=rlen)/Z,label='normalized Posterior')
                        dpost  = ud_distpost3(r=rplot, w=w, wsd=wsd, rlen=rlen)/Z
                    dprior = (1/(2*rlen**3))*np.exp(-rplot/rlen)*rplot**2
                    
                    
                if model == 'GGD':
                    Z = scipy.integrate.quad(ud_distpost4,rlo,rhi,args=(w,wsd,rlen,alpha,beta))[0]
                    if not (Z == 0 or 1/Z == np.inf):
                        ax[0].plot(r0,ud_distpost4(r=r0,w=w,wsd=wsd,rlen=rlen,alpha=alpha,beta=beta)/Z,label='normalized Posterior')
                        dpost  = ud_distpost4(r=rplot, w=w, wsd=wsd, rlen=rlen,alpha=alpha,beta=beta)/Z
                    dprior = 1/gamma((beta+1)/alpha)*alpha/(rlen**(beta+1))*rplot**beta*np.exp(-(rplot/rlen)**alpha)
                    
                # with the photogeometric model, assume that the normalization factor of the posterior 1/Z is approximately 1/max(posterior)
                
                if model == 'Photogeometric':
                    dpost = r_ud_distpost5_photogeo(r=rplot,parallax=w,parallax_error=wsd,p=hp,phot_g_mean_mag=phot_g_mean_mag,bp_rp=bp_rp,pseudocolour = pseudocolour)
                    dpost = dpost/max(dpost)
                    dprior = r_photogeo_dist_prior(r=rplot,rlen=rlen,beta=beta,alpha=alpha,p=hp,phot_g_mean_mag=phot_g_mean_mag,bp_rp=bp_rp,pseudocolour = pseudocolour)
                    dprior = dprior/max(dprior)
                    Z = max(dpost)
                    
                ax[0].set_xlabel('distance [pc]')
                ax[0].legend()
                ax[2].plot(1e-3*rplot,dprior,color='green',label = f'{model} Prior' )
                ax[2].set_xlim([1e-3*rplotlo,1e-3*rplothi])
                ax[2].set_xlabel('distance [kpc]')
                if not (Z == 0 or 1/Z == np.inf):
                    #ax[2].set_ylim([0,1.05*max(np.array([dprior,dpost]).flatten())])
                    ax[2].plot(1e-3*rplot,dpost,color = 'black',label= f'{model} Posterior')
                    
                ax[2].axvline(1e-3*rRes[0],label ='$r_{med}$ (quantile 0.5)')
                ax[2].axvline(1e-3*rRes[1], linestyle ='--',label='$r_{lo}$ (quantile 0.159)')
                ax[2].axvline(1e-3*rRes[2], linestyle ='--',label='$r_{hi}$ (quantile 0.841)')
                                
            #Compare to distances published in gaiaedr3_distance
                
                print('')
                print('Comparison to published distances in gaiaedr3_distance:')
                print('')
                job2 = Gaia.launch_job("select top 10 *"
                                    "from external.gaiaedr3_distance " 
                                    "where source_id={}".format(int(source_id)))
                r2 = job2.get_results()
                
                r_med_geo = float(r2['r_med_geo'])
                r_lo_geo = float(r2['r_lo_geo'])
                r_hi_geo = float(r2['r_hi_geo'])
                
                r_med_photogeo = float(r2['r_med_photogeo'])
                r_lo_photogeo = float(r2['r_lo_photogeo'])
                r_hi_photogeo = float(r2['r_hi_photogeo'])
                
                
                if model != 'Photogeometric':
                    print('geometric posterior:')
                    print('')
                    print('r_med_geo [pc]:',r_med_geo)
                    print('r_lo_geo [pc]:',r_lo_geo)
                    print('r_hi_geo [pc]:',r_hi_geo)
                    
                    sigma_x = ((rRes[2]-rRes[0]) + (rRes[0]-rRes[1]))/2
                    sigma_y = ((r_hi_geo-r_med_geo) + (r_med_geo-r_lo_geo))/2
                    
                    comp = abs(r_med_geo-rRes[0])/np.sqrt(sigma_x**2+sigma_y**2)
                    
                    print('Deviation:',comp)
                    
                    ax[2].plot([], [], ' ', label='Distances from catalogue:')
                    ax[2].axvline(1e-3*r_med_geo,label ='$r_{med_{cat}}$ (quantile 0.5)',color='cornflowerblue', alpha=0.7)
                    ax[2].axvline(1e-3*r_lo_geo, linestyle ='--',label='$r_{lo_{cat}}$ (quantile 0.159)',color='cornflowerblue', alpha=0.7)
                    ax[2].axvline(1e-3*r_hi_geo, linestyle ='--',label='$r_{hi_{cat}}$ (quantile 0.841) ',color='cornflowerblue', alpha=0.7)
                else:
                    print('photogeometric posterior:')
                    print('')
                    print('r_med_photogeo [pc]:',r_med_photogeo)
                    print('r_lo_photogeo [pc]:',r_lo_photogeo)
                    print('r_hi_photogeo [pc]:',r_hi_photogeo)
                    
                    sigma_x = ((rRes[2]-rRes[0]) + (rRes[0]-rRes[1]))/2
                    sigma_y = ((r_hi_photogeo-r_med_photogeo) + (r_med_photogeo-r_lo_photogeo))/2
                    
                    comp = abs(r_med_photogeo-rRes[0])/np.sqrt(sigma_x**2+sigma_y**2)
                    
                    print('Deviation:',comp)
                    
                    ax[2].plot([], [],' ' ,label='Distances from catalogue:')
                    ax[2].axvline(1e-3*r_med_photogeo,label ='$r_{med_{cat}}$ (quantile 0.5) ',color='cornflowerblue', alpha=0.7)
                    ax[2].axvline(1e-3*r_lo_photogeo, linestyle ='--',label='$r_{lo{cat}}$ (quantile 0.159) ',color='cornflowerblue', alpha=0.7)
                    ax[2].axvline(1e-3*r_hi_photogeo, linestyle ='--',label='$r_{hi_{cat}}$ (quantile 0.841) ',color='cornflowerblue', alpha=0.7)
                    
                    
                ax[2].plot([], [], ' ', label=f'w [mas] = {1e3*w:.3f}')
                ax[2].plot([], [], ' ', label=f'wsd [mas] = {1e3*wsd:.3f}')       
                ax[2].plot([], [], ' ', label=f'wsd/w = {wsd/w:.3f}')
                
                ax[2].legend(fontsize='6')
                
                fig.tight_layout()
                
                plt.savefig("./results/{}.pdf".format(filename_out), format="pdf", bbox_inches="tight")
                plt.show()

# function that gets called when the input are custom values for the parallaxes and prior parameters. With this mode, the sample number, stepsize etc. of the metropolis algorithm can be set manually. 
def input_custom(model,filename_out,rows_prior_summary,HDIprob,rlo,rhi,rplotlo,w,wsd,rlen,alpha,beta,metrop_start,metrop_step,metrop_Nsamp,metrop_Nburnin):
    
    w=w
    wsd=wsd
    rlen = rlen
    alpha = alpha 
    beta = beta
    glon = np.nan # Galactic longitude in degrees (0 to 360)
    glat =  np.nan # Galactic latitude (-90 to +90)
    rInit = metrop_start
    rStep = metrop_step
    Nsamp = metrop_Nsamp
    Nburnin = metrop_Nburnin
    source_id = np.nan
    hp = np.nan
    
    
    # Convert inputs        
            
    w = 1e-3*w
    wsd = 1e-3*wsd
    
    if model == 'EDSD':
        failMessage = np.nan
        rQuant0 = quantile_distpost3(w=w,wsd=wsd,rlen=rlen,rInit=rInit,rStep=rStep,Nburnin=Nburnin,Nsamp=Nsamp,probs=np.array([0.5,0.159,0.841])) #(1-HDIprob)/2,(1+HDIprob)/2
        rQuant = rQuant0[0]
        
    if model == 'GGD':
        failMessage = np.nan
        rQuant0 = quantile_distpost4(w=w,wsd=wsd,rlen=rlen,alpha=alpha,beta=beta,rStep=rStep,rInit=rInit,Nburnin=Nburnin,Nsamp=Nsamp,probs=np.array([0.5,0.159,0.841]))
        rQuant = rQuant0[0]
        
    if model == 'Photogeometric':
        rQuant0 = ({'code':0,'val':np.nan,'message':'Photogeometric model only available for modes "Single" and "CSVfile(source_id only)"' },np.nan)
        rQuant = rQuant0[0]
    if rQuant['code'] == 1:
        rRes = np.array([rQuant['val'][0],rQuant['val'][1],rQuant['val'][2],2])
        
    else:
        
        failMessage = rQuant['message']
        
    if failMessage == rQuant['message']:
        rRes = np.empty((6,))
        rRes[:] = np.nan
        rRes = np.array([rRes[0],rRes[1],rRes[2],rRes[3],rRes[4],0])    
        
    #if rRes[3] == 1:
     #   print("HDI: probability contained, #steps to find:",HDIinfo
    if type(failMessage) == str :
        
        print('')
        print(f"\033[1;31m {failMessage} \033[0m")
        print ('')
        
    if rRes[3] == 0:
         print('result_flag = 0.Attempting to plot ...')
            
    if np.isnan(rRes[0]) or np.isnan(rRes[1]) or np.isnan(rRes[2]) or np.isnan(rRes[3]) :
        rplothi = 5*rlen
    else:
        rplothi = 2*rRes[2]
        
    # Print summary statisti
    print('w [mas]',1e3*w)
    print('wsd [mas]',1e3*wsd)
    print('wsd/w', wsd/w)
    
    if model == 'GGD' or model == 'Photogeometric':
        print('alpha',alpha)
        print('beta',beta)
            
    print('')
    print(f'glon [deg] for HEALpixel {hp}:',glon)
    print(f'glat [deg] for HEALpixel {hp}:',glat)
    print('')
    print('rest [pc]:',rRes[0])
    print('rlo [pc]:', rRes[1])
    print('rhi [pc]:', rRes[2])
    print('rlen [pc]:', rlen)
    #print('result_flag', rRes[3]) 
    print('')
    print('MCMC initialization [pc]:', rInit)
    print('MCMC stepsize [pc]:',rStep)
    #print('MCMC number of iterations:',Nsamp+Nburnin)
    print('MCMC number of burn-in samples:',Nburnin)
    if not np.isnan(rQuant0[1]).all():
        print('MCMC number of retained iterations:',len(rQuant0[1]))
    print('')
    
    
    header = ['source_id','HEALpixel level 5','w [mas]', 'wsd [mas]', 'wsd/w','glon [deg]','glat [deg]','rest [pc]','rlo [pc]','rhi [pc]','rlen [pc]','result_flag']
    header_MC = ['source_id','samples'] 
    
    # save summary + MCMC chains in file
    
    with open('./results/{}_summary.csv'.format(filename_out), 'w', encoding='UTF8') as f2:
        
        header = ['source_id','HEALpixel level 5','w [mas]', 'wsd [mas]', 'wsd/w','glon [deg]','glat [deg]','rest [pc]','rlo [pc]','rhi [pc]','rlen [pc]','result_flag']
        writer = csv.writer(f2)
        # write the header
        writer.writerow(header)
        data = [source_id,hp,1e3*w,1e3*wsd,wsd/w,glon,glat,rRes[0],rRes[1],rRes[2],rlen,rRes[3]]               
        writer.writerow(data)      
        
    with open('./results/{}_MCMCsamples.csv'.format(filename_out), 'w', encoding='UTF8') as f3:
        writer3 = csv.writer(f3)
        writer3.writerow(header_MC)
        if type(rQuant0[1]) == float:   
            writer3.writerow([source_id,rQuant0[1]])
        else:
            writer3.writerow([source_id]+rQuant0[1].tolist())
    
    # Plots 
    
    samp = rQuant0[1]
    
    if not np.isnan(samp).all():
    
        fig,ax = plt.subplots(3,1,figsize=(7,7), gridspec_kw={'height_ratios': [1,1,2]})
        
        r1 = np.arange(0,len(samp))
        
        ax[1].plot(r1[::50],samp[::50])
        ax[1].set_xlabel('samplenumber')
        ax[0].hist(samp,bins=100,density=True,label='MCMC samples')
        r0 = np.linspace(min(samp),max(samp),num=100)
        
        Nplot = 1e3
        s = np.arange(1/(2*Nplot),1/Nplot*(Nplot+1),1/Nplot)
        rplot = s*(rplothi-rplotlo) + rplotlo
        
        
        if model == 'EDSD':
            Z = scipy.integrate.quad(ud_distpost3,rlo,rhi,args=(w,wsd,rlen))[0]
            if not (Z == 0 or 1/Z == np.inf):
                ax[0].plot(r0,ud_distpost3(r=r0,w=w,wsd=wsd,rlen=rlen)/Z,label='normalized Posterior')
                dpost  = ud_distpost3(r=rplot, w=w, wsd=wsd, rlen=rlen)/Z
            dprior = (1/(2*rlen**3))*np.exp(-rplot/rlen)*rplot**2
            
            
        if model == 'GGD':
            Z = scipy.integrate.quad(ud_distpost4,rlo,rhi,args=(w,wsd,rlen,alpha,beta))[0]
            if not (Z == 0 or 1/Z == np.inf):
                ax[0].plot(r0,ud_distpost4(r=r0,w=w,wsd=wsd,rlen=rlen,alpha=alpha,beta=beta)/Z,label='normalized Posterior')
                dpost  = ud_distpost4(r=rplot, w=w, wsd=wsd, rlen=rlen,alpha=alpha,beta=beta)/Z
            dprior = 1/gamma((beta+1)/alpha)*alpha/(rlen**(beta+1))*rplot**beta*np.exp(-(rplot/rlen)**alpha)
            
            
        ax[0].set_xlabel('distance [pc]')
        ax[0].legend()
        ax[2].plot(1e-3*rplot,dprior,color='green',label = f'{model} Prior' )
        ax[2].set_xlim([1e-3*rplotlo,1e-3*rplothi])
        ax[2].set_xlabel('distance [kpc]')
        if not (Z == 0 or 1/Z == np.inf):
            #ax[2].set_ylim([0,1.05*max(np.array([dprior,dpost]).flatten())])
            ax[2].plot(1e-3*rplot,dpost,color = 'black',label= f'{model} Posterior')
            
        ax[2].axvline(1e-3*rRes[0],label ='$r_{med}$ (quantile 0.5)')
        ax[2].axvline(1e-3*rRes[1], linestyle ='--',label='$r_{lo}$ (quantile 0.159)')
        ax[2].axvline(1e-3*rRes[2], linestyle ='--',label='$r_{hi}$ (quantile 0.841)')
                        
            
        ax[2].plot([], [], ' ', label=f'w [mas] = {1e3*w:.3f}')
        ax[2].plot([], [], ' ', label=f'wsd [mas] = {1e3*wsd:.3f}')       
        ax[2].plot([], [], ' ', label=f'wsd/w = {wsd/w:.3f}')
        
        ax[2].legend(fontsize='6')
        
        fig.tight_layout()
        plt.savefig("./results/{}.pdf".format(filename_out), format="pdf", bbox_inches="tight")
        plt.show()
    
    
        
# Define the interactive function. This function gets called each time, the start button is pressed. The input variables are connected to the widgets.     
    
def interactive_function(custom,source_id,filename_in,filename_out,model,w,wsd,rlen,alpha,beta,metrop_start,metrop_step,metrop_Nsamp,metrop_Nburnin,seed):
    np.seterr(divide = 'ignore') #ignore divide by 0 errors
    plt.style.use('ggplot')
    
    # general setups:
    
    HDIprob = 0.6827
    
    np.random.seed(seed) #set the seed for all np.random functions
    
    #Range for normalizing Posterior
    
    rlo = 0
    rhi = 1e5
    
    #Plotting range: lower limit
    
    rplotlo = 0
    
    # open csv-file with information on EDSD prior variables for each healpixel and read in the data into an array where the 
    # rownumber corresponds to the healpixel number
    
    with open('prior_summary.csv', newline='') as prior_summary:
        reader_prior_summary = csv.reader(prior_summary)
        rows_prior_summary = []
        for row in reader_prior_summary:
            rows_prior_summary.append(row)

                                                  
#Input CSVfile (name, parallax, parallax_error, ra, dec)-----------------------------------------------------------------


        if custom == 'CSVfile (name, parallax, parallax_error, ra, dec)':
            
            input_csv_parallax(model,filename_out,filename_in,rows_prior_summary,HDIprob,rlo,rhi,rplotlo)

                    
# Input: CSVfile (only source_id)------------------------------------------------------------------------------------------------------------
        
        if custom == 'CSVfile (source_id only)':
                                
            input_csv_sourceid(model,filename_out,filename_in,rows_prior_summary,HDIprob,rlo,rhi,rplotlo)
            
# Input: single source_id -----------------------------------------------------------------------------------------------------

        if custom == 'Single':
                          
            input_single(model,source_id,filename_out,rows_prior_summary,HDIprob,rlo,rhi,rplotlo)
    
# Input: Custom Parameters: ------------------------------------------------------------------------------------------------------      

        if custom == 'Custom':
                        input_custom(model,filename_out,rows_prior_summary,HDIprob,rlo,rhi,rplotlo,w,wsd,rlen,alpha,beta,metrop_start,metrop_step,metrop_Nsamp,metrop_Nburnin)
