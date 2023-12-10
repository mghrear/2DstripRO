import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import stripROtools

# Location of data
file_loc = ["/Users/majdghrear/data/VMM3a_SRS/AUG23/UH_DLC/Po210/540Vmesh_bc_40p000_tac_60_ccs_6_cs_3_dt_200_mst_15_spc_1500_dp_400_cr_0p20-2p00_coin_center-of-mass_test.root", "/Users/majdghrear/data/VMM3a_SRS/AUG23/UH_DLC/Po210/540Vmesh_2_bc_40p000_tac_60_ccs_6_cs_3_dt_200_mst_15_spc_1500_dp_400_cr_0p20-2p00_coin_center-of-mass_test.root"]

#file_loc = ["/Users/majdghrear/data/VMM3a_SRS/AUG23/UH_NoDLC/Po210/540Vmesh_bc_40p000_tac_60_ccs_6_cs_3_dt_200_mst_15_spc_1500_dp_400_cr_0p20-2p00_coin_center-of-mass_test.root","/Users/majdghrear/data/VMM3a_SRS/AUG23/UH_NoDLC/Po210/540Vmesh_2_bc_40p000_tac_60_ccs_6_cs_3_dt_200_mst_15_spc_1500_dp_400_cr_0p20-2p00_coin_center-of-mass_test.root"]

# Create pandas data frame of the cluster info
df_cluster = stripROtools.read_root(file_loc, clusters=True, hits=False)

# Define additional columns

# number of hits
df_cluster["nhits"]=df_cluster.apply(lambda row: len(row.strips0)+len(row.strips1) ,axis=1)
df_cluster["electrons_x"] = df_cluster.adc0.apply(lambda x: 6240 * ( x  / 9.0 ) ) # 9 mV/fC is VMM gain setting for x channels, 170mV is the pedestal, 1200mV is th operating voltage 1024 is the number of possible ADC values
df_cluster["electrons_y"] = df_cluster.adc1.apply(lambda x: 6240 * ( x / 4.5 ) ) # 4.5 mV/fC is VMM gain setting for y channels, 170mV is the pedestal, 1200mV is th operating voltage 1024 is the number of possible ADC values
df_cluster["electrons"] = df_cluster.electrons_x + df_cluster.electrons_y

# length on x/y plane in units of strip lengths
df_cluster["L"]=np.sqrt((df_cluster.strips0.apply(np.max) - df_cluster.strips0.apply(np.min))**2 + (df_cluster.strips1.apply(np.max) - df_cluster.strips1.apply(np.min))**2)

# time range on x strips
df_cluster["TR0"] = df_cluster.times1.apply(np.max)-df_cluster.times1.apply(np.min)

# time range on y strips
df_cluster["TR1"] = df_cluster.times1.apply(np.max)-df_cluster.times1.apply(np.min)

#fiducialize clusters on s ingle vmmm combo in xHyL
df_cluster = stripROtools.fiducializeVMM(df_cluster, n_vmm_x=5, n_vmm_y=10, min_hits=5, map="UH")

# Cut for throughgoing events  with over 8 hits that are within 25 degrees of vertical
df_cut = df_cluster.loc[ (df_cluster.TR0 > 1200) & (df_cluster.TR1 > 1200) & (df_cluster.nhits > 8) & (df_cluster.L <= 28) ]
df_cut=df_cut.reset_index(drop=True)

print("Number of remaining events: ", len(df_cut) )


def get_ErrVecs(x_vals,y_vals,z_vals,charges):
	
	X = np.array([x_vals,y_vals,z_vals]).T
	# 1) Center on barycenter
	# Barycenter is the charge-weighted mean position
	x_b = np.sum(X*(charges.reshape(len(charges),1)),axis=0)/np.sum(charges)
	# Shift data to barycenter
	X = X-x_b
	# 2) Find principle axis
	# Use charges for weights
	W = charges.reshape(len(charges),1)
	# Compute weighted covariance matrix
	WCM = ( (W*X).T @ X ) / np.sum(W)
	U1,S1,D1 =  np.linalg.svd(WCM)
	v_PA = np.array([D1[0][0],D1[0][1],D1[0][2]])
	v_PA = np.sign(v_PA[2]) * v_PA
	
	# projection of mean-centered position onto principle axis
	proj = np.array([(X@v_PA)*v_PA[0],(X@v_PA)*v_PA[1],(X@v_PA)*v_PA[2]]).T
	# Mismeasurement vectors
	# The distribution of the x and y values gives us sigma x and sigma y
	err =X-proj
	return x_vals,y_vals,z_vals,charges,v_PA,x_b,err


def fit_gauss(offsets):
	
	# Histogram the charge distribution for fe55 events in the specified time period
        xmin = offsets.min()
        xmax = offsets.max()
        nbins = 20

        hist, bin_edges = np.histogram(offsets,nbins,(xmin,xmax))
        bin_centers = (bin_edges[1:]+bin_edges[:-1])/2.

        # Find non-zero bins in Histogram
        nz = hist>0
        # Get posssion error bars for non-zero bins
        n_err = np.sqrt(hist[nz])

        # Fit Gaussian to binned data
        coeff, covar = curve_fit(gaus, bin_centers[nz], hist[nz], sigma=n_err, absolute_sigma=True, p0=(0,100,0,20))
        # Compute fit (statistical) errors
        perr = np.sqrt(np.diag(covar))

        # If the uncertainty is too high, the fit has failed
        if ( np.absolute(perr[2]) > 0.25*np.absolute(coeff[2])) or ( np.absolute(perr[3]) > 0.25*np.absolute(coeff[3]) ) or np.isnan(perr[2]) or np.isnan(perr[3]) :
            raise Exception("Poor fit")

        if plot == True:
            plt.figure()
            hist, bin_edges,patches = plt.hist(offsets,nbins,(xmin,xmax), color='g',alpha=0.6)
            plt.xlabel("Time Offset")
            plt.ylabel("Count")
            f_opti = gaus(bin_centers,*coeff)
            plt.plot(bin_centers, f_opti, 'r--', linewidth=2, label='curve_fit')
            plt.axvline(3*coeff[3]+coeff[2])
            plt.show()



        return coeff[2], perr[2], coeff[3], perr[3]
    
    except:

        print("-fit failed-")
        return np.nan, np.nan, np.nan, np.nan

all_z = np.array([])
all_x_err = np.array([])
all_y_err = np.array([])


for indx in range(len(df_cut)):

	try:

		event = df_cut.iloc[indx]

		# Make a TrackTools object for the event
		dsp = stripROtools.TrackTools(event = event, gain_x=9, gain_y=4.5, v_drift=8.0)

		# Remove delayed hits (as discussed in the slides)
		dsp.prune_track(T_L = 40, T_H=250)
		dsp.prune_track2(gap=2)
		x,y,z,c = dsp.Reconst3D_v1( mu = -6.68, sigma = 16.4 , n_sigma = 3, plot = False)

		# Only consider tracks with 5 or more points
		if len(x) > 5:

			# Get error vectors
			x_vals,y_vals,z_vals,weights,direction,start, err = get_ErrVecs( x,y,z,c)

			x_errs = err[:,0]
			y_errs = err[:,1]


			all_z = np.append(all_z,z_vals)
			all_x_err = np.append(all_x_err,x_errs)
			all_y_err = np.append(all_y_err,y_errs)
            

        
	except:
		pass

# Switch to cm
all_z = all_z*1E-4


for z_low in np.arange(0,1.2,0.2):

	z_high = z_low + 0.2

	#make data cut
	data_cut = (all_z > z_low) & (all_z < z_high)

	print("abs. z = ", (z_low+z_high)/2.0 )
	print("x std = ", np.std(all_x_err[data_cut]))
	print("y std = ", np.std(all_y_err[data_cut]))
	print("----------------------------------")
