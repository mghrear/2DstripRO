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

			print("here")
			print(x_errs)
			print(y_errs)
			print(z_vals)

			all_z = np.append(all_z,z_vals)
			all_x_err = np.append(all_x_err,x_errs)
			all_y_err = np.append(all_y_err,y_errs)
            

        
	except:
		pass

# Switch to cm
all_z = all_z*1E-4

print(all_x_err)


for z_low in np.arange(0,1.2,0.2):

	z_high = z_low + 0.2

	#make data cut
	data_cut = (all_z > z_low) & (all_z < z_high)

	print(all_x_err[data_cut]) 
	print("abs. z = ", (z_low+z_high)/2.0 )
	print("x std = ", np.std(all_x_err[data_cut]))
	print("y std = ", np.std(all_y_err[data_cut]))
	print("----------------------------------")
