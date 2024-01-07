import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from scipy import stats
from scipy.optimize import curve_fit
from scipy.stats import crystalball
import json


# Read the root file outputed by vmmsdat and converts it to a pandas dataframe
# Input is a list containing the locations of the root files and flags for whether cluster and/or hit information should be included
def read_root(files, clusters=True, hits=False):

    if (clusters == True) and (hits == False):

        df_cluster = read_clusters(files)

        return df_cluster

    if (clusters == True) and (hits == True):

        df_cluster = read_clusters(files)
        df_hits = read_hits(files)

        return df_cluster, df_hits

    
    if (clusters == False) and (hits == True):

        df_hits = read_hits(files)

        return df_hits


# Read the cluster info from a root file output from vmmsdat
def read_cluster(file_loc):

	file = uproot.open(file_loc)

	clusters_detector = file['clusters_detector']['clusters_detector']

	dict = {
	'id' : clusters_detector['id'].array(),
	'id0' : clusters_detector['id0'].array(),
	'id1' : clusters_detector['id1'].array(),
	'size0' : clusters_detector['size0'].array(),
	'size1' : clusters_detector['size1'].array(),
	'adc0' : clusters_detector['adc0'].array(),
	'adc1' : clusters_detector['adc1'].array(),
	'pos0' : clusters_detector['pos0'].array(),
	'pos1' : clusters_detector['pos1'].array(),
	'time0' : clusters_detector['time0'].array(),
	'time1' : clusters_detector['time1'].array(),
	'dt0' : clusters_detector['dt0'].array(),
	'dt1' : clusters_detector['dt1'].array(),
	'delta_plane' : clusters_detector['delta_plane'].array(),
	'span_cluster0' : clusters_detector['span_cluster0'].array(),
	'span_cluster1' : clusters_detector['span_cluster1'].array(),
	'strips0' : clusters_detector['strips0'].array(),
	'strips1' : clusters_detector['strips1'].array(),
	'adcs0' : clusters_detector['adcs0'].array(),
	'adcs1' : clusters_detector['adcs1'].array(),
	'times0' : clusters_detector['times0'].array(),
	'times1' : clusters_detector['times1'].array()
	}

	df = pd.DataFrame(data = dict)

	return df

# Read the hit info from a root file output from vmmsdat
def read_hit(file_loc):
    
    file = uproot.open(file_loc)
    hits = file['hits']['hits']
    dict = {
        'id' : hits['id'].array(),
        'det' : hits['det'].array(),
        'plane' : hits['plane'].array(),
        'fec' : hits['fec'].array(),
        'vmm' : hits['vmm'].array(),
        'readout_time' : hits['readout_time'].array(),
        'time' : hits['time'].array(),
        'ch' : hits['ch'].array(),
        'pos' : hits['pos'].array(),
        'bcid' : hits['bcid'].array(),
        'tdc' : hits['tdc'].array(),
        'adc' : hits['adc'].array(),
        'over_threshold' : hits['over_threshold'].array(),
        'chip_time' : hits['chip_time'].array()
        }
    df = pd.DataFrame(data = dict)
    
    return df

# Same as read_cluster and read_hit but these read several root files
read_clusters = lambda  files: pd.concat( [read_cluster(file) for file in files] ,ignore_index=True)
read_hits = lambda  files: pd.concat( [read_hit(file) for file in files] ,ignore_index=True)


# Gaussian Function
def gaus(x, y_off, const, mu, sigma):
    return y_off + const* np.exp(-0.5*((x - mu)/sigma)**2)


# Horizontal fit Function
def horizontal(x, H):
    return H


# Exponential plateau fit function
def exp_plat(x, a, b, c):
	return a * (1.0 - np.e**(-b * (x - c)))


# exponential decay fit function
def invs(x, a, b):
	return (1.0/(a*x)) + b

# linear fit function
def linear(x, a, b):
	return (a * x) + b


# townsend gain eqn fit function
def townsend(X, A, B, m):

    V,t = X # V is the mesh voltage [V] and t is the amplifcation gap [m]

    p = 101325 # p is pressure in Pa
    
    return A * ( (V/t) / p )**m * np.exp( -B * ( ( p / (V/t) )**(1-m) ) )
    
# Creates an object which manages all config and calib informatuion
class VMMconfig:
    
    def __init__(self, strip_map_loc = None, pedestal_loc = None, THL_DAC_loc = None, PLSR_DAC_loc = None, THL_loc = None):
        
        # Store file locations
        self.strip_map_loc = strip_map_loc  # Location of the strip to channel mapping file
        self.pedestal_loc = pedestal_loc  # Location of the pedestal scan output file
        self.THL_DAC_loc = THL_DAC_loc # Location of the global threshold DAC calibration output file
        self.PLSR_DAC_loc = PLSR_DAC_loc # Location of the global Pulser DAC calibration output file
        self.THL_loc = THL_loc # Location of the measured threshold output file
        
        if self.strip_map_loc == None:
            raise Exception("You must provide a strip to channel mapping")
            
        #store strip mapping info into dataframe
        self.StripInfo = pd.DataFrame(columns = ['det', 'plane', 'fec', 'vmm','ch','pos'])
        mapping = json.load(open(self.strip_map_loc))
        
        # Add mapping info to the data frame
        for vmm in mapping["vmm_geometry"]:
            for ch, pos in enumerate(vmm['id0']):
                self.StripInfo = self.StripInfo.append({'det' : vmm['detector'], 'plane' : vmm['plane'], 'fec' : vmm['fec'], 'vmm' : vmm['vmm'], 'ch' :  ch, 'pos' : pos }, ignore_index = True)
        
        # Add pedestal info if it is available
        if self.pedestal_loc != None:
            # Read pedestal info [mV]
            df_ped = pd.read_csv(self.pedestal_loc,usecols = [' fec', 'vmm', 'ch', ' pedestal [mV]']).rename(columns={" fec": "fec", " pedestal [mV]": "pedestal"})
            # Remove rows with header info, switch data type to int
            df_ped = df_ped.loc[df_ped.fec != ' fec'].astype('int32')
            # Convert the VMM slow control fecID to fecID set in firmware based on IP address
            df_ped["fec"] = df_ped["fec"].apply(self.fecIDmap)
            # Add pedestal info to df
            self.StripInfo = self.StripInfo.merge(df_ped, on=['vmm', 'ch', 'fec'], how='left')
            # Create an offline mask keeping only channels with 140mV < pedestal < 200mV
            # These channels are  damaged and do not fire
            self.StripInfo["mask"] = self.StripInfo["pedestal"].apply(lambda x: (140<x) and (x<200) )
            
        # Add threshold info if it is available
        if self.THL_loc != None:
            # Read threshold info [mV]
            df_thres = pd.read_csv(self.THL_loc, usecols = [' fec', 'vmm', 'ch', ' threshold [mV]']).rename(columns={" fec": "fec", " threshold [mV]": "threshold"})
            # Remove rows with header info, switch data type to int
            df_thres = df_thres.loc[df_thres.fec != ' fec'].astype('int32')
            # Convert the VMM slow control fecID to fecID set in firmware based on IP address
            df_thres["fec"] = df_thres["fec"].apply(self.fecIDmap)
            # Add threshold info to df
            self.StripInfo = self.StripInfo.merge(df_thres, on=['vmm', 'ch', 'fec'], how='left')
        
        # Store threshold DAC info if it is given
        if self.THL_DAC_loc != None:
            # read specified columns of the Threshold DAC csv and rename the columns
            self.THL_DAC = pd.read_csv(self.THL_DAC_loc, usecols= [' fec', 'vmm', ' threshold dac setting',' threshold dac measured']).rename(columns={" fec": "fec", " threshold dac setting": "THL_DAC", " threshold dac measured": "THL_mV"})
            # Remove rows with header info, switch data type to int
            self.THL_DAC = self.THL_DAC.loc[(self.THL_DAC.fec != ' fec') ].astype('int32')
            # Remove rows with DAC < 200 as we do not fit those (See Lucian's thesis section 3.1) 
            self.THL_DAC = self.THL_DAC.loc[(self.THL_DAC.THL_DAC >= 200) ]
            # Convert the VMM slow control fecID to fecID set in firmware based on IP address
            self.THL_DAC["fec"] = self.THL_DAC["fec"].apply(self.fecIDmap)

            # Collect DAC and corresponding mV values for each fec/VMM combination as an array
            THL_DACs = self.THL_DAC.groupby(['fec', 'vmm'])['THL_DAC'].apply(np.array).reset_index()
            THL_mVs = self.THL_DAC.groupby(['fec', 'vmm'])['THL_mV'].apply(np.array).reset_index()

            # Merge it togather
            self.THL_DAC = THL_DACs.merge(THL_mVs, on=['vmm', 'fec'], how='left')

            # Compute the slope and offset
            self.THL_DAC["slope"] = self.THL_DAC.apply(lambda row: np.polyfit(row.THL_DAC, row.THL_mV, 1)[0] ,axis=1)
            self.THL_DAC["offset"] = self.THL_DAC.apply(lambda row: np.polyfit(row.THL_DAC, row.THL_mV, 1)[1] ,axis=1)

        # Store threshold DAC info if it is given
        if self.PLSR_DAC_loc != None:         
            # read specified columns of the Pulser DAC csv and rename the columns
            self.PLSR_DAC = pd.read_csv(self.PLSR_DAC_loc, usecols= [' fec', 'vmm', ' pulser dac setting',' pulser dac measured']).rename(columns={" fec": "fec", " pulser dac setting": "PLSR_DAC", " pulser dac measured": "PLSR_mV"})
            # Remove rows with header info, switch data type to int
            self.PLSR_DAC = self.PLSR_DAC.loc[(self.PLSR_DAC.fec != ' fec') ].astype('int32')
            # Remove rows with DAC < 200 as we do not fit those (See Lucian's thesis section 3.1) 
            self.PLSR_DAC = self.PLSR_DAC.loc[(self.PLSR_DAC.PLSR_DAC >= 200) ]
            # Convert the VMM slow control fecID to fecID set in firmware based on IP addres
            self.PLSR_DAC["fec"] = self.PLSR_DAC["fec"].apply(self.fecIDmap)

            # Collect DAC and corresponding mV values for each fec/VMM combination as an array
            PLSR_DACs = self.PLSR_DAC.groupby(['fec', 'vmm'])['PLSR_DAC'].apply(np.array).reset_index()
            PLSR_mVs = self.PLSR_DAC.groupby(['fec', 'vmm'])['PLSR_mV'].apply(np.array).reset_index()

            # Merge it togather
            self.PLSR_DAC = PLSR_DACs.merge(PLSR_mVs, on=['vmm', 'fec'], how='left')

            # Compute the slope and offset
            self.PLSR_DAC["slope"] = self.PLSR_DAC.apply(lambda row: np.polyfit(row.PLSR_DAC, row.PLSR_mV, 1)[0] ,axis=1)
            self.PLSR_DAC["offset"] = self.PLSR_DAC.apply(lambda row: np.polyfit(row.PLSR_DAC, row.PLSR_mV, 1)[1] ,axis=1)

        
    # This method converts the VMMslowcontrol fecID to the fecID set in firmware (based on IP address)
    # This must be set manually, currently there is only 1 fec and it's IP is 2
    def fecIDmap(self,x):
        if x == 1:
            return 2
        else:
            raise Exception("Invalid fecID. Update the fecIDmap method")
            
    # suggests threshold DAC values such that there is equal distance in mV between pedestal and threshold for each VMM
    # target_from_pedestal is the desired distance from pedestal
    def THL_DAC_settings(self, target_from_pedestal):

        if self.THL_DAC_loc == None:
            raise Exception("A threshold DAC calibration must be provided")
            
        if self.pedestal_loc == None:
            raise Exception("A pedestal scan must be provided")

        # Collect pedestal info for each channel
        test_ped = self.StripInfo[["fec","vmm","pedestal","mask"]].copy()
        # Remove damaged and masked channels (channels where pedestal is too high or low)
        test_ped = test_ped.loc[ test_ped["mask"] ].reset_index(drop=True)
        # Get the mean pedestal value per VMM
        test_ped = test_ped.groupby(["fec","vmm"]).mean().reset_index()

        # Collect DAC calibration info VMM
        test_THL = self.THL_DAC.copy()
        # Add the mean pedestal value per VMM to this data frame
        test_THL = test_THL.merge(test_ped, on=['fec','vmm'], how='left')

        # Compute the target threshold as the pedestal + the provided target from pedestal
        test_THL["target_thres_mV"] = test_THL["pedestal"] + target_from_pedestal
        # Use the DAC calibrations to get suggested DAC value for each VMM
        test_THL["threshold_DAC"] = ( test_THL["target_thres_mV"] - self.THL_DAC["offset"])/self.THL_DAC["slope"]
        
        return test_THL[["fec","vmm","pedestal","target_thres_mV","threshold_DAC"]]
    
    
    # suggests threshold DAC values given a target in mV
    def PLSR_DAC_settings(self, target):

        if self.PLSR_DAC_loc == None:
            raise Exception("A Pulser DAC calibration must be provided")
            
        df_suggested_DAC = self.PLSR_DAC[["fec","vmm"]].copy()
        df_suggested_DAC["PLSR_DAC"] = (target - self.PLSR_DAC["offset"])/self.PLSR_DAC["slope"]
        
        return df_suggested_DAC


# Track-level class for reconstructing and visualizing
class TrackTools:
    
    def __init__(self, event, gain_x=9.0, gain_y=4.5, n_strips_x = 512, n_strips_y=512, v_drift=8.0, pitch_x=200, pitch_y=200):
        
        # Hit level information
        self.strips_x = event.strips0
        self.strips_y = event.strips1
        self.ADC_x = event.adcs0
        self.ADC_y = event.adcs1
        # Time is in ns, set lowest time as t=0
        self.times_x = event.times0 - min(min(event.times0),min(event.times1))
        self.times_y = event.times1 - min(min(event.times0),min(event.times1))
        
        # Shaping amplifier gain in mV/fC
        self.gain_x = gain_x
        self.gain_y = gain_y
        
        # Total number of strips in x and y
        self.n_strips_x = n_strips_x
        self.n_strips_y = n_strips_y
        
        # Drift speed as computed by Magboltz in um/ns
        self.v_drift = v_drift

        # Stip pitch is always 200 um
        self.pitch_x = pitch_x
        self.pitch_y = pitch_y
        
    def TimeHistView(self, t_bin = 10):

        t_max = int(max(np.max(self.times_x),np.max(self.times_y)))


        time_edges = np.arange(0, t_max, t_bin)

        print(self.times_x)
        
        plt.figure()
        plt.hist(self.times_x, bins = time_edges, histtype="step",color='k',label= "x strips")
        plt.hist(self.times_y, bins = time_edges, histtype="step", ls='--',color='r',label= "y strips")
        plt.legend()
        plt.xlabel("Time [ns]")
        plt.ylabel("No. Hits")
        plt.show()

    def Strip2DView(self,fullview = True):
                
        # Define empty arrays to store 2D histogram information
        x_vals = np.array([],dtype=int)
        y_vals = np.array([],dtype=int)
        weights = np.array([])
        
        # Loop through x strip data
        for x_hit,adc_x in zip(self.strips_x,self.ADC_x):
            
            x_vals = np.append( x_vals, np.ones(self.n_strips_y,dtype=int)*x_hit )
            y_vals = np.append( y_vals, np.arange(self.n_strips_y) )
            weights = np.append(weights, np.ones(self.n_strips_y,dtype=int)*(6242*adc_x / self.gain_x) )
            
        # Loop through y strip data
        for y_hit,adc_y in zip(self.strips_y,self.ADC_y):
            
            x_vals = np.append( x_vals, np.arange(self.n_strips_x) )
            y_vals = np.append( y_vals, np.ones(self.n_strips_x,dtype=int)*y_hit )
            weights = np.append(weights, np.ones(self.n_strips_x,dtype=int)*(6242*adc_y / self.gain_y) )
            
        
        # Define 2D histogram edges
        x_edges = np.arange(-0.5,self.n_strips_x + 1.5,1)
        y_edges = np.arange(-0.5,self.n_strips_y + 1.5,1)
        
        # Plot the 2D histogram
        plt.figure()
        plt.hist2d(x_vals, y_vals, bins=(x_edges, y_edges), weights=weights, cmap=plt.cm.jet)
        plt.colorbar(label="No. electrons")
        plt.xlabel("Strips x")
        plt.ylabel("Strips y")

        #fullview shows entire readout, otherwise the display zooms into the event
        if fullview == False:            
            Mrange = max(max(self.strips_x)-min(self.strips_x), max(self.strips_y)-min(self.strips_y) )
            xmin = min(self.strips_x)- Mrange
            xmax = max(self.strips_x)+ Mrange
            ymin = min(self.strips_y)- Mrange
            ymax = max(self.strips_y)+ Mrange
            plt.axis([xmin, xmax, ymin, ymax])

        plt.show()

    def Strip2DView_times(self,fullview = True):
                
        # Define empty arrays to store 2D histogram information
        x_vals = np.array([],dtype=int)
        y_vals = np.array([],dtype=int)
        weights = np.array([])
        
        # Loop through x strip data
        for x_hit,time_x in zip(self.strips_x,self.times_x):
            for y_hit,time_y in zip(self.strips_y,self.times_y):

                x_vals = np.append( x_vals, x_hit )
                y_vals = np.append( y_vals, y_hit )
                weights = np.append(weights, (time_x+time_y)/2.0 )
            
        
        # Define 2D histogram edges
        x_edges = np.arange(-0.5,self.n_strips_x + 1.5,1)
        y_edges = np.arange(-0.5,self.n_strips_y + 1.5,1)
        
        # Plot the 2D histogram
        plt.figure()
        plt.hist2d(x_vals, y_vals, bins=(x_edges, y_edges), weights=weights, cmap=plt.cm.jet)
        plt.colorbar(label="Mean time [ns]")
        plt.xlabel("Strips x")
        plt.ylabel("Strips y")

        #fullview shows entire readout, otherwise the display zooms into the event
        if fullview == False:            
            Mrange = max(max(self.strips_x)-min(self.strips_x), max(self.strips_y)-min(self.strips_y) )
            xmin = min(self.strips_x)- Mrange
            xmax = max(self.strips_x)+ Mrange
            ymin = min(self.strips_y)- Mrange
            ymax = max(self.strips_y)+ Mrange
            plt.axis([xmin, xmax, ymin, ymax])

        plt.show()

    def prune_track(self, T_L = 40, T_H=250):

        # For each x hit, compute the time difference to its neighboring strip on the left and right hand side
        Txdiff = np.absolute(np.diff(self.times_x))
        Txdiff_L = np.append(Txdiff[0],Txdiff)
        Txdiff_R = np.append(Txdiff,Txdiff[-1])
        # For each x hit, compute the distance (in strips) to the neighboring strip on the left and right hand side
        Sxdiff = np.absolute(np.diff(self.strips_x))
        Sxdiff_L = np.append(Sxdiff[0],Sxdiff)
        Sxdiff_R = np.append(Sxdiff,Sxdiff[-1])

        # Now find the minimum of time/strips on the left vs right hand side
        xdiff = np.minimum(Txdiff_L/Sxdiff_L,Txdiff_R/Sxdiff_R)

        # For each y hit, compute the time difference to its neighboring strip on the left and right hand side
        Tydiff = np.absolute(np.diff(self.times_y))
        Tydiff_L = np.append(Tydiff[0],Tydiff)
        Tydiff_R = np.append(Tydiff,Tydiff[-1])
        # For each x hit, compute the distance (in strips) to the neighboring strip on the left and right hand side
        Sydiff = np.absolute(np.diff(self.strips_y))
        Sydiff_L = np.append(Sydiff[0],Sydiff)
        Sydiff_R = np.append(Sydiff,Sydiff[-1])

        # Now find the minimum of time/strips on the left vs right hand side
        ydiff = np.minimum(Tydiff_L/Sydiff_L,Tydiff_R/Sydiff_R)

        # Use min time/strip to remove bad hits
        self.strips_x = self.strips_x[ (xdiff>T_L) & (xdiff<T_H)]
        self.strips_y = self.strips_y[ (ydiff>T_L) & (ydiff<T_H)]
        self.ADC_x = self.ADC_x[ (xdiff>T_L) & (xdiff<T_H)]
        self.ADC_y = self.ADC_y[ (ydiff>T_L) & (ydiff<T_H)]
        self.times_x = self.times_x[ (xdiff>T_L) & (xdiff<T_H)]
        self.times_y = self.times_y[ (ydiff>T_L) & (ydiff<T_H)]

    def prune_track2(self, gap=2):
        
        # For each x hit, compute the distance (in strips) to the neighboring strip on the left and right hand side
        xdiffs = np.diff(self.strips_x)
        xdiffs_L = np.append(xdiffs[0],xdiffs)
        xdiffs_R = np.append(xdiffs,xdiffs[-1])

        # Now find the minimum distance on the left vs right hand side
        min_xdiffs = np.minimum(xdiffs_L,xdiffs_R)

        # For each y hit, compute the distance (in strips) to the neighboring strip on the left and right hand side
        ydiffs = np.diff(self.strips_y)
        ydiffs_L = np.append(ydiffs[0],ydiffs)
        ydiffs_R = np.append(ydiffs,ydiffs[-1])

        # Now find the minimum distance on the left vs right hand side
        min_ydiffs = np.minimum(ydiffs_L,ydiffs_R)

        # Use this to remove isolated hits which are usually noise
        self.strips_x = self.strips_x[ min_xdiffs <= gap ]
        self.strips_y = self.strips_y[ min_ydiffs <= gap ]
        self.ADC_x = self.ADC_x[ min_xdiffs <= gap ]
        self.ADC_y = self.ADC_y[ min_ydiffs <= gap ]
        self.times_x = self.times_x[ min_xdiffs <= gap ]
        self.times_y = self.times_y[ min_ydiffs <= gap ]


    def Reconst3D_v0 (self, plot = True):


        # Return Error if any strip fires twice, this method is not meant to handle that
        u,c = np.unique(self.strips_x, return_counts=True)
        repeats = True in (c>1)
        if repeats == True:
            raise Exception("There is a repeated x strip, try another method")

        u,c = np.unique(self.strips_y, return_counts=True)
        repeats = True in (c>1)
        if repeats == True:
            raise Exception("There is a repeated y strip, try another method")
        

        # Define empty arrays to store 2D histogram information
        x_vals = []
        y_vals = []
        weights = []
        avg_times = []


        # Loop through every x and y hit combination
        for x_hit, adc_x, time_x  in zip(self.strips_x, self.ADC_x, self.times_x):
            for y_hit, adc_y, time_y  in zip(self.strips_y, self.ADC_y, self.times_y):

                x_vals += [x_hit]
                y_vals += [y_hit]
                weights += [ (6242.0*adc_x / self.gain_x) + (6242.0*adc_y / self.gain_y) ]
                avg_times += [ (time_x+time_y)/2.0 ]

        # Convert to physical quatities
        x_vals = np.array(x_vals)*self.pitch_x            # Multiply by pitch for physical distance
        y_vals = np.array(y_vals)*self.pitch_y            # Multiply by pitch for physical distance
        weights = np.array(weights)                     # Normalize weights
        weights = weights * ( np.sum(self.ADC_x * 6242.0 / self.gain_x) + np.sum(self.ADC_y * 6242.0 / self.gain_y) ) / np.sum(weights)
        z_vals  = np.array(avg_times) * self.v_drift    # multiply by drift speed for z

        if plot == True:
            # Plot the 3D scatter
            fig = plt.figure()
            ax = Axes3D(fig)

            #set color map
            cm = plt.get_cmap('jet')
            cNorm = matplotlib.colors.Normalize(vmin=min(weights), vmax=max(weights))
            scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
            
            ax.scatter(x_vals, y_vals, z_vals, c=scalarMap.to_rgba(weights),s=300)
            scalarMap.set_array(weights)
            fig.colorbar(scalarMap,label="No. electrons")

            # Force all axis to have equal limits
            extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
            sz = extents[:,1] - extents[:,0]
            centers = np.mean(extents, axis=1)
            maxsize = max(abs(sz))
            r = maxsize/2
            for ctr, dim in zip(centers, 'xyz'):
                getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)

            # Set labels
            ax.set_xlabel('x [um]')
            ax.set_ylabel('y [um]')
            ax.set_zlabel('z [um]')
            
            plt.show()

        return x_vals, y_vals, z_vals, weights
        

    def Reconst3D_v1 (self, mu = 4.1, sigma = 18.4 , n_sigma = 3, plot = True):
        # This 3D reconstruction algorithim only matches x and y hits if they are within a time window specified by mu, sigma, n_sigma
        # After x and y hits are matched, the x ADCs are spread evenly among all matched y hits and vice versa
        # The time is the average of the x hit time and y hit time
        # Unmatched hits are spread along all matched vertices via a time-weighted spread

        # Truth array - contains truth value for x and y hits that fire within the time gap window.
        # i.e. if Tarray_{ij} = True, then the ith x hit and the jth y hit are within the gap window
        # and should be combined, this constitutes an xy-hit
        Tarray = np.abs((np.subtract.outer(self.times_x,self.times_y)-mu) / sigma) < n_sigma

        # If no strips are matched use v0 reconstruction method
        if (True in Tarray) == False:
            print("Warning: None of the hits are matched, running Strip3D_v0 instead")
            return self.Reconst3D_v0(plot)

        else:
            # This counts the number of simultaniously triggering y hits for each x hit
            TCol = np.sum(Tarray,axis=1)*1.0
            # This counts the number of simultaniously triggering x hits for each y hit
            TRow = np.sum(Tarray,axis=0)*1.0
            
            # Throw an error if there are unmatched hits
            # This can be updated later
            if (0 in TCol) or (0 in TRow):
                print("Warning: Unmatched hits. Performing time-weighted spread")

            # Collect unmatched hit info
            # Convert ADC to electron count units
            unmatched_ADCs = np.append(self.ADC_x[ TCol == 0 ] * (6242.0 / self.gain_x), self.ADC_y[ TRow == 0 ] * (6242.0 / self.gain_y))
            # Shift x and y times based on mean offset
            unmatched_times = np.append( self.times_x[ TCol == 0 ] - (mu/2.0) ,  self.times_y[ TRow == 0 ] + (mu/2.0) )

            # Rebuild arrays, ommiting unmatched hits
            # Convert ADC to electron count units
            x_times = self.times_x[ TCol > 0 ]
            ADC_x = self.ADC_x[ TCol > 0 ] * (6242.0 / self.gain_x)
            strips_x = self.strips_x[ TCol > 0 ]
            y_times = self.times_y[ TRow > 0 ]
            ADC_y = self.ADC_y[ TRow > 0 ] * (6242.0 / self.gain_y)
            strips_y = self.strips_y[ TRow > 0 ]
            Tarray = np.abs((np.subtract.outer(x_times,y_times)-mu) / sigma) < n_sigma
            TCol = np.sum(Tarray,axis=1)*1.0
            TRow = np.sum(Tarray,axis=0)*1.0

            # This divides the ADC of the x hit by the number of simultaniously triggering y hits
            ADCx_V = np.divide(ADC_x,TCol)

            # This is a matrix of the x ADC contribution to all xy-hits
            elecx_M = np.multiply(ADCx_V[..., None],Tarray)

            # This divides the ADC of the y hit by the number of simultaniously triggering x hits
            ADCy_V = np.divide(ADC_y,TRow)

            # This is a matrix of the y ADC contribution to all xy-hits
            elecy_M = np.multiply(ADCy_V,Tarray)

            # This is the total ADC assigned to each xy-hit
            elec_M = elecx_M+elecy_M

            # This holds the x strip position for each xy-hit
            Stripx_M = np.multiply(strips_x[..., None],Tarray)

            # This holds the y strip position for each xy-hit
            Stripy_M = np.multiply(strips_y,Tarray)

            # This holds the x time measurment for each xy-hit
            Timex_M = np.multiply(x_times[..., None],Tarray)

            # This holds the y time measurment for each xy-hit
            Timey_M = np.multiply(y_times,Tarray)

            # This holds the average time measurment for each xy-hit
            Time_M = (Timex_M + Timey_M) / 2.0

            # absolute time offsets between matched vertices and unmatched hits
            abs_t_off = np.abs( Time_M-np.tensordot(unmatched_times, Tarray, axes=0) )
            # Really we want to weight by the inverse time difference 
            abs_t_off = np.reciprocal(abs_t_off,where= abs_t_off!=0)

            # Corresponding umatched ADC and time offset normalization factor
            ADC_norm = unmatched_ADCs/abs_t_off.sum(axis=1).sum(axis=1)

            # Multiply togather and sum to get total unmatched ADC contribution for each vertex
            unmatched_contrib = (abs_t_off*np.tensordot(ADC_norm, Tarray, axes=0)).sum(axis=0)

            # Add to ADC matrix
            elec_M += unmatched_contrib

            # Convert to physical quatities
            x_vals = Stripx_M[Tarray]*self.pitch_x              # Multiply by pitch for physical distance
            y_vals = Stripy_M[Tarray]*self.pitch_y              # Multiply by pitch for physical distance
            weights = elec_M[Tarray]                            # Weight is number of electrons
            z_vals  = Time_M[Tarray] * self.v_drift             # multiply by drift speed for z

            if plot == True:
                # Plot the 3D scatter
                fig = plt.figure()
                ax = Axes3D(fig)

                #set color map
                cm = plt.get_cmap('jet')
                cNorm = matplotlib.colors.Normalize(vmin=min(weights), vmax=max(weights))
                scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
                
                ax.scatter(x_vals, y_vals, z_vals, c=scalarMap.to_rgba(weights),s=300)
                scalarMap.set_array(weights)
                fig.colorbar(scalarMap,label="No. electrons")

                # Force all axis to have equal limits
                extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
                sz = extents[:,1] - extents[:,0]
                centers = np.mean(extents, axis=1)
                maxsize = max(abs(sz))
                r = maxsize/2
                for ctr, dim in zip(centers, 'xyz'):
                    getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)

                # Set labels
                ax.set_xlabel('x [um]')
                ax.set_ylabel('y [um]')
                ax.set_zlabel('z [um]')
                
                plt.show()

            return x_vals, y_vals, z_vals, weights
        

# Source: https://matplotlib.org/stable/gallery/images_contours_and_fields/image_annotated_heatmap.html
def heatmap(data, row_labels, col_labels, ax=None, cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=False, bottom=True,
                   labeltop=False, labelbottom=True)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=0, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar

# Source: https://matplotlib.org/stable/gallery/images_contours_and_fields/image_annotated_heatmap.html
def annotate_heatmap(im, data=None, valfmt="{x:.0f}", textcolors=("black", "white"), threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts

# Fiducializes dataframe so that clusters are contained in a specified vmm in x and y
# Must select a map ('UH' or 'UoS')
def fiducializeVMM(df_cluster, n_vmm_x,n_vmm_y, min_hits, map):

    # fiducialize events based on selected vmm combo and map
    if (n_vmm_x == 2) and (n_vmm_y == 10) and (map == "UH"):
        df_cluster["flag"]= df_cluster.apply(lambda row: (min(row.strips0) >= 156 ) & (max(row.strips0) <= 217  ) & (min(row.strips1) >= 156  ) & (max(row.strips1) <=  217 )  ,axis = 1)
    elif (n_vmm_x == 2) and (n_vmm_y == 13) and (map == "UH"):
        df_cluster["flag"]= df_cluster.apply(lambda row: (min(row.strips0) >= 156 ) & (max(row.strips0) <= 217  ) & (min(row.strips1) >= 280  ) & (max(row.strips1) <=  342 )  ,axis = 1)
    elif (n_vmm_x == 5) and (n_vmm_y == 10) and (map == "UH"):
        df_cluster["flag"]= df_cluster.apply(lambda row: (min(row.strips0) >= 230 ) & (max(row.strips0) <= 342  ) & (min(row.strips1) >= 156  ) & (max(row.strips1) <=  217 )  ,axis = 1)
    elif (n_vmm_x == 5) and (n_vmm_y == 13) and (map == "UH"):
        df_cluster["flag"]= df_cluster.apply(lambda row: (min(row.strips0) >= 280 ) & (max(row.strips0) <= 342  ) & (min(row.strips1) >= 280  ) & (max(row.strips1) <=  342 )  ,axis = 1)




    elif (n_vmm_x == 2) and (n_vmm_y == 10) and (map == "UoS"):
        df_cluster["flag"]= df_cluster.apply(lambda row: (min(row.strips0) >= 179 ) & (max(row.strips0) <= 233  ) & (min(row.strips1) >= 64  ) & (max(row.strips1) <=  120 )  ,axis = 1)
    elif (n_vmm_x == 2) and (n_vmm_y == 13) and (map == "UoS"):
        df_cluster["flag"]= df_cluster.apply(lambda row: (min(row.strips0) >= 179 ) & (max(row.strips0) <= 233  ) & (min(row.strips1) >= 122  ) & (max(row.strips1) <=  178 )  ,axis = 1)
    elif (n_vmm_x == 5) and (n_vmm_y == 10) and (map == "UoS"):
        df_cluster["flag"]= df_cluster.apply(lambda row: (min(row.strips0) >= 236 ) & (max(row.strips0) <= 293  ) & (min(row.strips1) >= 64  ) & (max(row.strips1) <=  120 )  ,axis = 1)
    elif (n_vmm_x == 5) and (n_vmm_y == 13) and (map == "UoS"):
        df_cluster["flag"]= df_cluster.apply(lambda row: (min(row.strips0) >= 236 ) & (max(row.strips0) <= 293  ) & (min(row.strips1) >= 122  ) & (max(row.strips1) <=  178 )  ,axis = 1)
    else:
        raise Exception("provide valid map / vmm combo")

    df_fid = df_cluster.loc[ (df_cluster.flag == True) & (df_cluster.nhits >= min_hits) ].reset_index(drop=True)

    return df_fid

# Fiducializes dataframe so that clusters are contained in a specified quadrant
def fiducializeQuadrant(df_cluster, x_loc,y_loc, min_hits):

    # fiducialize events based on selected quadrant
    if (x_loc == "xL") and (y_loc == "yL"):
        df_cluster["flag"]= df_cluster.apply(lambda row: (min(row.strips0) >= 125 ) & (max(row.strips0) <= 249  ) & (min(row.strips1) >= 125  ) & (max(row.strips1) <=  249 )  ,axis = 1)
    elif (x_loc == "xL") and (y_loc == "yH"):
        df_cluster["flag"]= df_cluster.apply(lambda row: (min(row.strips0) >= 125 ) & (max(row.strips0) <= 249  ) & (min(row.strips1) >= 250  ) & (max(row.strips1) <=  374 )  ,axis = 1)
    elif (x_loc == "xH") and (y_loc == "yL"):
        df_cluster["flag"]= df_cluster.apply(lambda row: (min(row.strips0) >= 250 ) & (max(row.strips0) <= 374  ) & (min(row.strips1) >= 125  ) & (max(row.strips1) <=  249 )  ,axis = 1)
    elif (x_loc == "xH") and (y_loc == "yH"):
        df_cluster["flag"]= df_cluster.apply(lambda row: (min(row.strips0) >= 250 ) & (max(row.strips0) <= 374  ) & (min(row.strips1) >= 250  ) & (max(row.strips1) <=  374 )  ,axis = 1)
    else:
        raise Exception("provide valid quadrant")

    df_fid = df_cluster.loc[ (df_cluster.flag == True) & (df_cluster.nhits >= min_hits) ]

    return df_fid

        
# Fit a Crystal Ball function to fe55 events
def fitCB(df, plot=True):

    try:
        # Get gain values
        gain = df.gain
        # Keep only gain entries with z-score < 3 (exclude outlier which may be cosmic tracks or nuclear recoils)
        gain =  gain[(np.abs(stats.zscore(gain)) < 3)]

        xmin = 0
        xmax = gain.max()
        nbins = 50

        hist, bin_edges = np.histogram(gain,nbins,(xmin,xmax))
        bin_centers = (bin_edges[1:]+bin_edges[:-1])/2.
        # Find Non-zero bins in Histogram
        nz = hist>0
        # Get error bars for bins
        n_err = np.sqrt(hist[nz])

        # Create numpy Histogram, use density this time
        hist2, bin_edges2 = np.histogram(gain,nbins,(xmin,xmax), density = True)
        bin_centers2 = (bin_edges2[1:]+bin_edges2[:-1])/2.

        # Guess mu as bin_center with most hits
        mu_guess = bin_centers2[np.argmax(hist2)]

        # Find Non-zero bins in Histogram
        nz2 = hist2>0
        # Get error bars for bins
        n_err2 = (np.sqrt(hist[nz])/hist[nz]) * hist2[nz2] # Fractional error times hist value

        # Define Range and Fit :
        coeff, covar = curve_fit(crystalball.pdf, bin_centers2[nz2], hist2[nz2], sigma=n_err2, absolute_sigma=True, p0=(1, 2,mu_guess,1600))

        f_opti = crystalball.pdf(bin_centers,*coeff)

        perr = np.sqrt(np.diag(covar))

        if np.absolute(perr[2]) > np.absolute(coeff[2]):
            raise Exception("Poor fit")


        if plot == True:
            plt.figure()
            hist2, bin_edges2, patches2 = plt.hist(gain,nbins,(xmin,xmax), density = True, color='g',alpha=0.6)
            bin_centers2 = (bin_edges2[1:]+bin_edges2[:-1])/2.
            plt.xlabel("Gain")
            plt.ylabel("Count")
            plt.plot(bin_centers, f_opti, 'r--', linewidth=2, label='curve_fit')
            plt.show()

        charge_sharing = 1.0*np.sum(df.electrons_x)/np.sum(df.electrons_y)


        return coeff[0], perr[0], coeff[1], perr[1], coeff[2], perr[2], coeff[3], perr[3], charge_sharing

    except:
        print("-fit failed-")
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan


# Fit a Gaussian to x/y time offset distribution 
# There are two methods for obtaining the x/y time offset distribution "max_ADC" and "mean_time"
def fit_offset(df, method = "max_ADC", plot=True):

    try:
        # collect time offsets
        if method == "max_ADC":
            offsets = df.maxADC_offset
        elif method == "mean_time":
            offsets = df.mean_offset
        else:
            raise Exception("Select Method: max_ADC or mean_time")


        # Histogram the charge distribution for fe55 events in the specified time period
        xmin = offsets.min()
        xmax = offsets.max()
        nbins = 50

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



# Fit a Gaussian to x/y time offset distribution 
# There are two methods for obtaining the x/y time offset distribution "max_ADC" and "mean_time"
def fit_horizontal(x_vals, y_vals, y_errs):

    
    try:
        # Only fit the non nan entries
        x_vals = np.array(x_vals)[~np.isnan( np.array(y_vals))]
        y_errs = np.array(y_errs)[~np.isnan( np.array(y_vals))]
        y_vals = np.array(y_vals)[~np.isnan( np.array(y_vals))]

        # Fit Gaussian to binned data
        coeff, covar = curve_fit(horizontal, x_vals, y_vals, sigma=y_errs, absolute_sigma=True, p0=(0))
        # Compute fit (statistical) errors
        perr = np.sqrt(np.diag(covar))

        return coeff[0], perr[0]
    
    except:
        print("-fit failed-")
        return np.nan, np.nan