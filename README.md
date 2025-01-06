# 2DstripRO

A repository of analysis for the comparison of x/y strip readouts instrumented with a PHA setup and VMM/SRS DAQ system

## stripROtools.py

Contains several classes and functions that are generally useful for analysing digitized x/y strip data.

The TrackTools class has methods for 3D reconsting and plotting tracks.

The VMMconfig class manages VMM calibration and configurtion data

## PHA_analysis.ipynb

Contains the analysis using PHA data

## VMM_config.ipynb

A notebook demonstrating how VMM config and calib data are managed

## S_curves.ipynb

Organizes and plots data from VMM S-curve scans

## strip_matching.ipynb

Uses Fe-55 data to obtain info use for matching x and y strips in 3D reconstruction algorithm

## Fe55.ipynb

Contains the analysis of VMM data using Fe-55 source

## Po210.ipynb

Contains the analysis of VMM data using Po-210 source

## Po210_sim.ipynb

Contains the analysis of simulated Po-210 data

# point_res.ipynb

Combines and plots data from Po210.ipynb and Po210_sim.ipynb

## plot_Po210_sim.ipynb and plot_Po210.ipynb
Notebooks used to make plots of simulated and experimental Po-210 emitted alpha tracks

## Po210_sim_v_expected.ipynb

A notebook where simulation parameters can be varied and compared to expeimental results

## test_3D_recon.ipynb

A notebook demonstrating the diffusion suppresion effect from scratch on straight lines of ionization

