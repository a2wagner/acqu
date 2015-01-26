#------------------------------------------------------------------------------
#                 Configuration file for MC data processing
#             (AcquRoot 4v4 user package demonstration analysis)
#------------------------------------------------------------------------------
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Set-up of branches in MC ROOT file as produced by GEANT4 A2 simulation
#         Tree Name     No. of branches   Process Type
Name:     h12           32                MC
#         Branch Name   Buffer size (bytes)
Branch:   npart         4
Branch:   dircos        4096
Branch:   plab          2048
Branch:   elab          2048
Branch:   idpart        2048
Branch:   vertex        12
Branch:   beam          16
Branch:   nhits         4
Branch:   ecryst        2048
Branch:   icryst        2048
Branch:   enai          4
Branch:   vhits         4
Branch:   eveto         2048
Branch:   iveto         2048
Branch:   etot          4
Branch:   eleak         4
Branch:   ntaps         4
Branch:   ictaps        2048
Branch:   ectapsl       2048
Branch:   ectapfs       2048
Branch:   tctaps        2048
Branch:   evtaps        2048
Branch:   nvtaps        4
Branch:   ivtaps        2048
Branch:   tcryst        2048
Branch:   tveto         2048
Branch:   nmwpc         4
Branch:   imwpc         1024
Branch:   mposx         1024
Branch:   mposy         1024
Branch:   mposz         1024
Branch:   emwpc         1024
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ROOT file(s) (GEANT4 output) to be analysed. Add more lines for more files
##Liste am besten erstellen, indem ins Verzeichnis gewechselt wird und die Dateien mit ls -1 angezeigt werden (evtl. mit grep etc. filtern)
#
#TreeFile:	/run/media/wagners/787c785e-bdb0-428f-8160-7c732ab5f090/home/wagners/Desktop/g4_sim.root
TreeFile:	/home/wagners/g4_sim_pi0_gg_00.root
#TreeFile:	/home/wagners/g4_sim_eta_gg_00.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_01.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_02.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_03.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_04.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_05.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_06.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_07.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_08.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_09.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_10.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_11.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_12.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_13.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_14.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_15.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_16.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_17.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_18.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_19.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_20.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_21.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_22.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_23.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_24.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_25.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_26.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_27.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_28.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_29.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_30.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_31.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_32.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_33.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_34.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_35.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_36.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_37.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_38.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_39.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_40.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_41.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_42.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_43.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_44.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_45.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_46.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_47.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_48.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_49.root
#TreeFile:	/data/simulation/background/channels/g4_sim/g4_sim_etap_e+e-g_50.root
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Batch mode operation...redirect output directory for log files etc.
BatchDir: out/
#
#Directory to save any created ROOT files (include final '/')
Directory:  /data/analysis/
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Manual setup of ADC numbers range, scalers, bits-in-ADC numbers. Usually
#needed only where experiment header information is not given, e.g. if the
#input tree file is generated by GEANT4.
#        ADCs    Scaler  Bits
SetADC:  40000   640	 12
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Class name for analysis class to be used in further processing
Analysis:      TA2UserAnalysis
#Configuration file for selected analysis class
AnalysisSetup: AR-Analysis.MC
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
