import ROOT, collections,sys,os
sys.path.append('./')
from collections import OrderedDict
from TIMBER.Analyzer import analyzer, HistGroup
from TIMBER.Tools.Common import CompileCpp
from TIMBER.Tools.Plot import *
import helpers
ROOT.gROOT.SetBatch(True)

###########################################
# Establish some global variables for use #
###########################################
plotdir = 'plots/'

# US, Europe/Asia, global
#redirector = 'root://cmsxrootd.fnal.gov/'
#redirector = 'root://xrootd-cms.infn.it/'
#redirector = 'root://cms-xrd-global.cern.ch/'

# Use this to access FNAL EOS files at CERN
redirector = 'root://cmseos.fnal.gov/'
rootfile_path = '/store/user/cmsdas/2026/long_exercises/BstarTW/rootfiles'
config = 'bstar_config.json'
if not os.path.exists(plotdir):
    os.makedirs(plotdir)

#################################################
# Create structs for MET flags and trigger info #
#################################################
# MET flags - https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
flags = ["Flag_goodVertices",
        "Flag_globalSuperTightHalo2016Filter",
        "Flag_HBHENoiseFilter",
        "Flag_HBHENoiseIsoFilter",
        "Flag_EcalDeadCellTriggerPrimitiveFilter",
        "Flag_BadPFMuonFilter"
        #"Flag_ecalBadCalibReducedMINIAODFilter"  # Still work in progress flag, may not be used
        ]
triggers = {
        '16': ["HLT_PFHT800","HLT_PFHT900","HLT_PFJet450"],
        '17': ["HLT_PFHT1050","HLT_PFJet500","HLT_AK8PFJet380_TrimMass30","HLT_AK8PFJet400_TrimMass30"],
        '18': ["HLT_PFHT1050","HLT_PFJet500","HLT_AK8PFJet380_TrimMass30","HLT_AK8PFJet400_TrimMass30"]
}

########################################
# Variables we want to create and plot #
########################################
# Dictionary will have form {'RDataFrame column name' : 'LaTeX title'}
varnames = {
        'lead_tau32'             : '#tau_{32}^{jet0}',
        'sublead_tau32'          : '#tau_{32}^{jet1}',
        'lead_tau21'             : '#tau_{21}^{jet0}',
        'sublead_tau21'          : '#tau_{21}^{jet1}',
        'nbjet_loose'            : 'loosebjets',
        'nbjet_medium'           : 'mediumbjets',
        'nbjet_tight'            : 'tightbjets',
        'lead_jetPt'             : 'p_{T}^{jet0}',
        'sublead_jetPt'          : 'p_{T}^{jet1}',
        'lead_eta'               : '#eta^{jet0}',
        'sublead_eta'            : '#eta^{jet1}',
        'lead_phi'               : '#phi^{jet0}',
        'sublead_phi'            : '#phi^{jet1}',
        'deltaphi'               : '#Delta#phi_{jet0,jet1}',
        'deltaeta'               : '#Delta#eta_{jet0,jet1}',
        'invariantMass'          : 'm_{tW} [GeV]',
        'lead_softdrop_mass'     : 'm_{SD}^{jet0}',
        'sublead_softdrop_mass'  : 'm_{SD}^{jet1}',
        'lead_deepAK8_TvsQCD'    : 'Deep AK8 TvsQCD^{jet0}',
        'sublead_deepAK8_TvsQCD' : 'Deep AK8 TvsQCD^{jet1}',
        'lead_deepAK8_WvsQCD'    : 'Deep AK8 WvsQCD^{jet0}',
        'sublead_deepAK8_WvsQCD' : 'Deep AK8 WvsQCD^{jet1}',
}

############################################
# Define functions for the event selection #
############################################
def select(setname, year):
    '''Function to perform the event selection on a specified dataset by: 
         (1) Applying MET filters and trigger selection to dataset
         (2) Identifying events with at least two back-to-back FatJets
         (3) Performing kinematic cuts on the resulting candidate jets
       After performing selection, the results will be saved to a ROOT snapshot and plots will be generated.
    Args:
        setname  (str): name of input dataset (signal, background)
        year     (str): 16, 17, 18
    '''
    # Initialize TIMBER analyzer
    file_path = '{redirector}{rootfile_path}/{setname}_bstar{year}.root'.format(
        redirector=redirector, rootfile_path=rootfile_path, setname=setname, year=year
    )
    a = analyzer(file_path)

    # Determine normalization weight
    if not a.isData:
        norm = helpers.getNormFactor(setname,year,config)
    else:
        norm = 1

    # Book actions on the RDataFrame
    a.Cut('filters',a.GetFlagString(flags))
    a.Cut('trigger',a.GetTriggerString(triggers[year]))         # Apply different triggers based on the year
    a.Define('jetIdx','hemispherize(FatJet_phi, FatJet_jetId)') # need to calculate if we have two jets (with Id) that are back-to-back
    a.Cut('nFatJets_cut','nFatJet > max(jetIdx[0],jetIdx[1])')  # If we don't do this, we may try to access variables of jets that don't exist! (leads to seg fault)
    a.Cut("hemis","(jetIdx[0] != -1)&&(jetIdx[1] != -1)")       # cut on that calculation

    # Having determined which events have two candidate jets meeting our criteria, let's make a collection specific to them
    # Then, every event will have a two-element long column Dijet_<variable> corresponding to the values of the jets which passed our back-to-back criteria
    a.SubCollection('Dijet', 'FatJet', 'jetIdx', useTake=True)

    # We can now cut on the values of the Dijet collection and simply index 0 (lead) or 1 (sublead) to get the information about the appropriate jet
    a.Define('deltaphi','abs(hardware::DeltaPhi(Dijet_phi[0], Dijet_phi[1]))')
    a.Define('deltaeta','Dijet_eta[0] - Dijet_eta[1]')
    a.Cut('pt_cut',   'Dijet_pt[0] > 400 && Dijet_pt[1] > 400')
    a.Cut('eta_cut',  'abs(Dijet_eta[0]) < 2.4 && abs(Dijet_eta[1]) < 2.4')
    a.Cut('deltaphi_cut', 'deltaphi > 1.57079632679')
    a.Cut('mjet_cut', 'Dijet_msoftdrop[0] > 50 && Dijet_msoftdrop[1] > 50')
    a.Define('lead_vector',    'hardware::TLvector(Dijet_pt[0], Dijet_eta[0], Dijet_phi[0], Dijet_msoftdrop[0])')
    a.Define('sublead_vector', 'hardware::TLvector(Dijet_pt[1], Dijet_eta[1], Dijet_phi[1], Dijet_msoftdrop[1])')
    a.Define('invariantMass','hardware::InvariantMass({lead_vector,sublead_vector})')
    a.Cut('mtw_cut','invariantMass > 1200')

    # Now, we can define the variables we're interested in plotting (see varnames dictionary in global definitions above)
    a.Define('lead_tau32',    'Dijet_tau2[0] > 0 ? Dijet_tau3[0]/Dijet_tau2[0] : -1') # Conditional to make sure tau2 != 0 for division
    a.Define('sublead_tau32', 'Dijet_tau2[1] > 0 ? Dijet_tau3[1]/Dijet_tau2[1] : -1') # condition ? <do if true> : <do if false>
    a.Define('lead_tau21',    'Dijet_tau1[0] > 0 ? Dijet_tau2[0]/Dijet_tau1[0] : -1') # Conditional to make sure tau2 != 0 for division
    a.Define('sublead_tau21', 'Dijet_tau1[1] > 0 ? Dijet_tau2[1]/Dijet_tau1[1] : -1') # condition ? <do if true> : <do if false>
    a.Define('lead_deepAK8_TvsQCD',    'Dijet_deepTag_TvsQCD[0]')
    a.Define('sublead_deepAK8_TvsQCD', 'Dijet_deepTag_TvsQCD[1]')
    a.Define('lead_deepAK8_WvsQCD',    'Dijet_deepTag_WvsQCD[0]')
    a.Define('sublead_deepAK8_WvsQCD', 'Dijet_deepTag_WvsQCD[1]')
    bcut = []
    if year == '16' :
        bcut = [0.2217,0.6321,0.8953]
    elif year == '17' :
        bcut = [0.1522,0.4941,0.8001]
    elif year == '18' :
        bcut = [0.1241,0.4184,0.7571]
    a.Define('nbjet_loose',  'Sum(Jet_btagDeepB > '+str(bcut[0])+')') # DeepCSV loose WP
    a.Define('nbjet_medium', 'Sum(Jet_btagDeepB > '+str(bcut[1])+')') # DeepCSV medium WP
    a.Define('nbjet_tight',  'Sum(Jet_btagDeepB > '+str(bcut[2])+')') # DeepCSV tight WP
    a.Define('lead_jetPt',   'Dijet_pt[0]')
    a.Define('sublead_jetPt','Dijet_pt[1]')
    a.Define('lead_eta',     'Dijet_eta[0]')
    a.Define('sublead_eta',  'Dijet_eta[1]')
    a.Define('lead_phi',     'Dijet_phi[0]')
    a.Define('sublead_phi',  'Dijet_phi[1]')
    a.Define('lead_softdrop_mass',   'Dijet_msoftdrop[0]')
    a.Define('sublead_softdrop_mass','Dijet_msoftdrop[1]')
    a.Define('norm',str(norm))

    # Before finishing up, create plots of the variables stored in the varnames dictionary
    print('Plotting the following variables:')
    outFile = ROOT.TFile.Open('rootfiles/{}_{}_selection.root'.format(setname, year),'RECREATE')
    outFile.cd()
    # Book a group to save the histograms
    hists = HistGroup('{}_{}'.format(setname, year))
    for varname in varnames.keys():
        print('\t{}'.format(varname))
        histname = '{}_{}_{}'.format(setname, year, varname)
        # Arguments for binning that you would normally pass to a TH1
        if "nbjet" in varname :
            hist_tuple = (histname,histname, 10,0,10)
        elif "tau" in varname :
            hist_tuple = (histname,histname,20,0,1)
        elif "Pt" in varname :
            hist_tuple = (histname,histname,30,400,2000)
        elif varname == "deltaphi" :
            hist_tuple = (histname,histname,30,0,3.2)
        elif "phi" in varname :
            hist_tuple = (histname,histname,30,-3.2,3.2)
        elif "eta" in varname :
            hist_tuple = (histname,histname,30,-5.0,5.0)
        elif "softdrop_mass" in varname :
            hist_tuple = (histname,histname,30,0,300)
        elif varname == "invariantMass" :
            hist_tuple = (histname,histname,30,1200,4000)
        else:
            hist_tuple = (histname,histname,20,0,1)
        # Project dataframe into a histogram (hist name/binning tuple, variable to plot from dataframe, weight)
        hist = a.GetActiveNode().DataFrame.Histo1D(hist_tuple,varname,'norm')
        hist.GetValue()                 # this gets the actual TH1 instead of a pointer to the TH1. Here is when all the booked actions are performed, so may take a while for larger datasets (e.g. QCD)
        hists.Add(varname, hist)        # add the TH1 to our HistGroup

    # Now, perform TH1.Write() on all TH1s in our HistGroup
    hists.Do('Write')

    # For fun, print the TIMBER node tree for a visualization of your selection
    a.PrintNodeTree('{}/{}_{}_selection_tree.dot'.format(plotdir,setname, year), verbose=True)

    # Finally, safely close the analyzer and rootfile
    outFile.Close()
    a.Close()

if __name__ == "__main__":
    # Take in user input from the command line, specifically the input dataset name and year
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-s', type=str, dest='setname',
                        action='store', required=True,
                        help='Setname to process. E.g. ttbar, signalLH2000, singletop_tW, QCDHT700, etc...')
    parser.add_argument('-y', type=str, dest='year',
                        action='store', required=True,
                        help='Year of set (16, 17, 18).')
    args = parser.parse_args()

    # Compile some of the C++ macros we'll need for our selection
    CompileCpp("TIMBER/Framework/include/common.h")
    CompileCpp('bstar.cc')      # Contains hemispherize() function for identifying back-to-back jets

    # Run our selection script.
    select(args.setname, args.year)
