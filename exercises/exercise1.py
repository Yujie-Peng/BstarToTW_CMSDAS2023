'''
   Exercise 3: Event pre-selection.
   Applies loose preselection cuts motivated by detector acceptance and the
   b*->tW signal topology, then plots the relevant signal MC variables.
'''
import ROOT, sys, os
sys.path.append('./')
from optparse import OptionParser
from collections import OrderedDict

from TIMBER.Analyzer import analyzer, HistGroup
from TIMBER.Tools.Common import CompileCpp
from TIMBER.Tools.Plot import *
import helpers

ROOT.gROOT.SetBatch(True)

# ── Command-line options ──────────────────────────────────────────────────────
parser = OptionParser()
parser.add_option('-y', '--year', metavar='YEAR', type='string', action='store',
                  default='', dest='year', help='Year (16, 17, 18)')
parser.add_option('--select', metavar='BOOL', action='store_true',
                  default=False, dest='select',
                  help='Run the selection and save histograms. Without this flag, '
                       'read saved histograms and make plots.')
(options, args) = parser.parse_args()

# ── Global configuration ──────────────────────────────────────────────────────
plotdir      = 'plots/'
redirector   = 'root://cmseos.fnal.gov/'
rootfile_path = '{}/store/user/cmsdas/2026/long_exercises/BstarTW/rootfiles'.format(redirector)
config       = 'bstar_config.json'

if not os.path.exists(plotdir):
    os.makedirs(plotdir)

# Compile C++ helpers once at import time
CompileCpp("TIMBER/Framework/include/common.h")
CompileCpp('bstar.cc')

# ── Sample definitions ────────────────────────────────────────────────────────
signal_names = ['signalLH2000']
names  = {'signalLH2000': 'b*_{LH} 2000 GeV'}
colors = {'signalLH2000': ROOT.kBlue}

# ── MET filters ───────────────────────────────────────────────────────────────
flags = [
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
]

# ── Triggers (year-dependent) ─────────────────────────────────────────────────
triggers = {
    '16': ["HLT_PFHT800", "HLT_PFHT900", "HLT_PFJet450"],
    '17': ["HLT_PFHT1050", "HLT_PFJet500",
           "HLT_AK8PFJet380_TrimMass30", "HLT_AK8PFJet400_TrimMass30"],
    '18': ["HLT_PFHT1050", "HLT_PFJet500",
           "HLT_AK8PFJet380_TrimMass30", "HLT_AK8PFJet400_TrimMass30"],
}

# ── Variables to plot ─────────────────────────────────────────────────────────
varnames = OrderedDict([
    ('lead_jetPt',             'p_{T}^{jet0} [GeV]'),
    ('sublead_jetPt',          'p_{T}^{jet1} [GeV]'),
    ('lead_eta',               '#eta^{jet0}'),
    ('sublead_eta',            '#eta^{jet1}'),
    ('lead_phi',               '#phi^{jet0}'),
    ('sublead_phi',            '#phi^{jet1}'),
    ('lead_softdrop_mass',     'm_{SD}^{jet0} [GeV]'),
    ('sublead_softdrop_mass',  'm_{SD}^{jet1} [GeV]'),
    ('lead_tau21',             '#tau_{21}^{jet0}'),
    ('sublead_tau21',          '#tau_{21}^{jet1}'),
    ('lead_tau32',             '#tau_{32}^{jet0}'),
    ('sublead_tau32',          '#tau_{32}^{jet1}'),
    ('lead_deepAK8_WvsQCD',    'Deep AK8 WvsQCD^{jet0}'),
    ('sublead_deepAK8_WvsQCD', 'Deep AK8 WvsQCD^{jet1}'),
    ('lead_deepAK8_TvsQCD',    'Deep AK8 TvsQCD^{jet0}'),
    ('sublead_deepAK8_TvsQCD', 'Deep AK8 TvsQCD^{jet1}'),
    ('nbjet_loose',            'Loose b-jets'),
    ('nbjet_medium',           'Medium b-jets'),
    ('nbjet_tight',            'Tight b-jets'),
    ('deltaphi',               '#Delta#phi_{jet0,jet1}'),
    ('deltaeta',               '#Delta#eta_{jet0,jet1}'),
    ('invariantMass',          'm_{tW} [GeV]'),
])

# ── Selection function ────────────────────────────────────────────────────────
def select(setname, year):
    ROOT.ROOT.EnableImplicitMT(2)

    file_path = '%s/%s_bstar%s.root' % (rootfile_path, setname, year)
    a = analyzer(file_path)

    # Normalization weight (xsec * lumi / N_gen for MC, 1 for data)
    norm = helpers.getNormFactor(setname, year, config) if not a.isData else 1.

    # ── 1. Event-level quality filters and trigger ────────────────────────────
    a.Cut('filters', a.GetFlagString(flags))
    a.Cut('trigger', a.GetTriggerString(triggers[year]))

    # ── 2. Identify two back-to-back fat jets ─────────────────────────────────
    # hemispherize() (bstar.cc) returns ordered pT indices of two fat jets in
    # opposite phi hemispheres that also pass jetId requirements
    a.Define('jetIdx', 'hemispherize(FatJet_phi, FatJet_jetId)')
    a.Cut('nFatJets_cut', 'nFatJet > max(jetIdx[0],jetIdx[1])')
    a.Cut('hemis', '(jetIdx[0] != -1)&&(jetIdx[1] != -1)')

    # Build a two-element Dijet sub-collection (Dijet_pt, Dijet_eta, ...)
    # useTake=True tells TIMBER to use ROOT::VecOps::Take(FatJet_*, jetIdx)
    a.SubCollection('Dijet', 'FatJet', 'jetIdx', useTake=True)

    # ── 3. Kinematic pre-selection ────────────────────────────────────────────
    a.Define('deltaphi', 'abs(hardware::DeltaPhi(Dijet_phi[0],Dijet_phi[1]))')
    a.Define('deltaeta', 'Dijet_eta[0] - Dijet_eta[1]')
    a.Cut('pt_cut',  'Dijet_pt[0] > 400 && Dijet_pt[1] > 400')
    a.Cut('eta_cut', 'abs(Dijet_eta[0]) < 2.4 && abs(Dijet_eta[1]) < 2.4')
    a.Cut('deltaphi_cut', 'deltaphi > 1.57079632679')
    a.Cut('msd_cut', 'Dijet_msoftdrop[0] > 50 && Dijet_msoftdrop[1] > 50')

    # ── 4. Invariant mass of the tW system ────────────────────────────────────
    a.Define('lead_vector',    'hardware::TLvector(Dijet_pt[0],Dijet_eta[0],Dijet_phi[0],Dijet_msoftdrop[0])')
    a.Define('sublead_vector', 'hardware::TLvector(Dijet_pt[1],Dijet_eta[1],Dijet_phi[1],Dijet_msoftdrop[1])')
    a.Define('invariantMass',  'hardware::InvariantMass({lead_vector,sublead_vector})')
    a.Cut('mtw_cut', 'invariantMass > 1200')

    # ── 5. N-subjettiness ratios ─────────────────────────────────────────────
    # tau32 = tau3/tau2: 3-prong (top) tagger
    a.Define('lead_tau32',    'Dijet_tau2[0]>0 ? Dijet_tau3[0]/Dijet_tau2[0] : -1')
    a.Define('sublead_tau32', 'Dijet_tau2[1]>0 ? Dijet_tau3[1]/Dijet_tau2[1] : -1')
    # tau21 = tau2/tau1: 2-prong (W) tagger
    a.Define('lead_tau21',    'Dijet_tau1[0]>0 ? Dijet_tau2[0]/Dijet_tau1[0] : -1')
    a.Define('sublead_tau21', 'Dijet_tau1[1]>0 ? Dijet_tau2[1]/Dijet_tau1[1] : -1')

    # ── 6. DeepAK8 scores ────────────────────────────────────────────────────  
    a.Define('lead_deepAK8_TvsQCD',    'Dijet_deepTag_TvsQCD[0]')
    a.Define('sublead_deepAK8_TvsQCD', 'Dijet_deepTag_TvsQCD[1]')
    a.Define('lead_deepAK8_WvsQCD',    'Dijet_deepTag_WvsQCD[0]')
    a.Define('sublead_deepAK8_WvsQCD', 'Dijet_deepTag_WvsQCD[1]')

    # ── 7. B-jet counting (AK4 DeepCSV) ─────────────────────────────────────
    bcut = {'16': [0.2217, 0.6321, 0.8953],
            '17': [0.1522, 0.4941, 0.8001],
            '18': [0.1241, 0.4184, 0.7571]}[year]
    a.Define('nbjet_loose',  'Sum(Jet_btagDeepB > %s)' % bcut[0])
    a.Define('nbjet_medium', 'Sum(Jet_btagDeepB > %s)' % bcut[1])
    a.Define('nbjet_tight',  'Sum(Jet_btagDeepB > %s)' % bcut[2])

    # ── 8. Remaining kinematic variables ────────────────────────────────────
    a.Define('lead_jetPt',            'Dijet_pt[0]')
    a.Define('sublead_jetPt',         'Dijet_pt[1]')
    a.Define('lead_eta',              'Dijet_eta[0]')
    a.Define('sublead_eta',           'Dijet_eta[1]')
    a.Define('lead_phi',              'Dijet_phi[0]')
    a.Define('sublead_phi',           'Dijet_phi[1]')
    a.Define('lead_softdrop_mass',    'Dijet_msoftdrop[0]')
    a.Define('sublead_softdrop_mass', 'Dijet_msoftdrop[1]')
    a.Define('norm', str(norm))

    # Print TIMBER cut-flow tree (visualize with graphviz)
    try:
        a.PrintNodeTree(plotdir + 'signal_tree_ex3.dot', verbose=True)
    except Exception as e:
        print('[warning] PrintNodeTree skipped: %s' % e)

    # ── 9. Book histograms ───────────────────────────────────────────────────
    out = HistGroup('%s_%s' % (setname, year))
    for varname in varnames.keys():
        histname = '%s_%s_%s' % (setname, year, varname)
        if 'nbjet' in varname:
            hist_tuple = (histname, histname, 10, 0, 10)
        elif 'tau' in varname:
            hist_tuple = (histname, histname, 20, 0, 1)
        elif 'Pt' in varname:
            hist_tuple = (histname, histname, 30, 400, 2000)
        elif varname == 'deltaphi':
            hist_tuple = (histname, histname, 30, 0, 3.2)
        elif 'phi' in varname:
            hist_tuple = (histname, histname, 30, -3.2, 3.2)
        elif 'eta' in varname:
            hist_tuple = (histname, histname, 30, -5.0, 5.0)
        elif 'softdrop' in varname:
            hist_tuple = (histname, histname, 30, 0, 300)
        elif 'invariant' in varname.lower() or varname == 'invariantMass':
            hist_tuple = (histname, histname, 30, 1200, 4000)
        else:
            hist_tuple = (histname, histname, 20, 0, 1)
        hist = a.GetActiveNode().DataFrame.Histo1D(hist_tuple, varname, 'norm')
        hist.GetValue()
        out.Add(varname, hist)

    return out


# ── Entry point ───────────────────────────────────────────────────────────────
if __name__ == '__main__':
    histgroups = {}

    for setname in signal_names:
        print('Processing %s...' % setname)

        if options.select:
            histgroup = select(setname, options.year)
            outfile = ROOT.TFile.Open('rootfiles/exercise3.root', 'RECREATE')
            outfile.cd()
            histgroup.Do('Write')
            outfile.Close()
            del histgroup

        infile = ROOT.TFile.Open('rootfiles/exercise3.root')
        if infile is None:
            raise TypeError('rootfiles/exercise3.root not found - run with --select first')

        histgroups[setname] = HistGroup(setname)
        for key in infile.GetListOfKeys():
            keyname = key.GetName()
            if setname not in keyname:
                continue
            varname = keyname[len('%s_%s_' % (setname, options.year)):]
            inhist = infile.Get(keyname)
            inhist.SetDirectory(0)
            histgroups[setname].Add(varname, inhist)

    for varname in varnames.keys():
        if varname not in histgroups[signal_names[0]].keys():
            print('[skip] %s not in histgroups' % varname)
            continue
        plot_filename = plotdir + 'exercise3_{}.png'.format(varname)
        sig = signal_names[0]
        EasyPlots(
            name      = plot_filename,
            histlist  = [histgroups[sig][varname]],
            xtitle    = varnames[varname],
            ytitle    = 'Events',
            datastyle = 'hist',
        )
        print('Saved: %s' % plot_filename)
