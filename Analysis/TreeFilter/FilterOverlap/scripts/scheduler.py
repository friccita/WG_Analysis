#!/usr/bin/env python
import os
from argparse import ArgumentParser

import scheduler_base
from scheduler_base import JobConf

p = ArgumentParser()
p.add_argument( '--run', dest='run', default=True, action='store_true', help='Run filtering' )
p.add_argument( '--check', dest='check', default=False, action='store_true', help='Run check of completion' )
p.add_argument( '--resubmit', dest='resubmit', default=False, action='store_true', help='Only submit missing output' )
p.add_argument( '--local', dest='local', default=False, action='store_true', help='Run locally' )
p.add_argument( '--step' , dest='step' , default=1    , type=int, help='Run overlap removal steps')
p.add_argument( '--year', dest='year', help='Specify the year', type=int)
p.add_argument( '--test', dest='test', default=False, action='store_true', help='Run a test job' )
options = p.parse_args()

if not options.check :
    options.run = True
else :
    options.run = False

options.batch = ( not options.local )
### ATTENTION! Here you specify the directory containing the processed ntuples, on which you want to run FilterOverlap.
#base = '/data/users/fengyb/WGammaNtuple'
#base = '/data2/users/kakw/Resonances%i/' %options.year
base = '/data/users/friccita/WGammaNtuple%i/' %options.year
### ATTENTION! Here you list the ntuple types (from RecoResonance) that you want to process, which is also the name of the subdirectory containing them.
input_dirs = [
#               'LepGamma_elg',
#               'LepLep_elel',
#               'LepGamma_mug',
#               'LepLep_mumu',
#    'LepGammaNoPhId_elg','LepGammaNoPhId_mug'
    'SingleLep_mu','SingleLep_el'
]
jobtag = '_2019_12_04'

jobs2016 = [
    # WJets HT bins, remove photon overlap with WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 (ISR+FSR),
    JobConf(base, 'WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'               , tags=['HT_pholap']      , suffix='TrueHTOlapPhOlap',year=2016),
    JobConf(base, 'WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'   , tags=['pholap']  , suffix = 'PhOlap'  ,year=2016),
    JobConf(base, 'WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'   , tags=['pholap']  , suffix = 'PhOlap'  ,year=2016),
    JobConf(base, 'WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'   , tags=['pholap']  , suffix = 'PhOlap'  ,year=2016),
    JobConf(base, 'WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'   , tags=['pholap']  , suffix = 'PhOlap'  ,year=2016),
    JobConf(base, 'WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'  , tags=['pholap']  , suffix = 'PhOlap'  ,year=2016),
    JobConf(base, 'WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8' , tags=['pholap']  , suffix = 'PhOlap'  ,year=2016),
    JobConf(base, 'WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'  , tags=['pholap']  , suffix = 'PhOlap'  ,year=2016),
    JobConf(base, 'WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'               , tags=['pholap']  , suffix = 'PhOlap'  ,year=2016),
    JobConf(base, 'WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'               , tags=['pholap'] , suffix = 'PhOlap'  ,year=2016),
    # DYJets, remove photon overlap with ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 (ISR+FSR)
    JobConf(base, 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'          , tags=['pholap']  , suffix = 'PhOlap'  ,year=2016),
    JobConf(base, 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'         , tags=['pholap'] , suffix = 'PhOlap'  ,year=2016),
    # TTJets, remove photon overlap with TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8 (ISR-only)
    JobConf(base, 'TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'            , tags=['pholapISR']  , suffix = 'PhOlap'  ,year=2016),
    JobConf(base, 'TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8', tags=['pholapISR']  , suffix = 'PhOlap'  ,year=2016),
    JobConf(base, 'TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'   , tags=['pholapISR']  , suffix = 'PhOlap'  ,year=2016),
    # WGamma photon pt bins
    JobConf(base, 'WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'                , tags=['PtMaxInc'], suffix = 'PhCutMax',year=2016),
    JobConf(base, 'WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'                 , tags=['PtMaxInc'], suffix = 'PhCutMax',year=2016),
    JobConf(base, 'WGToLNuG_PtG-130_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'        , tags=['PtRange130'], suffix = 'PhCutRange',year=2016),
    JobConf(base, 'WGToLNuG_PtG-130_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'         , tags=['PtRange130'], suffix = 'PhCutRange',year=2016),
    JobConf(base, 'WGToLNuG_PtG-500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'        , tags=['PtMin500'], suffix = 'PhCutMin',year=2016),
    JobConf(base, 'WGToLNuG_PtG-500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'         , tags=['PtMin500'], suffix = 'PhCutMin',year=2016),
]

jobs2017 = [
    # WJets HT bins, remove photon overlap with WGToLNuG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8 (ISR+FSR),
    JobConf(base, 'WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8', tags=['HT_pholap15'], suffix='TrueHTOlapPhOlap', year =2017),
    JobConf(base, 'WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8'   , tags=['pholap15']  , suffix = 'PhOlap' ,year=2017),
    JobConf(base, 'WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8'   , tags=['pholap15']  , suffix = 'PhOlap' ,year=2017),
    JobConf(base, 'WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8'   , tags=['pholap15']  , suffix = 'PhOlap' ,year=2017),
    JobConf(base, 'WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8'   , tags=['pholap15']  , suffix = 'PhOlap' ,year=2017),
    JobConf(base, 'WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8'  , tags=['pholap15']  , suffix = 'PhOlap' ,year=2017),
    JobConf(base, 'WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8' , tags=['pholap15']  , suffix = 'PhOlap' ,year=2017),
    JobConf(base, 'WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8'  , tags=['pholap15']  , suffix = 'PhOlap' ,year=2017),
    JobConf(base, 'WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8'               , tags=['pholap15']  , suffix = 'PhOlap' ,year=2017),
    JobConf(base, 'WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8'              , tags=['pholap15']  , suffix = 'PhOlap' ,year=2017),
    # DYJets, remove photon overlap with ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8 (ISR+FSR)
    JobConf(base, 'DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8'          , tags=['pholap15']  , suffix = 'PhOlap' ,year=2017),
    JobConf(base, 'DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8'         , tags=['pholap15']  , suffix = 'PhOlap' ,year=2017),
    # TTJets, remove photon overlap with TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8 (ISR-only, 10 GeV!)
    JobConf(base, 'TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8'            , tags=['pholapISR']  , suffix = 'PhOlap'  ,year=2017),
    JobConf(base, 'TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8', tags=['pholapISR']  , suffix = 'PhOlap'  ,year=2017),
    JobConf(base, 'TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8'   , tags=['pholapISR']  , suffix = 'PhOlap'  ,year=2017),
]

jobs2018 = [
    # WJets HT bins, remove photon overlap with WGToLNuG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8 (ISR+FSR),
    JobConf(base, 'WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8', tags=['HT_pholap15'], suffix='TrueHTOlapPhOlap', year =2018),
    JobConf(base, 'WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8'   , tags=['pholap15']  , suffix = 'PhOlap' ,year=2018),
    JobConf(base, 'WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8'   , tags=['pholap15']  , suffix = 'PhOlap' ,year=2018),
    JobConf(base, 'WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8'   , tags=['pholap15']  , suffix = 'PhOlap' ,year=2018),
    JobConf(base, 'WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8'   , tags=['pholap15']  , suffix = 'PhOlap' ,year=2018),
    JobConf(base, 'WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8'  , tags=['pholap15']  , suffix = 'PhOlap' ,year=2018),
    JobConf(base, 'WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8' , tags=['pholap15']  , suffix = 'PhOlap' ,year=2018),
    JobConf(base, 'WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8'  , tags=['pholap15']  , suffix = 'PhOlap' ,year=2018),
    JobConf(base, 'WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8'               , tags=['pholap15']  , suffix = 'PhOlap' ,year=2018),
    JobConf(base, 'WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8'              , tags=['pholap15']  , suffix = 'PhOlap' ,year=2018),
    # DYJets, remove photon overlap with ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8 (ISR+FSR)
    JobConf(base, 'DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8'          , tags=['pholap15']  , suffix = 'PhOlap' ,year=2018),
    JobConf(base, 'DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8'         , tags=['pholap15']  , suffix = 'PhOlap' ,year=2018),
    # TTJets overlapping with TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8 (ISR-only, 10 GeV!)
    JobConf(base, 'TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8'            , tags=['pholapISR']  , suffix = 'PhOlap'  ,year=2018),
    JobConf(base, 'TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8', tags=['pholapISR']  , suffix = 'PhOlap'  ,year=2018),
    JobConf(base, 'TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8'   , tags=['pholapISR']  , suffix = 'PhOlap'  ,year=2018),
    # Info for WWTo2L2Nu, do NOT remove photon overlap with WWG_TuneCP5_13TeV-amcatnlo-pythia8 (which has different decay: WWAToLNu2jA, photon pt>20)
]

if options.year==2016: jobs=jobs2016
if options.year==2017: jobs=jobs2017
if options.year==2018: jobs=jobs2018 

options.nFilesPerJob = 1
options.nproc = 2
options.treename='UMDNTuple/EventTree'
options.exename='RunAnalysis'
options.copyInputFiles=True
options.enableKeepFilter=False
if options.test :
    options.nproc = 1
    options.nFilesPerJob = 10
    options.totalEvents = 1001
    options.nJobs = 1
    options.batch = False
    options.local = True


module = 'Conf.py'

configs = []

for input_dir in input_dirs : 

    configs.append( {
                     'module' : module,
                     'args' : {},
                     'input' : input_dir+jobtag,
                     'output' : base + '/' + input_dir+jobtag,
                     'tag' : 'olap', 
                     'args_tag_pholap'           : { 'function' : 'filter_photon', 'pt_cut' : ' > 10 ', 'nph_cut' : ' == 0 ', 'isPromptFS_cut' : ' == True ', 'aeta_cut' : ' < 2.6 ', 'dr_cut' : ' > 0.05 ' },
                     'args_tag_pholapISR'         : { 'function' : 'filter_photon', 'pt_cut' : ' > 10 ', 'nph_cut' : ' == 0 ', 'aeta_cut' : ' < 2.6 ', 'dr_cut' : ' > 0.0 ', 'isr_cut' : ' == True ' },
                     'args_tag_PtMaxInc'         : { 'function' : 'filter_photon', 'pt_cut' : ' > 180 ', 'nph_cut' : ' == 0' },
                     'args_tag_PtRange130'         : { 'function' : 'filter_photon', 'pt_cut' : ' < 180 | > 550 ', 'nph_cut' : ' == 0' },
                     'args_tag_PtMin500'         : { 'function' : 'filter_photon', 'pt_cut' : ' > 550 ', 'nph_cut' : ' > 0 ' },
                     'args_tag_HT'               : { 'function' : 'filter_genht' , 'trueht_cut' :  ' < 100 ' },
                     'args_tag_mtMax600'         : { 'function' : 'filter_mtres' , 'mtres_cut' : ' < 600 ' },
                     'args_tag_mtMax1300'        : { 'function' : 'filter_mtres' , 'mtres_cut' : ' < 1300 ' },
                     'args_tag_mtMax500'         : { 'function' : 'filter_mtres' , 'mtres_cut' : ' < 500 ' },
                     'args_tag_mtMax400'         : { 'function' : 'filter_mtres' , 'mtres_cut' : ' < 400 ' },
                     'args_tag_mtMin400Max1200'  : { 'function' : 'filter_mtres' , 'mtres_cut' : ' >= 400 && < 1200 ' },
                     'args_tag_mtMin1200'        : { 'function' : 'filter_mtres' , 'mtres_cut' : ' >= 1200 ' },
                     'args_tag_HT_pholap'        : { 'function' : 'filter_combined', 'pt_cut' : ' > 10 ', 'nph_cut' : ' == 0 ', 'isPromptFS_cut' : ' == True ', 'aeta_cut' : ' < 2.6 ', 'dr_cut' : ' > 0.05 ', 'trueht_cut' :  ' < 100 ' },
                     # 2017/2018
                     'args_tag_HT_pholap15'        : { 'function' : 'filter_combined', 'pt_cut' : ' > 15 ', 'nph_cut' : ' == 0 ', 'isPromptFS_cut' : ' == True ', 'aeta_cut' : ' < 2.6 ', 'dr_cut' : ' > 0.05 ', 'trueht_cut' :  ' < 100 ' },
                     'args_tag_pholap15'           : { 'function' : 'filter_photon', 'pt_cut' : ' > 15 ', 'nph_cut' : ' == 0 ', 'isPromptFS_cut' : ' == True ', 'aeta_cut' : ' < 2.6 ', 'dr_cut' : ' > 0.05 ' },
    })


scheduler_base.RunJobs( jobs, configs, options)

