from make_sample_json import make_sample_json
import os

jsons_2016APV = {
    'UL16APV_DY10to50':{
        'xsec': 18610.0,
        'year': '2016APV',
        'histAxisName': 'DY10to50_centralUL16APV',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODAPVv9/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/',
        },
    'UL16APV_DY50':{
        'xsec': 6025.2,
        'year': '2016APV',
        'histAxisName': 'DY50_centralUL16APV',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODAPVv9/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/',
        },
    'UL16APV_TTGJets':{
        'xsec': 3.697,
        'year': '2016APV',
        'histAxisName': 'TTGJets_centralUL16APV',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODAPVv9/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/NANOAODSIM/',
        },
    'UL16APV_TTTo2L2Nu':{
        'xsec': 87.31483776,
        'year': '2016APV',
        'histAxisName': 'TTTo2L2Nu_centralUL16APV',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODAPVv9/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/',
        },
    'UL16APV_TTWJetsToLNu':{
        'xsec': 0.2043,
        'year': '2016APV',
        'histAxisName': 'TTWJetsToLNu_centralUL16APV',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODAPVv9/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/NANOAODSIM/',
        },
    'UL16APV_TTZToLLNuNu_M_10':{
        'xsec': 0.281,
        'year': '2016APV',
        'histAxisName': 'TTZToLLNuNu_centralUL16APV',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODAPVv9/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/NANOAODSIM/',
        },
    'UL16APV_TW_antitop_5f_NoFullyHadronicDecays':{
        'xsec': 19.47,
        'year': '2016APV',
        'histAxisName': 'TW_NoFullyHadronicDecays_centralUL16APV',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODAPVv9/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/',
        },
    'UL16APV_TW_top_5f_NoFullyHadronicDecays':{
        'xsec': 19.47,
        'year': '2016APV',
        'histAxisName': 'TW_NoFullyHadronicDecays_centralUL16APV',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODAPVv9/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/',
        },
    'UL16APV_WJetsToLNu':{
        'xsec': 61526.7,
        'year': '2016APV',
        'histAxisName': 'WJetsToLNu_centralUL16APV',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODAPVv9/WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/',
        },
    'UL16APV_WWTo2L2Nu':{
        'xsec': 12.178,
        'year': '2016APV',
        'histAxisName': 'WWTo2L2Nu_centralUL16APV',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODAPVv9/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/',
        },
    'UL16APV_WWW_4F':{
        'xsec': 0.2086,
        'year': '2016APV',
        'histAxisName': 'WWW_centralUL16APV',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODAPVv9/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/NANOAODSIM/',
        },
    'UL16APV_WWZ_4F':{
        'xsec': 0.1651,
        'year': '2016APV',
        'histAxisName': 'WWZ_4F_centralUL16APV',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODAPVv9/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/NANOAODSIM/',
        },
    'UL16APV_WZTo3LNu':{
        'xsec': 5.284311882352942,
        'year': '2016APV',
        'histAxisName': 'WZTo3LNu_centralUL16APV',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODAPVv9/WZTo3LNu_mllmin4p0_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/',
        },
    'UL16APV_WZZ':{
        'xsec': 0.05565,
        'year': '2016APV',
        'histAxisName': 'WZZ_centralUL16APV',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODAPVv9/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/NANOAODSIM/',
        },
    'UL16APV_ZZTo4L':{
        'xsec': 1.256,
        'year': '2016APV',
        'histAxisName': 'ZZTo4L_centralUL16APV',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODAPVv9/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/NANOAODSIM/106X_mcRun2_asymptotic_preVFP_v11-v1/',
        },
    'UL16APV_ZZZ':{
        'xsec':0.01398,
        'year': '2016APV',
        'histAxisName': 'ZZZ_centralUL16APV',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODAPVv9/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/NANOAODSIM/',
        },
}

jsons_2016 = {
    'UL16_DY10to50':{
        'xsec': 18610.0,
        'year': '2016',
        'histAxisName': 'DY10to50_centralUL16',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODv9/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/',
        },
    'UL16_DY50':{
        'xsec': 6025.2,
        'year': '2016',
        'histAxisName': 'DY50_centralUL16',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODv9/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/',
        },
    'UL16_TTGJets':{
        'xsec': 3.697,
        'year': '2016',
        'histAxisName': 'TTGJets_centralUL16',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODv9/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/NANOAODSIM/',
        },
    'UL16_TTTo2L2Nu':{
        'xsec': 87.31483776,
        'year': '2016',
        'histAxisName': 'TTTo2L2Nu_centralUL16',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODv9/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/',
        },
    'UL16_TTWJetsToLNu':{
        'xsec': 0.2043,
        'year': '2016',
        'histAxisName': 'ttW_centralUL16',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODv9/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/NANOAODSIM/',
        },
    'UL16_TTZToLLNuNu_M_10':{
        'xsec': 0.281,
        'year': '2016',
        'histAxisName': 'ttZ_centralUL16',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODv9/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/NANOAODSIM/',
        },
    'UL16_TW_antitop_5f_NoFullyHadronicDecays':{
        'xsec': 19.47,
        'year': '2016',
        'histAxisName': 'TW_NoFullyHadronicDecays_centralUL16',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODv9/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/',
        },
    'UL16_TW_top_5f_NoFullyHadronicDecays':{
        'xsec': 19.47,
        'year': '2016',
        'histAxisName': 'TW_NoFullyHadronicDecays_centralUL16',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODv9/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/',
        },
    'UL16_WJetsToLNu':{
        'xsec': 61526.7,
        'year': '2016',
        'histAxisName': 'WJetsToLNu_centralUL16',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODv9/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/',
        },
    'UL16_WWTo2L2Nu':{
        'xsec': 12.178,
        'year': '2016',
        'histAxisName': 'WWTo2L2Nu_centralUL16',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODv9/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/',
        },
    'UL16_WWW_4F':{
        'xsec': 0.2086,
        'year': '2016',
        'histAxisName': 'WWW_centralUL16',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODv9/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/NANOAODSIM/',
        },
    'UL16_WWZ_4F':{
        'xsec': 0.1651,
        'year': '2016',
        'histAxisName': 'WWZ_4F_centralUL16',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODv9/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/NANOAODSIM/',
        },
    'UL16_WZTo3LNu':{
        'xsec': 5.284311882352942,
        'year': '2016',
        'histAxisName': 'WZTo3LNu_centralUL16',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODv9/WZTo3LNu_mllmin4p0_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/',
        },
    'UL16_WZZ':{
        'xsec': 0.05565,
        'year': '2016',
        'histAxisName': 'WZZ_centralUL16',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODv9/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/NANOAODSIM/',
        },
    'UL16_ZZTo4L':{
        'xsec': 1.256,
        'year': '2016',
        'histAxisName': 'ZZTo4L_centralUL16',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODv9/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/NANOAODSIM/',
        },
    'UL16_ZZZ':{
        'xsec':0.01398,
        'year': '2016',
        'histAxisName': 'ZZZ_centralUL16',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL16NanoAODAPVv9/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/NANOAODSIM/',
        },
}

jsons_2017 = {
    'UL17_DY10to50':{
        'xsec': 18610.0,
        'year': '2017',
        'histAxisName': 'DY10to50_centralUL17',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL17NanoAODv9/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/',
        },
    'UL17_DY50':{
        'xsec': 6025.2,
        'year': '2017',
        'histAxisName': 'DY50_centralUL17',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL17NanoAODv9/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/',
        },
    'UL17_TTGJets':{
        'xsec': 3.697,
        'year': '2017',
        'histAxisName': 'TTGJets_centralUL17',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL17NanoAODv9/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v1/',
        },
    'UL17_TTTo2L2Nu':{
        'xsec': 87.31483776,
        'year': '2017',
        'histAxisName': 'TTTo2L2Nu_centralUL17',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL17NanoAODv9/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/',
        },
    'UL17_TTWJetsToLNu':{
        'xsec': 0.2043,
        'year': '2017',
        'histAxisName': 'ttW_centralUL17',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL17NanoAODv9/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v1/',
        },
    'UL17_TTZToLLNuNu_M_10':{
        'xsec': 0.281,
        'year': '2017',
        'histAxisName': 'ttZ_centralUL17',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL17NanoAODv9/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v1/',
        },
    'UL17_TW_antitop_5f_NoFullyHadronicDecays':{
        'xsec': 19.47,
        'year': '2017',
        'histAxisName': 'TW_NoFullyHadronicDecays_centralUL17',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL17NanoAODv9/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v1/',
        },
    'UL17_TW_top_5f_NoFullyHadronicDecays':{
        'xsec': 19.47,
        'year': '2017',
        'histAxisName': 'TW_NoFullyHadronicDecays_centralUL17',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL17NanoAODv9/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v1/',
        },
    'UL17_WJetsToLNu':{
        'xsec': 61526.7,
        'year': '2017',
        'histAxisName': 'WJetsToLNu_centralUL17',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL17NanoAODv9/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v1/',
        },
    'UL17_WWTo2L2Nu':{
        'xsec': 12.178,
        'year': '2017',
        'histAxisName': 'WWTo2L2Nu_centralUL17',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL17NanoAODv9/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v2/',
        },
    'UL17_WWW_4F':{
        'xsec': 0.2086,
        'year': '2017',
        'histAxisName': 'WWW_centralUL17',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL17NanoAODv9/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/NANOAODSIM/',
        },
    'UL17_WWZ_4F':{
        'xsec': 0.1651,
        'year': '2017',
        'histAxisName': 'WWZ_4F_centralUL17',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL17NanoAODv9/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/NANOAODSIM/',
        },
    'UL17_WZTo3LNu':{
        'xsec': 5.284311882352942,
        'year': '2017',
        'histAxisName': 'WZTo3LNu_centralUL17',
        'path': '/cms/cephfs/data//store/mc/RunIISummer20UL17NanoAODv9/WZTo3LNu_mllmin4p0_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/',
        },
    'UL17_WZZ':{
        'xsec': 0.05565,
        'year': '2017',
        'histAxisName': 'WZZ_centralUL17',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL17NanoAODv9/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/NANOAODSIM/',
        },
    'UL17_ZZTo4L':{
        'xsec': 1.256,
        'year': '2017',
        'histAxisName': 'ZZTo4L_centralUL17',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL17NanoAODv9/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v2/',
        },
    'UL17_ZZZ':{
        'xsec':0.01398,
        'year': '2017',
        'histAxisName': 'ZZZ_centralUL17',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL17NanoAODv9/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/NANOAODSIM/',
        },
}

jsons_2018 = {
    'UL18_DY10to50':{
        'xsec': 18610.0,
        'year': '2018',
        'histAxisName': 'DY10to50_centralUL18',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL18NanoAODv9/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/',
        },
    'UL18_DY50':{
        'xsec': 6025.2,
        'year': '2018',
        'histAxisName': 'DY50_centralUL18',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL18NanoAODv9/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/',
        },
    'UL18_TTGJets':{
        'xsec': 3.697,
        'year': '2018',
        'histAxisName': 'TTGJets_centralUL18',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL18NanoAODv9/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/NANOAODSIM/',
        },
    'UL18_TTTo2L2Nu':{
        'xsec': 87.31483776,
        'year': '2018',
        'histAxisName': 'TTTo2L2Nu_centralUL18',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL18NanoAODv9/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/',
        },
    'UL18_TTWJetsToLNu':{
        'xsec': 0.2043,
        'year': '2018',
        'histAxisName': 'ttW_centralUL18',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL18NanoAODv9/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/',
        },
    'UL18_TTZToLLNuNu_M_10':{
        'xsec': 0.281,
        'year': '2018',
        'histAxisName': 'ttZ_centralUL18',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL18NanoAODv9/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/',
        },
    'UL18_TW_antitop_5f_NoFullyHadronicDecays':{
        'xsec': 19.47,
        'year': '2018',
        'histAxisName': 'TW_NoFullyHadronicDecays_centralUL18',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL18NanoAODv9/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/',
        },
    'UL18_TW_top_5f_NoFullyHadronicDecays':{
        'xsec': 19.47,
        'year': '2018',
        'histAxisName': 'TW_NoFullyHadronicDecays_centralUL18',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL18NanoAODv9/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/',
        },
    'UL18_WJetsToLNu':{
        'xsec': 61526.7,
        'year': '2018',
        'histAxisName': 'WJetsToLNu_centralUL18',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL18NanoAODv9/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/',
        },
    'UL18_WWTo2L2Nu':{
        'xsec': 12.178,
        'year': '2018',
        'histAxisName': 'WWTo2L2Nu_centralUL18',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL18NanoAODv9/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/',
        },
    'UL18_WWW_4F':{
        'xsec': 0.2086,
        'year': '2018',
        'histAxisName': 'WWW_centralUL18',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL18NanoAODv9/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/NANOAODSIM/',
        },
    'UL18_WWZ_4F':{
        'xsec': 0.1651,
        'year': '2018',
        'histAxisName': 'WWZ_4F_centralUL18',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL18NanoAODv9/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/NANOAODSIM/',
        },
    'UL18_WZTo3LNu':{
        'xsec': 5.284311882352942,
        'year': '2018',
        'histAxisName': 'WZTo3LNu_centralUL18',
        'path': '/cms/cephfs/data//store/mc/RunIISummer20UL18NanoAODv9/WZTo3LNu_mllmin4p0_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/',
        },
    'UL18_WZZ':{
        'xsec': 0.05565,
        'year': '2018',
        'histAxisName': 'WZZ_centralUL18',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL18NanoAODv9/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/NANOAODSIM/',
        },
    'UL18_ZZTo4L':{
        'xsec': 1.256,
        'year': '2018',
        'histAxisName': 'ZZTo4L_centralUL18',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL18NanoAODv9/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/',
        },
    'UL18_ZZZ':{
        'xsec':0.01398,
        'year': '2018',
        'histAxisName': 'ZZZ_centralUL18',
        'path': '/cms/cephfs/data/store/mc/RunIISummer20UL18NanoAODv9/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/NANOAODSIM/',
        },
}

new_SMEFTsim_ttbar = {
    'UL16APV_TT01j2l_mtt0to700': {
        "xsec": 77.0803,
        "year": "2016APV",
        "histAxisName": "TT01j2l_UL16APV_mtt_0to700",
        "paths": ["/cms/cephfs/data/store/user/hnelson2/mc/ttbarEFT_Run2/UL16APV/postLHE/v1/NAOD_TTto2L2Nu_1Jets_smeft_MTT_0to700", 
                "/cms/cephfs/data/store/user/hnelson2/mc/ttbarEFT_Run2/UL16APV/postSIM/v1/NAOD_TTto2L2Nu_1Jets_smeft_MTT_0to700/"]
        },
    'UL16APV_TT01j2l_mtt700to900': {
        "xsec": 6.9496,
        "year": "2016APV",
        "histAxisName": "TT01j2l_UL16APV_mtt_700to900",
        "paths": ["/cms/cephfs/data/store/user/hnelson2/mc/ttbarEFT_Run2/UL16APV/postLHE/v1/NAOD_TTto2L2Nu_1Jets_smeft_MTT_700to900",],
        },  
    'UL16APV_TT01j2l_mtt900toInf': {
        "xsec": 3.5296,
        "year": "2016APV",
        "histAxisName": "TT01j2l_UL16APV_mtt_900toInf",
        "paths": ["/cms/cephfs/data/store/user/hnelson2/mc/ttbarEFT_Run2/UL16APV/postLHE/v1/NAOD_TTto2L2Nu_1Jets_smeft_MTT_900toInf"],
        },
    'UL16_TT01j2l_mtt0to700': {
        "xsec": 77.0803,
        "year": "2016",
        "histAxisName": "TT01j2l_UL16_mtt_0to700",
        "paths": ["/cms/cephfs/data/store/user/hnelson2/mc/ttbarEFT_Run2/UL16/postSIM/v1/NAOD_TTto2L2Nu_1Jets_smeft_MTT_0to700"],
        },
    'UL16_TT01j2l_mtt700to900': {
        "xsec": 6.9496,
        "year": "2016",
        "histAxisName": "TT01j2l_UL16_mtt_700to900",
        "paths": ["/cms/cephfs/data/store/user/hnelson2/mc/ttbarEFT_Run2/UL16/postSIM/v1/NAOD_TTto2L2Nu_1Jets_smeft_MTT_700to900"],
        },
    'UL16_TT01j2l_mtt900toInf': {
        "xsec": 3.5296,
        "year": "2016",
        "histAxisName": "TT01j2l_UL16_mtt_900toInf",
        "paths": ["/cms/cephfs/data/store/user/hnelson2/mc/ttbarEFT_Run2/UL16/postSIM/v1/NAOD_TTto2L2Nu_1Jets_smeft_MTT_900toInf"],
        },
    'UL17_TT01j2l_mtt0to700': {
        "xsec": 77.0803,
        "year": "2017",
        "histAxisName": "TT01j2l_UL17_mtt_0to700",
        "paths": ["/cms/cephfs/data/store/user/hnelson2/mc/ttbarEFT_Run2/UL17/postLHE/v3/NAOD_TTto2L2Nu_1Jets_smeft_MTT_0to700_v1",
            "/cms/cephfs/data/store/user/hnelson2/mc/ttbarEFT_Run2/UL17/postLHE/v3/NAOD_TTto2L2Nu_1Jets_smeft_MTT_0to700_v2"],
        },
    'UL17_TT01j2l_mtt700to900': {
        "xsec": 6.9496,
        "year": "2017",
        "histAxisName": "TT01j2l_UL17_mtt_700to900",
        "paths": ["/cms/cephfs/data/store/user/hnelson2/mc/ttbarEFT_Run2/UL17/postLHE/v3/NAOD_TTto2L2Nu_1Jets_smeft_MTT_700to900_v1",
            "/cms/cephfs/data/store/user/hnelson2/mc/ttbarEFT_Run2/UL17/postLHE/v3/NAOD_TTto2L2Nu_1Jets_smeft_MTT_700to900_v2"],
        },
    'UL17_TT01j2l_mtt900toInf': {
        "xsec": 3.5296,
        "year": "2017",
        "histAxisName": "TT01j2l_UL17_mtt_900toInf",
        "paths": ["/cms/cephfs/data/store/user/hnelson2/mc/ttbarEFT_Run2/UL17/postLHE/v3/NAOD_TTto2L2Nu_1Jets_smeft_MTT_900toInf_v1",
            "/cms/cephfs/data/store/user/hnelson2/mc/ttbarEFT_Run2/UL17/postLHE/v3/NAOD_TTto2L2Nu_1Jets_smeft_MTT_900toInf_v2"],
        },
    'UL18_TT01j2l_mtt0to700': {
        "xsec": 77.0803,
        "year": "2017",
        "histAxisName": "TT01j2l_UL17_mtt_0to700",
        "paths": ["/cms/cephfs/data/store/user/hnelson2/mc/ttbarEFT_Run2/UL18/postSIM/v4/NAOD_TTto2L2Nu_1Jets_smeft_MTT_0to700"],
        },
    'UL18_TT01j2l_mtt700to900': {
        "xsec": 6.9496,
        "year": "2017",
        "histAxisName": "TT01j2l_UL18_mtt_700to900",
        "paths": ["/cms/cephfs/data/store/user/hnelson2/mc/ttbarEFT_Run2/UL18/postSIM/v2/NAOD_TTto2L2Nu_1Jets_smeft_MTT_700to900"],
        },
    'UL18_TT01j2l_mtt900toInf': {
        "xsec": 3.5296,
        "year": "2017",
        "histAxisName": "TT01j2l_UL18_mtt_900toInf",
        "paths": ["/cms/cephfs/data/store/user/hnelson2/mc/ttbarEFT_Run2/UL18/postSIM/v3/NAOD_TTto2L2Nu_1Jets_smeft_MTT_900toInf"],
        },
    }


# for json, dict in jsons_2016APV.items():
#     make_sample_json(
#         xsec=dict['xsec'], 
#         year=dict['year'], 
#         treeName='Events', 
#         histAxisName=dict['histAxisName'], 
#         options='',
#         era='',
#         path=dict['path'], 
#         outname=json,
#     )

# for json, dict in jsons_2016.items():
#     make_sample_json(
#         xsec=dict['xsec'], 
#         year='2016', 
#         treeName='Events', 
#         histAxisName=dict['histAxisName'], 
#         options='',
#         era='',
#         path=dict['path'], 
#         outname=json,
#     )

# for json, dict in jsons_2017.items():
#     make_sample_json(
#         xsec=dict['xsec'], 
#         year='2017', 
#         treeName='Events', 
#         histAxisName=dict['histAxisName'], 
#         options='',
#         era='',
#         path=dict['path'], 
#         outname=json,
#     )

# for json, dict in jsons_2018.items():
#     make_sample_json(
#         xsec=dict['xsec'], 
#         year='2018', 
#         treeName='Events', 
#         histAxisName=dict['histAxisName'], 
#         options='',
#         era='',
#         path=dict['path'], 
#         outname=json,
#     )

for json, dict in new_SMEFTsim_ttbar.items():
    make_sample_json(
        xsec=dict['xsec'], 
        year=dict['year'], 
        treeName='Events', 
        histAxisName=dict['histAxisName'], 
        options='',
        era='',
        paths=dict['paths'], 
        outname=json,
    )