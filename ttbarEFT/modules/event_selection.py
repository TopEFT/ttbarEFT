import awkward as ak

triggers_dict = {
	"2016": {
		"SingleElectron":[
			"Ele27_WPTight_Gsf",
			"Photon175",	
		],
		"SingleMuon": [
			"Mu50", 
			"TkMu50",
		],
		"DoubleEG":[
			"Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
		],
	},
	"2017": {
		"SingleElectron":[
			"Ele35_WPTight_Gsf",
			"Photon200",
		],
		"SingleMuon":[
			"Mu50",
			"TkMu100",
			"OldMu100",	
		],
		"DoubleEG":[
			"Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
		],
	},
	"2018": {
		"EGamma":[
			"Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
			"Ele32_WPTight_Gsf",
			"Photon200",
		],
		"SingleMuon":[
			"Mu50", 
			"TkMu100",
			"OldMu100",
		],
	},
}

exclude_triggers_dict = {
	"2016": {
		"SingleMuon": [],
		"SingleElectron": triggers_dict["2016"]["SingleMuon"],
		"DoubleEG": triggers_dict["2016"]["SingleMuon"] + triggers_dict["2016"]["SingleElectron"],
	},	
	"2017": {
		"SingleMuon": [],
		"SingleElectron": triggers_dict["2017"]["SingleMuon"],
		"DoubleEG": triggers_dict["2017"]["SingleMuon"] + triggers_dict["2017"]["SingleElectron"],
	},
	"2018": {
		"SingleMuon": [],
		"EGamma": triggers_dict["2018"]["SingleMuon"],
	},
}

FCNC_triggers_dict = {
	"2016":{
		"SingleMuon":[
			"IsoMu24",
			"IsoTkMu24",
		],
		"SingleElectron":[
			"Ele27_WPTight_Gsf",
			"Ele25_eta2p1_WPTight_Gsf",
		],
		"DoubleMuon":[
			"Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
			"Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
			"Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",
			"Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
			"TripleMu_12_10_5",	
		],
		"DoubleEG":[
			"Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
			"Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
			"Ele16_Ele12_Ele8_CaloIdL_TrackIdL",
		],
		"MuonEG":[
			"Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
			"Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
			"Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
			"Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
			"Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
			"Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL",
			"Mu8_DiEle12_CaloIdL_TrackIdL",	
			"DiMu9_Ele9_CaloIdL_TrackIdL",
		],	
	},
	"2017": {
		"SingleMuon":[
			"IsoMu24",
			"IsoMu27",
		],
		"SingleElectron":[
			"Ele32_WPTight_Gsf_L1DoubleEG",
			"Ele35_WPTight_Gsf",
		],
		"DoubleMuon":[
			"Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
			"Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
			"TripleMu_12_10_5",
		],
		"DoubleEG":[
			"Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
			"Ele16_Ele12_Ele8_CaloIdL_TrackIdL",
		],
		"MuonEG":[
			"Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
			"Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
			"Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
			"Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
			"Mu8_DiEle12_CaloIdL_TrackIdL",
			"Mu8_DiEle12_CaloIdL_TrackIdL_DZ",
			"DiMu9_Ele9_CaloIdL_TrackIdL_DZ",	
		],
	},
	"2018":{
		"SingleMuon":[
			"IsoMu24",
			"IsoMu27",
		],
		"EGamma":[
			"Ele32_WPTight_Gsf",
			"Ele35_WPTight_Gsf",
			"Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
			"Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
			"Ele16_Ele12_Ele8_CaloIdL_TrackIdL",
		],
		"DoubleMuon":[
			"Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
			"Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
			"TripleMu_12_10_5",	
		],
		"MuonEG":[
			"Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
			"Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
			"Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",	
			"Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
			"Mu8_DiEle12_CaloIdL_TrackIdL",
			"Mu8_DiEle12_CaloIdL_TrackIdL_DZ",
			"DiMu9_Ele9_CaloIdL_TrackIdL_DZ",
		],
	},
}

FCNC_exclude_triggers_dict = {
	"2016":{
		"SingleMuon":[],
		"DoubleMuon": FCNC_triggers_dict["2016"]["SingleMuon"],
		"SingleElectron": FCNC_triggers_dict["2016"]["SingleMuon"] + FCNC_triggers_dict["2016"]["DoubleMuon"],
		"DoubleEG": FCNC_triggers_dict["2016"]["SingleMuon"] + FCNC_triggers_dict["2016"]["DoubleMuon"] + FCNC_triggers_dict["2016"]["SingleElectron"],
		"MuonEG": FCNC_triggers_dict["2016"]["SingleMuon"] + FCNC_triggers_dict["2016"]["DoubleMuon"] + FCNC_triggers_dict["2016"]["SingleElectron"] + FCNC_triggers_dict["2016"]["DoubleEG"],
	},
	"2017": {
		"SingleMuon": [],
		"DoubleMuon": FCNC_triggers_dict["2017"]["SingleMuon"],
		"SingleElectron": FCNC_triggers_dict["2017"]["SingleMuon"] + FCNC_triggers_dict["2017"]["DoubleMuon"],
		"DoubleEG": FCNC_triggers_dict["2017"]["SingleMuon"] + FCNC_triggers_dict["2017"]["DoubleMuon"] + FCNC_triggers_dict["2017"]["SingleElectron"],
		"MuonEG":FCNC_triggers_dict["2017"]["SingleMuon"] + FCNC_triggers_dict["2017"]["DoubleMuon"] + FCNC_triggers_dict["2017"]["SingleElectron"] + FCNC_triggers_dict["2017"]["DoubleEG"],	
	},	
	"2018":{
		"SingleMuon":[],
		"DoubleMuon": FCNC_triggers_dict["2018"]["SingleMuon"],
		"EGamma": FCNC_triggers_dict["2018"]["SingleMuon"] + FCNC_triggers_dict["2018"]["DoubleMuon"], 
		"MuonEG": FCNC_triggers_dict["2018"]["SingleMuon"] + FCNC_triggers_dict["2018"]["DoubleMuon"] + FCNC_triggers_dict["2018"]["EGamma"], 
	},
}