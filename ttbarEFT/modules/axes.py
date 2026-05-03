info = {
    "njets":{
        "regular": (10, 0, 10),
        "label": r"njets",
    },
    "nleps":{
        "regular": (10, 0, 10),
        "label": r"nleps",
    },
    "ntops":{
        "regular": (10, 0, 10),
        "label": r"ntops",
    },
    "mll":{
        "regular": (40, 0, 800),
        "label": r"invariant mass of the leptons [GeV]",
    },
    "dr_leps":{
        "regular":(40, 0, 8),
        "label":r"$\Delta R$ (leading lepton, subleading lepton)",
    },
    "l0pt":{
        "regular": (40, 0, 400),
        "label":r"leading lepton $p_T$ [GeV]",
    },
    "top_pt":{
        "regular": (35, 0, 700),
        "label": r"top $p_T$ [GeV]",
    },
    "mt2":{
        "variable": [0,10,20,30,40,50,60,70,80,90,100,110,120,140,160,180,220,500],
        "label": r"$m_{T2}$ [GeV]",
    },
    "sow":{
        "regular": (1, 0, 2),
        "label": r"sun of weights",
    },
}

CR_axes = {
    "njets": {
        "regular": (10, 0, 10),
        "label": r"njets",
    },
    "nleps":{
        "regular": (10, 0, 10),
        "label": r"nleps",
    },
    "nbjets":{
        "regular": (10, 0, 10),
        "label": r"nbjets",
    },
    "l0pt":{
        "regular": (8, 0, 400),
        "label":r"leading lepton $p_T$ [GeV]",
    },
    "l0eta":{
        "regular": (30, -3, 3),
        "label":r"leading lepton $\eta$",
    },
    "l0phi":{
        "regular": (30, -3, 3),
        "label":r"leading lepton $phi$",
    },
    "l1pt":{
        "regular": (8, 0, 400),
        "label":r"subleading lepton $p_T$ [GeV]",
    },
    "l1eta":{
        "regular": (30, -3, 3),
        "label":r"subleading lepton $\eta$",
    },
    "l1phi":{
        "regular": (30, -3, 3),
        "label":r"subleading lepton $phi$",
    },
    "j0pt": {
        "regular": (10, 0, 500),
        "label": r"leading jet $p_T$ [GeV]" 
    }
}