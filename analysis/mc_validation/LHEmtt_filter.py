import numpy as np
import awkward as ak

# silence warnings due to using NanoGEN instead of full NanoAOD
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
NanoAODSchema.warn_missing_crossrefs = False

import hist
from hist import Hist
from topcoffea.modules import utils
from topcoffea.modules.histEFT import HistEFT
import topcoffea.modules.eft_helper as efth

fname = '/project01/ndcms/hnelson2/ttbarEFT/nanoGen/LHCEFT/TT01j2lCARef/nanoGen_660.root'

# Load in events from root file
events = NanoEventsFactory.from_root(
    fname,
    schemaclass=NanoAODSchema.v6,
    metadata={"dataset": "TT01j2l"},
).events()

wc_lst = utils.get_list_of_wc_names(fname)

# genpart = events.GenPart
# print("type(genpart)", type(genpart))

# test = events.LHEPart[0]
# print("test length: ", len(test))
# print("test type: ", type(test))

# for item in test:
# 	print(item)

# for i in range(10):
# 	# top = events.LHEPart[i][2]+events.LHEPart[i][3]+events.LHEPart[i][4]
# 	ttbar = events.LHEPart[i][2]+events.LHEPart[i][3]+events.LHEPart[i][4]+events.LHEPart[i][5]+events.LHEPart[i][6]+events.LHEPart[i][7]

# 	# print("particle2 : ", events.LHEPart[i][2])
# 	# print("particle3 : ", events.LHEPart[i][3])
# 	# print("particle4 : ", events.LHEPart[i][4])

# 	# print("top type: ", type(top))
# 	# print("top fields: ", top.fields)
# 	# print("top inv mass: ", top.mass)
# 	# print("top pt: ", top.pt)
# 	print("ttbar inv mass: ", ttbar.mass)
# 	print("ttbar pt: ", ttbar.pt)

# print("events.GenPart.fields", events.GenPart.fields)
# print("GenPart type: ", type(events.GenPart))
# print("events.LHEPart.fields", events.LHEPart.fields)
# print("LHEPart type: ", type(events.LHEPart))


lhepart = events.LHEPart
gluons = lhepart[abs(lhepart.pdgId) == 21]
ngluons = ak.num(gluons)

bjets = lhepart[abs(lhepart.pdgId) == 5]
nbjets = ak.num(bjets)

nlhepart = ak.num(lhepart)

h = Hist(hist.axis.Regular(bins=5, start=0, stop=5, name="ngluons", label="ngluons"))
h.fill(ngluons)

h2 = Hist(hist.axis.Regular(bins=5, start=0, stop=5, name="nbjets", label="nbjets"))
h2.fill(nbjets)

h3 = Hist(hist.axis.Regular(bins=10, start=0, stop=10, name="nlhepart", label="nlhepart"))
h3.fill(nlhepart)

# print(h3)

# print("histogram: ", h)
# print(h.values())
# print(h2)
# print(h2.values())

# event_9part = lhepart[nlhepart==9]

# for i in range(50):
	# print(len(lhepart.pdgId[i]))
	# print(event_9part[i][-1].pdgId)



