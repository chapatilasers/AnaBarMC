import uproot
import awkward
import matplotlib.pyplot as plt
from scipy.stats import norm
import scipy

with uproot.open("gep_12GeV2.root") as f:
	tree = f["T"]
	events = tree.arrays()
'''
	sdtrid = tree["Earm.CDET_Scint.hit.sdtridx"].array()
	testz = tree["Earm.CDET_Scint.hit.zhitg"].array()
	zvert = tree["SDTrack.vz"].array()
	PID = tree["SDTrack.PID"].array()
	TID = tree["SDTrack.TID"].array()
	posx = tree["SDTrack.posx"].array()
	posy = tree["SDTrack.posy"].array()
	posz = tree["SDTrack.posz"].array()
	momx = tree["SDTrack.momx"].array()
	momy = tree["SDTrack.momy"].array()
	momz = tree["SDTrack.momz"].array()
	etot = tree["SDTrack.Etot"].array()
'''

dat_for_tree = {"X_vtx":[], "Y_vtx":[], "Z_vtx":[], "Px_p":[], "Py_p":[], "Pz_p":[], "En_p":[], "PDG":[]}


for event in events:
	cDetHits = []
	for jHit in range(len(event["Earm.CDET_Scint.hit.xhit"])):
		xHit = event["Earm.CDET_Scint.hit.xhit"][jHit]
		print(xHit)
		yHit = event["Earm.CDET_Scint.hit.yhit"][jHit]
		zHit = event["Earm.CDET_Scint.hit.zhit"][jHit]
		if xHit != 0 and yHit != 0 and zHit != 0:
			cDetHits.append(jHit)
	if len(cDetHits) > 0:
		for jhit in cDetHits:
			print("jhit:", jhit)
			print("length of event['SDTrack.posy']:", len(event["SDTrack.posy"]))
		# need to translate into the coordinate frame (and units) of the CDet simulation
			dat_for_tree["X_vtx"].append(-event["SDTrack.posy"][jhit] * 100)
			dat_for_tree["Y_vtx"].append(-event["SDTrack.posz"][jhit] * 100)
			dat_for_tree["Z_vtx"].append(event["SDTrack.posx"][jhit] * 100)
			dat_for_tree["Px_p"].append(-event["SDTrack.momy"][jhit] * 1000)
			dat_for_tree["Py_p"].append(-event["SDTrack.momz"][jhit] * 1000)
			dat_for_tree["Pz_p"].append(event["SDTrack.momx"][jhit] * 1000)
			dat_for_tree["En_p"].append(event["SDTrack.Etot"][jhit] * 1000)
			dat_for_tree["PDG"].append(event["SDTrack.PID"][jhit])


file = uproot.recreate("proccessed_demonstration.root")
file["h1"] = dat_for_tree
file.close()


