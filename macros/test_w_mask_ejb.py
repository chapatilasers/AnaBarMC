import uproot
import particle
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt

# Open the G4SBS root file and navigate to the tree containing the event data
file = uproot.open("gep_12Gev1000.root")
tree = file["T"]

# Get the arrays for the variables of interest
cdet_hit = tree["Earm.CDET_Scint.hit.nhits"].array()
track_indexes = tree["Earm.CDET_Scint.hit.sdtridx"].array()
pid = tree["SDTrack.PID"].array()
xpos = tree["SDTrack.posx"].array()
ypos = tree["SDTrack.posy"].array()
zpos = tree["SDTrack.posz"].array()
xmomentum = tree["SDTrack.momx"].array()
ymomentum = tree["SDTrack.momy"].array()
zmomentum = tree["SDTrack.momz"].array()
energy = tree["SDTrack.Etot"].array()


dat_for_tree = {"X_vtx": [], "Y_vtx": [], "Z_vtx": [], "Px_p": [], "Py_p": [], "Pz_p": [], "En_p": [], "Mass": [], "PDG": []}
for event in range(len(cdet_hit)):

    for hit in range(cdet_hit[event]):

        track_idx = track_indexes[event][hit]

        dat_for_tree["X_vtx"].append(xpos[event][track_idx])
        dat_for_tree["Y_vtx"].append(ypos[event][track_idx])
        dat_for_tree["Z_vtx"].append(zpos[event][track_idx])
        dat_for_tree["Px_p"].append(xmomentum[event][track_idx])
        dat_for_tree["Py_p"].append(ymomentum[event][track_idx])
        dat_for_tree["Pz_p"].append(zmomentum[event][track_idx])
        dat_for_tree["En_p"].append(energy[event][track_idx])
        dat_for_tree["Mass"].append(particle.Particle.from_pdgid(pid[event][track_idx]).mass)
        dat_for_tree["PDG"].append(pid[event][track_idx])

# Create a new ROOT file with the extracted information
#file = uproot.recreate("del.root")
#file["h1"] = dat_for_tree
#file.close()

fig = plt.figure()
gs = fig.add_gridspec(2, 3)
axs = gs.subplots()
plts = axs.flat

print("\n\n##########OLD POSITIONS##########\n")

plts[0].hist(dat_for_tree["X_vtx"], bins=50)
print(np.mean(dat_for_tree["X_vtx"]))
plts[0].set_xlabel('X_vtx')
plts[0].set_ylabel('counts')
plts[0].set_title('X_vtx')


plts[1].hist(dat_for_tree["Y_vtx"], bins=50)
print(np.mean(dat_for_tree["Y_vtx"]))
plts[1].set_xlabel('Y_vtx')
plts[1].set_ylabel('counts')
plts[1].set_title('Y_vtx')


plts[2].hist(dat_for_tree["Z_vtx"], bins=50)
print(np.mean(dat_for_tree["Z_vtx"]))
plts[2].set_xlabel('Z_vtx')
plts[2].set_ylabel('counts')
plts[2].set_title('Z_vtx')


print("\n\n##########OLD Momenta##########\n")

plts[3].hist(dat_for_tree["Px_p"], bins=50)
print(np.mean(dat_for_tree["Px_p"]))
plts[3].set_xlabel('Px_p')
plts[3].set_ylabel('counts')
plts[3].set_title('Px_p')


plts[4].hist(dat_for_tree["Py_p"], bins=50)
print(np.mean(dat_for_tree["Py_p"]))
plts[4].set_xlabel('Py_p')
plts[4].set_ylabel('counts')
plts[4].set_title('Py_p')


plts[5].hist(dat_for_tree["Pz_p"], bins=50)
print(np.mean(dat_for_tree["Pz_p"]))
plts[5].set_xlabel('Pz_p')
plts[5].set_ylabel('counts')
plts[5].set_title('Pz_p')


fig2 = plt.figure()
gs = fig2.add_gridspec(2, 3)
axs = gs.subplots()
plts = axs.flat

tx = np.array(dat_for_tree["Z_vtx"])
ty = np.array(dat_for_tree["X_vtx"])
angle = np.radians(27.0)

print("angle:", angle)

#newx = -(tx * np.cos(angle) + ty * np.sin(angle)) * 100
#newy = -((-tx * np.sin(angle) + ty * np.cos(angle)) - 4.0735) * 100
newx = -(-tx * np.sin(angle) + ty * np.cos(angle)) * 100
newy = -((tx * np.cos(angle) + ty * np.sin(angle)) - 4.0735) * 100
newz = -np.array(dat_for_tree["Y_vtx"]) * 100

print("\n\n##########New POSITIONS##########\n")

plts[0].hist(newx, bins=50)
print(np.mean(newx))
plts[0].set_xlabel('newx')
plts[0].set_ylabel('counts')
plts[0].set_title('newx')


plts[1].hist(newy, bins=50)
print(np.mean(newy))
plts[1].set_xlabel('newy')
plts[1].set_ylabel('counts')
plts[1].set_title('newy')


plts[2].hist(newz, bins=50)
print(np.mean(newz))
plts[2].set_xlabel('newz')
plts[2].set_ylabel('counts')
plts[2].set_title('newz')



txp = np.array(dat_for_tree["Pz_p"])
typ = np.array(dat_for_tree["Px_p"])

newxp = -(-txp * np.sin(angle) + typ * np.cos(angle))
newyp = -(txp * np.cos(angle) + typ * np.sin(angle))
newzp = -np.array(dat_for_tree["Py_p"])

print("\n\n##########New Momenta##########\n")

plts[3].hist(newxp, bins=50)
print(np.mean(newxp))
plts[3].set_xlabel('Px_p')
plts[3].set_ylabel('counts')
plts[3].set_title('Px_p')


plts[4].hist(newyp, bins=50)
print(np.mean(newyp))
plts[4].set_xlabel('Py_p')
plts[4].set_ylabel('counts')
plts[4].set_title('Py_p')


plts[5].hist(newzp, bins=50)
print(np.mean(newzp))
plts[5].set_xlabel('Pz_p')
plts[5].set_ylabel('counts')
plts[5].set_title('Pz_p')

plt.show()





