set defer_build_mode, 1
load trajectory-aligned.pdb
load_traj trajectory-aligned.dcd
hide all
select AurA, resi 123-388
deselect
intra_fit trajectory-aligned and AurA
show cartoon, all
# AurA
color green, AurA
# T288
show sticks, resi 288
# activation loop
color red, resi 281-293
distance resi 284 and name CB, resi 225 and name CB
show sticks, resi 284 or resi 225
# Tpx2
color yellow, resi 1-41
# ADP
show sticks, resn MOL
util.cbay('resn MOL')
