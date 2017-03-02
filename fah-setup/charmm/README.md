Starting with 1OL5-WT-pdbfixer.pdb
* (need to fill in loop of TPX2)

* Manually (VI or emacs) change all water atoms from ATOM to HETATM
* Load into charmm-gui
* pdb type: RCSB
* Next:
* select all chains (+TPX2)
* rename ADP to ADP
* rename MG to MG
* should not have to rename waters (will be renamed to TIP)
* preserve hydrogen coordinates

* repeat without chain B for -TPX2

will need to define box size (?)
put all that into parmed

parmed.github.io omm_charmm (might be old)

for coord in coords: find min and max to find box, then setBox
sometimes in the pb it'll be at the top
visualize to make sure box is in the right position

github.com/chayast/charmm36_ff
par_all36_prot.prm and top (rtf?)
stream/na/nad_ppi.str -- chaya is emailing
using this file, also need parent file (all linked)

should be able to follow example

parmed should be able to set up water properly, was "officially" fixed, might still be problematic
(issue was reopened, status unknown)
not exactly the same as openmm
might not need .str for waters, but just in case it's linked (parmed should do internally)

not CharmmCrdFile -- use openmm PDBFile (use to get positions)
setPositions(app.PDBFile.positions)


