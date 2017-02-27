CRO parameters from Laruen J. Webb and Jeremy T. First

### OpenMM

* `CNC.xml` is a hand-converted version of `CNC.{rtp,prm}` for use with OpenMM 7.0 or later `simk.openmm.app.ForceField`

### Citations

[1] Stafford AJ, Ensign DJ, and Webb LJ.
    Vibrational Stark effect spectroscopy at the interface of Ras and Rap1A bound to the Ras binding domain of RalGDS reveals an electrostatic mechanism for protein-protein interaction.
    JPC B 114:15331, 2010.
[2] Ensign DJ and Webb LJ.
    Factors determining the electrostatic fields in molecular dynamics simulations of the Ras/effector interface.
    Proteins 79:3511, 2011.

### Original correspondence:

```
---------- Forwarded message ----------
From: Jeremy T First <jeremy_first@utexas.edu>
Date: Wed, Aug 31, 2016 at 8:10 PM
Subject: Re: CN-modified cysteine parameters for AMBER?
To: John Chodera <john.chodera@choderalab.org>
Cc: "Webb, Lauren J" <lwebb@cm.utexas.edu>, "jeremy_first@utexas.edu" <jeremy_first@utmail.utexas.edu>, Julie Behr <julie.behr@choderalab.org>


Hi John and Julie,

Attached are the Amber parameters for they cyano-cysteine residue that we developed a few years ago. They are in Gromacs format, I hope that is okay.

The parameters were developed for use with the ff03 force-field, although they should play nice with any Amber force field. I've attached two papers that detail the methods of parameterization/validation. We would appreciate it if these papers were cited for the use of these parameters.

If you are working in Gromacs, the four CRO.* files need to dropped in the directory containing the force field, and the line '#include "CRO.prm" ' needs to be added to the forcefield.itp file. The cyano carbon should be named "CD", and the cyano nitrogen should be named "NE" in the structure file. The rest of the atoms are named in analogy to cysteine.

If you run into any problems with this, please don't hesitate to email me!


Very respectfully,

Jeremy Todd First
Webb Research Group
University of Texas at Austin
Jeremy_First@utexas.edu
(443) 243-1187
```
