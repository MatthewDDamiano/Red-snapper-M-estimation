Code for simulation testing the internal estimability of natural mortality for Southeast United States red snapper (Lutjanus campechanus)
Written by Matt Damiano and Kyle Shertzer

This repository contains two folders:
1) Base simulation framework, which has all of the code necessary to conduct a single simulation using the operating model (R) and estimation model (Beaufort Assessment Model; AD Model Builder). This can be run from the R project file "Red snapper M estimation."  
2) Profiling, which contains all of the code required to conduct single and bivariate likelihood profiling across multiple iterations of the simulation model. This can similarly be run from the R project "Profiling" in the folder.

Additional details: 
The Profiling folder contains a Results folder, which contains a copy of the base simulation framework, as the wrapper code for profiling needs to call the simulation framework from somewhere. For organization, I made individual results sub-folders for each configuration of the model.

Associated publication:
Can natural mortality be internally estimated within the stock assessment of Southeast United States Atlantic red snapper (Lutjanus campechanus)?
https://doi.org/10.1016/j.fishres.2026.107789
Reach out to Matt Damiano at matt.damiano@deq.nc.gov for a copy of the manuscript
