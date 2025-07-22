 # microneut-analysis
 ## general idea
 <p>The main plot of this is to recreate Excel's GRG non-linear solver in order to analyze large amounts of microneutralization assay data.</p>

 ## data
 <p>The general idea of this is to be able to analyze data from a microneutralization assay that is run on a 96-well plate
 and then typed into a standardized Excel sheet used at Sigmovir (see linked Excel called _____). Each sample from each animal is diluted four times (20, 80, 320, 1280) and run in duplicate. One animal sample takes up each 4x2 section of the 96-well plate and the bottom right 8 wells are for the 2x2 positive and negative controls. Each control 2x2 section is averaged for one mean positive and mean negative for each plate.</p>
 
 ## analysis
 <p>Each animal sample is then taken and analyzed on its own. The duplicates are averaged to create a mean observed value for each of the four dilutions (20, 80, 320, 1280). In the expression being minimized there are four variables: IP, a1, SP, and upwards. IP is the inflection point at which the cruve changes directions (it's S shaped so like the middle point of the S). a1 is the slope of the IP. SP is the positive control and upwards is the negative control, but those two are fixed and do not change during optimization. The log base 2 of all of the dilutions is also taken.
 The expression for the predicted value being optimized is e^(a1*(log2(dilution)-IP))/(1+e^(a1*(log2(dilution)-IP)))*SP+upwards
 It then takes the sum of least squares of the predicted values subtracted from the logs of the dilutions.<\p>
 