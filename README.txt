JOINT MAPMAKING AND POWER SPECTRUM ESTIMATION PIPELINE FOR 21 CM COSMOLOGY
by Josh Dillon

A work in progress.



Baselines:
Baseline vectors are written in ASCII and expressed 


Polarization Convention:
The sky polarization is expressed in terms of I, Q, U, and V. For a given normal vector n-hat, the transverse directions x and y (both of which are functions of n-hat) are given such that +x points towards the South, and +y points towards the East. That means that +Q is N/S, -Q is E/W, +U is SE/NW, and -U is NE/SW. This follows the HEALPix / Cosmology convetion: http://healpix.jpl.nasa.gov/html/intronode12.htm

Each feed is characterized by a complex beam given by the antenna index (0 indexed), the feed polarization (X or Y or something else), the sky polarization (x or y, see above), the pointing index (for drift-and-shift intruments), and the frequency.

beamFileFormat: InstrumentData/Beams/beam_[antIndex]_[antPol]_[skyPol]_[pointIndex]_[freq].npz




Possible Improvements:
-Right now, beams are linearly interpolated between frequencies of given beams. We might want to do better than that.
