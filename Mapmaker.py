import numpy as np
import os
from Source.Specifications import Specifications
from Source.PrimaryBeams import PrimaryBeams
from Source.VisibilitySimulator import VisibilitySimulator

freq = 120
#TODO: the frequency should be passed as a command line parameter

scriptDirectory = os.path.dirname(os.path.abspath(__file__))
s = Specifications(scriptDirectory, "/configuration.txt",freq)
PBs = PrimaryBeams(s)

#other info we need:
	#list of LSTs
	#perVisibilityNoise in jansky
	#field specifications

if s.simulateVisibilities:
	visibilities = VisibilitySimulator(s,PBs)
else 