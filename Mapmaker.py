import numpy as np
from Source.Specifications import Specifications
import os

freq = 150
#TODO: the frequency should be passed as a command line parameter

scriptDirectory = os.path.dirname(os.path.abspath(__file__))

s = Specifications(scriptDirectory + "/configuration.txt")

print s.test1


