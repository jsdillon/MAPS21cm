# SUPPORTING CLASS FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# by Josh Dillon

import ConfigParser

#This class takes the location of the configuration file and loads relevant specs as attributes to the object
class Specifications:
    def __init__(self,configLocation):
        config = ConfigParser.ConfigParser()
        config.read(configLocation)
        try: 
            self.test1 = config.get('Input Data Settings','test1')
        except:
            print "\nWARNING: Something is wrong with configuration.txt\n"
        
