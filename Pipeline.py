import os
from Mapmaker import Mapmaker
from PowerSpectrumEstimator import PowerSpectrumEstimator
from Source.Specifications import Specifications
import multiprocessing as mp

#Pick out frequency list
#Perform Mapmaking for Every Frequency
#Perform PSE---this may want to farm out jobs

#This function makes maps at all frequencies inside the specified range
def MakeMapsAtAllFrequencies(configFile = "configuration.txt"):
    s = Specifications(os.path.dirname(os.path.abspath(__file__)), "/" + configFile)
    cubeFreqs = s.frequencyList[(s.frequencyList >= s.frequencyRange[0]) * (s.frequencyList <= s.frequencyRange[1])]
    pool = mp.Pool()
    results = [pool.apply_async(Mapmaker, (), {'freq': freq, 'useLogFile': True}) for freq in cubeFreqs]
    return [p.get() for p in results]

def RunPipeline():
    MakeMapsAtAllFrequencies()
    PowerSpectrumEstimator()

if __name__ == "__main__":    
    RunPipeline()    
    

