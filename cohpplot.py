from aloha.cohp_analysis import *
from aloha import plot_setting
from aloha.cohp_plotter import CohpPlot
import os

"""
TODO:
    1. add a function to remove duplicate values (only unique values are kept)
"""
class Cohpplot:
    def __init__(self, filepath='.',
                sum_labels=False, 
                summed_spin=True,
                label=None, 
                orbital=None,
                sort_by=None, 
                index=None, 
                symbol=None):
        self.filepath = filepath
        self.sum_labels = sum_labels
        self.summed_spin = summed_spin
        self.label=label
        self.orbital=orbital
        self.sort_by=sort_by
        self.index=index
        self.symbol=symbol

        self.plot_cohp()
    
    def plot_cohp(self):
        COHP = Cohpout(filepath=self.filepath)
        if os.path.isfile('COHPCAR.lobster'):
            CohpPlot(COHP.pcohp(summed_spin=self.summed_spin,
                                label=self.label,
                                orbital=self.orbital,
                                sort_by=self.sort_by,
                                index=self.index,
                                symbol=self.symbol,
                                )
                     ,sum_labels=self.sum_labels
                    ).plot(dpi=250)
        else:
            print('COHPCAR.lobster does not exist')
 


if __name__ == "__main__":
    pass
