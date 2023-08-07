from aloha_cohp.cohp_analysis import *
from aloha_cohp import plot_setting
from aloha_cohp.cohp_plotter import CohpPlot
import os



# one
# COHP = Cohpout()
# if os.path.isfile('COHPCAR.lobster'):
#     CohpPlot(COHP.pcohp(),sum_labels=self.sum_labels
#             ).plot()
# else:
#     print('COHPCAR.lobster does not exist')



# two
class Cohpplot:
    def __init__(self, sum_labels=False, 
                summed_spin=True,
                label=None, 
                orbital=None,
                sort_by=None, 
                index=None, 
                symbol=None):
        self.sum_labels = sum_labels
        self.summed_spin = summed_spin
        self.label=label
        self.orbital=orbital
        self.sort_by=sort_by
        self.index=index
        self.symbol=symbol

        self.plot_cohp()
    
    def plot_cohp(self):
        COHP = Cohpout()
        if os.path.isfile('COHPCAR.lobster'):
            CohpPlot(COHP.pcohp(summed_spin=self.summed_spin,
                                label=self.label,
                                orbital=self.orbital,
                                sort_by=self.sort_by,
                                index=self.index,
                                symbol=self.symbol,
                                )
                     ,sum_labels=self.sum_labels
                    ).plot()
        else:
            print('COHPCAR.lobster does not exist')
 


if __name__ == "__main__":
    pass
