from mendeleev import element
import sys
import os
from pymatgen.electronic_structure.cohp import CompleteCohp, get_integrated_cohp_in_energy_range
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.plotter import CohpPlotter
import itertools
from ase.io import read
import pandas as pd
import numpy as np
from matplotlib import rcParams
import matplotlib.pyplot as plt
from aloha import plot_setting

pd.set_option('display.expand_frame_repr', False)
pd.set_option('display.max_columns', None)
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_rows', None)

orbital_order={'1s':0,'2s':1,'2p':2,'3s':3,'3p':4,'4s':5,'3d':6,'4p':7,'5s':8,'4d':9,'5p':10,'6s':11,'4f':12,
               '5d':13,'6p':14,'7s':15,'5f':16,'6d':17,'7p':18}
metal_elements = ['Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi']

def rm_num(string):
    new=''.join([i for i in string if not i.isdigit()])
    return new

class Cohpout:
    def __init__(self,
                 filepath='.',
                ):
        COHPCAR_path = filepath+"/COHPCAR.lobster"
        POSCAR_path = filepath+"/POSCAR"
        
        self.cohpcar=COHPCAR_path
        self.poscar=POSCAR_path
        self.completecohp=CompleteCohp.from_file(fmt="LOBSTER", filename=COHPCAR_path, structure_file=POSCAR_path)
        self.atoms=read(POSCAR_path)
        self.d, self.dataframe =self._dict()

    def _dict(self):
        """ create dataframe for COHP """
        data=[]
        d={}
        for i in range(len(self.completecohp.bonds)):
            label = str(i+1)
            row = {}

            """ finding element index based on POSCAR (therefore, started from 0) """
            elem1 = str(self.completecohp.bonds[label]['sites'][0]._species)[:-1]
            elem2 = str(self.completecohp.bonds[label]['sites'][1]._species)[:-1]
            x1, y1, z1 = self.completecohp.bonds[label]['sites'][0].coords
            elem1_idx = [atom.index for atom in self.atoms if x1-0.01<atom.x<x1+0.01 and y1-0.01<atom.y<y1+0.01 and z1-0.01<atom.z<z1+0.01]
            x2, y2, z2 = self.completecohp.bonds[label]['sites'][1].coords
            elem2_idx = [atom.index for atom in self.atoms if x2-0.01<atom.x<x2+0.01 and y2-0.01<atom.y<y2+0.01 and z2-0.01<atom.z<z2+0.01]

            row['label'] = label
            row['pair']=f"{elem1}({elem1_idx[0]})-{elem2}({elem2_idx[0]})"
            row['elem1']=elem1
            row['elem2']=elem2
            row['elem1_idx']=elem1_idx[0]
            row['elem2_idx']=elem2_idx[0]
            row['length']=float(self.completecohp.bonds[label]['length'])
            row['icohp']=-float(get_integrated_cohp_in_energy_range(self.completecohp,label))
            elements=[elem1,elem2]
            row['elements']=elements
            for elem in elements:
                ele_obj=element(elem).ec.conf
                ele_orb=list(ele_obj.keys())
                row[elem]={'tot_orb':ele_orb}
                for i in range(len(ele_orb)):
                    n = int(ele_orb[i][0])
                    n_orb = ele_orb[i][1]
                    if n_orb == 's':
                        row[elem]['s_orb']=[(n,"s")]
                    elif n_orb == 'p':
                        row[elem]['p_orb']=[(n,"px"),(n,"py"),(n,"pz")]
                    elif n_orb == 'd':
                        row[elem]['d_orb']=[(n,"dxy"),(n,"dyz"),(n,"dz2"),(n,"dxz"),(n,"dx2")]
                    elif n_orb == 'f':
                        row[elem]['f_orb']=[(n,"f_3"),(n,"f_2"),(n,"f_1"),(n,"f0"),(n,"f1"),(n,"f2"),(n,"f3")]
                row[elem]['tot_orb']=sorted(ele_orb, key=lambda x : orbital_order[str(x[0])+str(x[1])])
            d[label]=row
            data.append(row)
        return d, pd.DataFrame(data)

    def _get_order(self, elements, lm_orbital):
        order = 0
        if isinstance(lm_orbital, dict):
            lm_keys = lm_orbital.keys()
            if any(elem in lm_keys for elem in elements):
                if len(lm_keys) == 2:
                    order = 3
                elif elements[0] in lm_keys:
                    order = 1
                elif elements[1] in lm_keys:
                    order = 2
        return order
    
    def _get_default_orbital_list(self, elem, label):
        """Get default orbital list of given element."""
        orb = self.d[label][elem]["tot_orb"][-1][1]
        return self.d[label][elem][f'{orb}_orb']

    def _get_specific_orbital_list(self, elem, lm_list, label):
        """Get specific orbital list according to lm_orbital."""
        e_temp = []
        for lm in lm_list:
            if lm in ['s', 'p', 'd', 'f']:
                e_temp.extend(self.d[label][elem][f'{lm}_orb'])
            else:
                for orbital in self.d[label][elem].keys():
                    if orbital != 'tot_orb':
                        e_temp.extend([(n, orbital) for n, orbital in self.d[label][elem][orbital] if orbital == lm])
        return e_temp
     
    def _get_orbital_list(self, elem, lm_orbital, label):
        """Get orbital list of given element."""
        if isinstance(lm_orbital, dict):
            if elem in lm_orbital.keys():
                lm_list = lm_orbital[elem]
                return self._get_specific_orbital_list(elem, lm_list, label) 
            else:
                return self._get_default_orbital_list(elem, label) 
        else:
            return self._get_default_orbital_list(elem, label) 
          
    def _get_orbital_lists(self, elements, lm_orbital, label):
        """Get orbital lists of given elements."""
        e = []
        for elem in elements:
            e_temp = self._get_orbital_list(elem, lm_orbital, label)
            e.append(e_temp)
        return e

    def _get_icohps_list(self, label_list, orbital_list, summed_spin_channels):
        """Get ICOHPs list for given label list and orbital list."""
        icohps_list = []
        for i in range(len(label_list)):
            pcohp = self.completecohp.get_summed_cohp_by_label_and_orbital_list(label_list[i], orbital_list[i], summed_spin_channels=summed_spin_channels)
            icohps = pcohp.icohp
            cohpd = pcohp.as_dict()
            e_fermi = cohpd["efermi"]
            energies = cohpd["energies"]
            spins = [Spin.up] if summed_spin_channels else [Spin.up, Spin.down]
            for spin in spins:
                for j in range(len(energies)):
                    if energies[j] == e_fermi:
                        icohp = float(icohps[spin][j])
                        icohps_list.append(icohp)
        return icohps_list
    
    def _get_label_and_orbital_lists(self, elements, e, label, order):
        e1 = e[0]
        e2 = e[1]
        data_label = []
        orbital_list = []
        label_list = []
        porb = []
        idx1=self.d[label]["elem1_idx"]
        idx2=self.d[label]["elem2_idx"]
        for elem in elements:
            orb = self.d[label][elem]["tot_orb"][-1][1]
            porb.append(orb)

        if order == 0: # ex) {'O': 's', 'Mn': 'd'}
            orbital_list = [[[x, y] for x in e1 for y in e2]]
            label_list = [[label for _ in range(len(orbital_list[0]))]]
            data_label = [f'{elements[0]}({porb[0]})-{elements[1]}({porb[1]})']
        elif order == 1: # ex) {'O': ['s', 'p']}
            for x in e1:
                orbital_temp = [[x, y] for y in e2]
                orbital_list.append(orbital_temp)
                label_list.append([label for _ in range(len(orbital_temp))])
                data_label.append(f'{elements[0]}({x[1]})-{elements[1]}({porb[1]})')
            data_label.append(f'{elements[0]}({e1[0][1][0]})-{elements[1]}({porb[1]})')
        elif order == 2: # ex) {'Mn': ['dxy']}
            for y in e2:
                orbital_temp = [[x, y] for x in e1]
                orbital_list.append(orbital_temp)
                label_list.append([label for _ in range(len(orbital_temp))])
                data_label.append(f'{elements[0]}({porb[0]})-{elements[1]}({y[1]})')
            data_label.append(f'{elements[0]}({porb[0]})-{elements[1]}({e2[0][1][0]})')
        elif order == 3: # ex) {'O': ['px','pz'], 'Mn': ['dxy']}
            for x in e1:
                for y in e2:
                    orbital_temp = [[x, y]]
                    orbital_list.append(orbital_temp)
                    label_list.append([label for _ in range(len(orbital_temp))])
                    data_label.append(f'{elements[0]}({x[1]})-{elements[1]}({y[1]})')
            data_label.append(f'{elements[0]}({e1[0][1][0]})-{elements[1]}({e2[0][1][0]})')
        return data_label, orbital_list, label_list
        

    def _get_pcohp(self, label=None, lm_orbital=None, summed_spin_channels=True):
        df = pd.DataFrame()
        elements = self.d[label]["elements"]
        e = self._get_orbital_lists(elements, lm_orbital, label)
        order = self._get_order(elements, lm_orbital)
        data_label, orbital_list, label_list = self._get_label_and_orbital_lists(elements, e, label, order)
        self.d[label]['lm_orbital'] = {}
        icohps_list = []
        key = 1
        data_list= []
        for i in range(len(label_list)):
            pcohp=self.completecohp.get_summed_cohp_by_label_and_orbital_list(label_list[i],orbital_list[i],
                                                                      summed_spin_channels=summed_spin_channels)
            icohps=pcohp.icohp
            cohpd=pcohp.as_dict()
            e_fermi=cohpd["efermi"]
            energies=cohpd["energies"]
            spins=[Spin.up] if summed_spin_channels else [Spin.up,Spin.down]
            icohp_sum_temp=[]

            for spin in spins:
                for j in range(len(energies)):
                    if energies[j]==e_fermi:
                        icohp=float(icohps[spin][j])
                        icohp=round(icohp,5)
                        icohp_sum_temp.append(icohp)
            dat_label=data_label[i]

            icohps_list.extend(icohp_sum_temp)

            dat_label = data_label[i]
            zero = lambda x: abs(x) if abs(x) == 0.000 else x
            row = {
                'label': int(label),
                'ele1': elements[0],
                'idx1': self.d[label]["elem1_idx"],
                'ele2': elements[1],
                'idx2': self.d[label]["elem2_idx"],
                'pair': dat_label,
            }
            if summed_spin_channels:
                row['-ICOHP'] = -icohp_sum_temp[0]
            else:
                row['-ICOHP(up)'] = -icohp_sum_temp[0]
                row['-ICOHP(down)'] = -icohp_sum_temp[1]
            row['distance'] = self.d[label]['length']
            
            self.d[label]['lm_orbital'][key] = {
                "dat_label": dat_label,
                "pcohp": pcohp
            }
            if summed_spin_channels:
                self.d[label]['lm_orbital'][key]["-ICOHP"] = -icohp_sum_temp[0]
            else:
                self.d[label]['lm_orbital'][key]["-ICOHP(up)"] = -icohp_sum_temp[0]
                self.d[label]['lm_orbital'][key]["-ICOHP(down)"] = -icohp_sum_temp[1]
            data_list.append(row)
            key += 1
        icohp_sum=[]
        if summed_spin_channels:
            icohp_sum.append(-sum(icohps_list))
        else:
            icohp_sum.append(-sum(icohps_list[0::2]))
            icohp_sum.append(-sum(icohps_list[1::2]))
        df = pd.DataFrame(data_list)
        # add icohp_sum_all to dataframe
        for i in range(len(label_list)):
            if isinstance(lm_orbital, dict):
                if summed_spin_channels:
                    df.loc[i, '-ICOHP_sum'] = icohp_sum[0]
                else:
                    df.loc[i, '-ICOHP_sum(up)'] = icohp_sum[0]
                    df.loc[i, '-ICOHP_sum(down)'] = icohp_sum[1]
        
        # print(df)
        return df, icohp_sum

    """ print pCOHP """

    def pcohp(self, label=None, lm_orbital=None, summed_spin_channels=True, sort_by=None, index=None):
        df_list = []
        icohp_sum_list=[]
        if isinstance(index,int):
            # find label if index in self.d[label]["elem1_idx"] or self.d[label]["elem2_idx"]
            labels=[]
            for label in self.d.keys():
                if index in [self.d[label]["elem1_idx"], self.d[label]["elem2_idx"]]:
                    labels.append(label)
            label=labels
        if label is None:
            labels = list(self.d.keys())
        elif isinstance(label, list):
            labels = label
        else:
            labels = [label]

        for lb in labels:
            lb = str(lb)
            df, icohp_sum = self._get_pcohp(label=lb, lm_orbital=lm_orbital, summed_spin_channels=summed_spin_channels)
            df_list.append(df)
            # if not isinstance(lm_orbital, dict):
            icohp_sum_list.extend(icohp_sum)
            
        df_total = pd.concat(df_list)

        if sort_by is not None:
            df_total = df_total.sort_values(by=[sort_by,'label'])

        print(df_total.to_string(index=False))
        # if not isinstance(lm_orbital, dict):
        print(f"\t\t\t   -ICOHP sum: {sum(icohp_sum_list):5f}\n")
        return df_total
    
    def print_all(self, sort_by=None):
        """Print all labels and ICOHP values for all orbitals."""
        # change column name of self.dataframe
        self.dataframe = self.dataframe.rename(columns={'icohp': '-ICOHP', 'length': 'distance', 
                                                        'elem1': 'ele1', 'elem2': 'ele2',
                                                        'elem1_idx': 'idx1', 'elem2_idx': 'idx2'})
        new_columns = ['label', 'ele1', 'idx1', 'ele2', 'idx2', 'pair', '-ICOHP', 'distance']
        df = self.dataframe.reindex(columns=new_columns)
        if sort_by is not None:
            df = df.sort_values(by=[sort_by,'label'])
        print(df.to_string(index=False))
        print(f"ICOHP sum: {self.dataframe['-ICOHP'].sum():.5f}")


class CohpPlot:
    def __init__(self, d, sum_orbital=False):
        self.d = d
        plotter = CohpPlotter()
        self.energies_list = []
        self.cohp_list = []
        self.icohp_list = []
        self.dat_label_list = []
        for key in self.d.keys():
            if 'lm_orbital' in self.d[key]:
                for i in self.d[key]['lm_orbital'].keys():
                    self._cohp = self.d[key]['lm_orbital'][i]['pcohp']
                    self.dat_label = self.d[key]['lm_orbital'][i]['dat_label']
                    self.icohp = self.d[key]['lm_orbital'][i]['-ICOHP']
                    plotter.add_cohp(self.dat_label, self._cohp)
                    self.cohp = plotter
                    self.energies_list.append(self.cohp._cohps[self.dat_label]["energies"])
                    self.cohp_list.append(self.cohp._cohps[self.dat_label]["COHP"])
                    self.icohp_list.append(self.cohp._cohps[self.dat_label]["ICOHP"])
                    self.dat_label_list.append(self.dat_label)
        if sum_orbital:
            self._sum_orbital()

    def _sum_orbital(self):
        self.energies_list = self.energies_list[:1]
        for spin in [Spin.up, Spin.down]:
            for i in range(len(self.cohp_list)):
                if spin in self.cohp_list[i]:
                    self.cohp_list[0][spin] += self.cohp_list[i][spin]
            for i in range(len(self.icohp_list)):
                if spin in self.icohp_list[i]:
                    self.icohp_list[0][spin] += self.icohp_list[i][spin]
        self.dat_label_list = self.dat_label_list[:1]
        
    def plot(self,
                  subplot=True,
                  width=None,
                  height=None,
                  dpi=300,
                  xlim=10,
                  ylim=None,
                  fill=True,
                  spins = [Spin.up, Spin.down],
                  plot_negative=True,
                  colors=None,
                  legend_on=True,
                  plt=None,
                  linewidth=1.0,
                  legend_frame_on=False,
                  img_format="png",
                  save_files=True,
                  filename="cohp_plot",
                    ):
        
        if subplot:
            plt = plot_setting.pretty_subplot(1,2,width=width,height=height,dpi=dpi,gridspec_kw=({'wspace': 0.0}),
                                              sharey=True, plt=plt)
        else:
            plt = plot_setting.pretty_plot(width=width,height=height,dpi=dpi, plt=plt)

        if colors is None:
            cycler_obj = plt.rcParams['axes.prop_cycle']
            colors = [item['color'] for item in cycler_obj]

        fig = plt.gcf()

        for i, dat_label in enumerate(self.dat_label_list):
            print(dat_label)
            energies = self.energies_list[i]
            cohp = self.cohp_list[i]
            icohp = self.icohp_list[i]

            if subplot:
                ax1 = fig.axes[0]
                ax2 = fig.axes[1]
                cohps=[cohp, icohp]
                axs=[ax1,ax2]
            else:
                ax1 = fig.axes[0]
                cohps=[cohp]
                axs=[ax1]
            
            for spin in [Spin.up, Spin.down]:
                y = energies
                # print('spin:',spin)
                if spin in cohp:
                    x1 = -cohp[spin] if plot_negative else cohp[spin]
                    ax1.plot(x1, y, label=str(dat_label),color=colors[i],linewidth=linewidth)

                if len(cohps) == 2:
                    # print('icohp')
                    # print(spin)
                    if spin in icohp:
                        x2 = -icohp[spin] if plot_negative else icohp[spin]
                        ax2.plot(x2, y, label=str(dat_label), color=colors[i], linewidth=linewidth)
                        plot_setting.draw_themed_line(0, ax2, orientation="horizontal")

                        for i in range(len(energies)):
                            if energies[i] == 0.0:
                                icohpf=icohp[spin][i]
                        ax2.plot(-icohpf, 0, marker='o',fillstyle='full', markerfacecolor='white', markersize=4, color = 'k')
            
        loc = "upper right" if subplot else "upper right"
        ncol = 1 if subplot else num_columns    


        if legend_on:
            if subplot:
                ax2.legend(loc=loc, frameon=legend_frame_on, ncol=ncol)
            else:
                ax1.legend(loc=loc, frameon=legend_frame_on, ncol=ncol)


        # plot(ax1, xmin=-xmax1,xmax=xmax1)
        ax1.set_ylabel('${E - E}$$_\mathrm{F}$ / eV')
        ax1.set_xlabel('$-\mathrm{pCOHP}$')
        ax1.set_ylim(-xlim,xlim)
        
        fig.text(
            0.03,
            0.03,
            'antibonding',
            ha="left",
            color=rcParams["grid.color"],
            va="center",
            #rotation="vertical",
            transform=ax1.transAxes,
        )
        fig.text(
            0.97,
            0.03, #0.05
            'bonding',
            ha="right",
            color=rcParams["grid.color"],
            va="center",
            #rotation="vertical",
            transform=ax1.transAxes,
        )
        if subplot:
            #plot(ax2, xmin=0,xmax=xmax2+0.5)
            ax1.axvline(0, color='k',linewidth=rcParams["ytick.major.width"])
            ax2.set_xlabel('$-\mathrm{ICOHP}$ / eV')
            fig.subplots_adjust(hspace=0)
            leg=ax2.legend()
            for line in leg.get_lines():
                line.set_linewidth(2.0)
        plt.show()

        if save_files:
            fig.savefig(filename, dpi=dpi, bbox_inches="tight")
            print(f"Figure saved as {filename}")

            # data save as filename.tsv
            filename = filename.split(".")[0]
            dat_filename = filename + ".tsv"
            with open(dat_filename, "w") as f:
                f.write("# E-Ef / eV\tCOHP / eV\tICOHP / eV\n")
                for i, dat_label in enumerate(self.dat_label_list):
                    energies = self.energies_list[i]
                    cohp = self.cohp_list[i]
                    icohp = self.icohp_list[i]
                    for j in range(len(energies)):
                        if spin in cohp:
                            x1 = -cohp[spin] if plot_negative else cohp[spin]
                        if spin in icohp:
                            x2 = -icohp[spin] if plot_negative else icohp[spin]
                        f.write(f"{energies[j]:.5f}\t{x1[j]:.5f}\t{x2[j]:.5f}\n")
            print(f"Data saved as {dat_filename}")

if __name__ == "__main__":
    pass

# Example
# cohpout=Cohpout()
# metal_indices=[i for i in range(len(cohpout.atoms)) if cohpout.atoms[i].symbol in metal_elements]
# for i, metal_index in enumerate(metal_indices):
#     metal=cohpout.atoms[metal_index].symbol
#     print(f"{i}: {metal}({metal_index})")
#     cohpout.pcohp(lm_orbital={metal:'d'}, index=metal_index)


# Plot Example
# cohpout=Cohpout()
# CohpPlot(cohpout.pcohp(index=22),sum_orbital=True).plot(filename='cohp.png')