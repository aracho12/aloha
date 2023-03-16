from mendeleev import element
import sys
import os
from pymatgen.electronic_structure.cohp import CompleteCohp, get_integrated_cohp_in_energy_range
from pymatgen.electronic_structure.core import Spin
import itertools
from ase.io import read

orbital_order={'1s':0,'2s':1,'2p':2,'3s':3,'3p':4,'4s':5,'3d':6,'4p':7,'5s':8,'4d':9,'5p':10,'6s':11,'4f':12,
               '5d':13,'6p':14,'7s':15,'5f':16,'6d':17,'7p':18}

def rm_num(string):
    new=''.join([i for i in string if not i.isdigit()])
    return new

class Cohpout:
    def __init__(self,
                 filepath='.',
                 d={},
                 element_index=0,
                ):
        COHPCAR_path = filepath+"/COHPCAR.lobster"
        POSCAR_path = filepath+"/POSCAR"
        
        self.cohpcar=COHPCAR_path
        self.poscar=POSCAR_path
        self.completecohp=CompleteCohp.from_file(fmt="LOBSTER", filename=COHPCAR_path, structure_file=POSCAR_path)
        self.atoms=read(POSCAR_path)
        
        info={}
        completecohp=self.completecohp
        
        for i in range(len(self.completecohp.bonds)):
            label=str(i+1)
            info[label]={}

            """ finding element index based on POSCAR (therefore, started from 0) """
            elem1=str(self.completecohp.bonds[label]['sites'][0]._species)[:-1]
            elem2=str(self.completecohp.bonds[label]['sites'][1]._species)[:-1]
            x1,y1,z1=self.completecohp.bonds[label]['sites'][0].coords
            elem1_idx=[atom.index for atom in self.atoms if x1-0.01<atom.x<x1+0.01 and y1-0.01<atom.y<y1+0.01 and z1-0.01<atom.z<z1+0.01]
            x2,y2,z2=self.completecohp.bonds[label]['sites'][1].coords
            elem2_idx=[atom.index for atom in self.atoms if x2-0.01<atom.x<x2+0.01 and y2-0.01<atom.y<y2+0.01 and z2-0.01<atom.z<z2+0.01]
            info[label][elem1+str(elem1_idx[0])]={}
            info[label][elem2+str(elem2_idx[0])]={}
            info[label]['pair']="{}({})-{}({})".format(elem1,elem1_idx[0],elem2,elem2_idx[0])
            
            """ bond length """
            length=str(self.completecohp.bonds[label]['length'])
            info[label]['length']=float(length)

            """ ICOHP value """
            icohp=get_integrated_cohp_in_energy_range(self.completecohp,label) #relative_E_fermi=True, summed_spin_channels=True
            info[label]['icohp']=float(icohp)

            """ print whole label and ICOHP values for all orbital """ 
        db=[]
        print("label \t pair        \t -ICOHP  \t distance ")
        icohp_sum=[]
        for label in info.keys():
            print(f"{label} \t {info[label]['pair']} \t {float(-info[label]['icohp']):5f} \t {float(info[label]['length']):5f} \t")
            icohp_sum.append(round(float(-info[label]['icohp']),5))
        print(f"\t    -ICOHP sum:\t {sum(icohp_sum):5f} ")


        d=info.copy()
        
        """ """
        for label in info.keys():
            elements=[ele for ele in list(info[label].keys()) if ele not in ["icohp","pair","length"]]
            d[label]['elements']=elements
            for elem in elements:
                ele=element(rm_num(elem)).ec.conf
                ele_orb=list(ele.keys())
                d[label][elem]["tot_orb"]=ele_orb
                for i in range(len(ele_orb)):
                    n=int(ele_orb[i][0])
                    n_orb=ele_orb[i][1]
                    #print(n_orb)
                    if n_orb == 's':
                        d[label][elem]['s_orb']=[(n,"s")]
                    if n_orb == 'p':
                        d[label][elem]['p_orb']=[(n,"px"),(n,"py"),(n,"pz")]
                    if n_orb == 'd':
                        d[label][elem]['d_orb']=[(n,"dxy"),(n,"dyz"),(n,"dz2"),(n,"dxz"),(n,"dx2")]
                    if n_orb == 'f':
                        d[label][elem]['f_orb']=[(n,"f_3"),(n,"f_2"),(n,"f_3"),(n,"f0"),(n,"f1"),(n,"f2"),(n,"f3")]
                d[label][elem]["tot_orb"]=sorted(ele_orb, key=lambda x : orbital_order[str(x[0])+str(x[1])])
            self.d=d


    def get_pcohp(self,label=None,lm_orbital=None,summed_spin_channels=True):
        elements=self.d[label]["elements"] # i.g; [O79, Mn78]
        porb=[]
        label=str(label)
        label_list=[]
        orbital_list=[]
        e=[]
        pcohp_d={}
        
        """ 
        
        default pCOHP
        i.e) d-elements - p-elements : d-p orbitals
        
        """
        
        order=0
        for i in range(len(elements)):
            if isinstance(lm_orbital,dict):
                if rm_num(elements[i]) in lm_orbital.keys():
                    if len(lm_orbital.keys())==2:
                        order=3
                    elif i==0:
                        order=1
                    elif i==1:
                        order=2
                    pcohp_d[rm_num(elements[i])]=lm_orbital[rm_num(elements[i])]
                else:
                    pcohp_d[rm_num(elements[i])]="default"
            else:
                pcohp_d[rm_num(elements[i])]="default"
        # print(pcohp_d)

        for elem in elements:
            e_temp=[]
            if pcohp_d[rm_num(elem)]=="default":
                orb=self.d[label][elem]["tot_orb"][-1][1]
                if orb == 's':
                    e_temp.extend(self.d[label][elem]['s_orb'])
                    porb.append('s')
                elif orb == 'p': 
                    e_temp.extend(self.d[label][elem]['p_orb'])
                    porb.append('p')
                elif orb == 'd':
                    e_temp.extend(self.d[label][elem]['d_orb'])
                    porb.append('d')
                elif orb == 'f':
                    e_temp.extend(self.d[label][elem]['f_orb'])
                    porb.append('f')
            else:
                if isinstance(pcohp_d[rm_num(elem)],str):
                    lm_list=[pcohp_d[rm_num(elem)]]
                else:
                    lm_list=pcohp_d[rm_num(elem)]
                for lm in lm_list:
                    if lm in ['s','p','d','f']:
                        num=0
                        if lm == 's':
                            e_temp.extend(self.d[label][elem]['s_orb'])
                        elif lm =='p':
                            e_temp.extend(self.d[label][elem]['p_orb'])
                        elif lm == 'd':
                            e_temp.extend(self.d[label][elem]['d_orb'])
                        elif lm == 'f':
                            e_temp.extend(self.d[label][elem]['f_orb'])
                    else:
                        for orbital in self.d[label][elem].keys():
                            if orbital == 'tot_orb':
                                pass
                            else:
                                e_temp.extend([(n, orbital) for n, orbital in self.d[label][elem][orbital] if orbital == lm])
            e.append(e_temp)
        # print(e)
        e1=e[0]
        e2=e[1]
        
        """ get label list and orbital list for summed cohp """
        data_label=[]
        if order == 0:
            orbital_list=[[[x,y] for x in e1 for y in e2]]
            label_list=[[label for i in range(len(orbital_list[0]))]]
            data_label=[f'{rm_num(elements[0])}({porb[0]})-{rm_num(elements[1])}({porb[1]})']
        elif order == 1:
            for x in e1:
                orbital_temp=[[x,y] for y in e2]
                orbital_list.append(orbital_temp)
                label_list.append([label for i in range(len(orbital_temp))])
                data_label.append(f'{rm_num(elements[0])}({x[1]})-{rm_num(elements[1])}({porb[0]})')
            data_label.append(f'{rm_num(elements[0])}({e1[0][1][0]})-{rm_num(elements[1])}({porb[0]})')
        elif order == 2:
            for y in e2:
                orbital_temp=[[x,y] for x in e1]
                orbital_list.append(orbital_temp)
                label_list.append([label for i in range(len(orbital_temp))])
                data_label.append(f'{rm_num(elements[0])}({porb[0]})-{rm_num(elements[1])}({y[1]})')
            data_label.append(f'{rm_num(elements[0])}({porb[0]})-{rm_num(elements[1])}({e2[0][1][0]})')
        elif order ==3:
            for x in e1:
                for y in e2:
                    orbital_temp=[[x,y]]
                    orbital_temp2=[x,y]
                    orbital_list.append(orbital_temp)
                    label_list.append([label for i in range(len(orbital_temp))])
                    data_label.append(f'{rm_num(elements[0])}({x[1]})-{rm_num(elements[1])}({y[1]})')
            data_label.append(f'{rm_num(elements[0])}({e1[0][1][0]})-{rm_num(elements[1])}({e2[0][1][0]})')
        
        # print(orbital_list)
        # print(label_list)
        # print(data_label)
        
        self.d[label]['lm_orbital']={}
        
        
        """ label for plot """
        
        icohp_sum=[]
        key=1
        for i in range(len(label_list)):
            pcohp=self.completecohp.get_summed_cohp_by_label_and_orbital_list(label_list[i],orbital_list[i],
                                                                      summed_spin_channels=summed_spin_channels)        
            icohps=pcohp.icohp
            cohpd=pcohp.as_dict()
            e_fermi=cohpd["efermi"]
            energies=cohpd["energies"]
            
            if summed_spin_channels:
                spins=[Spin.up]
            else:
                spins=[Spin.up, Spin.down]
            icohps_list=[]
            
            for spin in spins:
                for j in range(len(energies)):
                    if energies[j] == e_fermi:
                        icohp=float(icohps[spin][j])
                        icohp_sum.append(icohp)
                        icohps_list.append(icohp)
            
            dat_label=data_label[i]
            
            zero = lambda x : abs(x) if abs(x)==0.000 else x
        
            if isinstance(lm_orbital,dict):
                if summed_spin_channels:
                    print(f"{label} \t {dat_label} \t {zero(float(-icohps_list[0])):5f} ")
                else:
                    print(f"{label} \t {dat_label} \t {zero(float(-icohps_list[0])):5f}     \t {zero(float(-icohps_list[1])):5f}")
            else:
                if summed_spin_channels:
                    print(f"{label} \t {dat_label} \t {zero(float(-icohps_list[0])):5f}    \t {self.d[label]['length']:5f}")
                else:
                    print(f"{label} \t {dat_label} \t {zero(float(-icohps_list[0])):5f}    \t {zero(float(-icohps_list[1])):5f}     \t {round(self.d[label]['length'],5)}")
            if summed_spin_channels:
                self.d[label]['lm_orbital'][key]={
                    "dat_label":dat_label,
                    "icohp":[float(icohps_list[0])],
                    "pcohp":pcohp}
            else:
                self.d[label]['lm_orbital'][key]={
                    "dat_label":dat_label,
                    "icohp":[float(icohps_list[0]),float(icohps_list[1])],
                    "pcohp":pcohp}
            key+=1
        
        if isinstance(lm_orbital,dict):
            if summed_spin_channels:
                print(f"\t    -ICOHP sum:\t {-sum(icohp_sum):5f}\t distance:  {self.d[label]['length']:5f}\n")
            else:
                print(f"\t    -ICOHP sum:\t {-sum(icohp_sum[0::2]):5f}\t {-sum(icohp_sum[1::2]):5f}\t distance:  {self.d[label]['length']:5f}\n")
                
        return self.d[label]['lm_orbital']

    def pcohp(self, label=None, lm_orbital=None, summed_spin_channels=True):
        # Header 
        if isinstance(lm_orbital,dict):
            if summed_spin_channels:
                print("label \t pair        \t -ICOHP     \t")
            else:
                print("label \t pair        \t -ICOHP(up) \t -ICOHP(down) \t")
        else:
            if summed_spin_channels:
                print("label \t pair        \t -ICOHP     \t distance \t")
            else:
                print("label \t pair        \t -ICOHP(up) \t -ICOHP(down) \t distance \t")
        
        if label==None:
            icohp_sum=[]
            for label in self.d.keys():
                label=str(label)
                pcohp=self.get_pcohp(label=label,lm_orbital=lm_orbital,summed_spin_channels=summed_spin_channels)
                if not isinstance(lm_orbital, dict):
                    icohp_sum.extend(pcohp[1]["icohp"])
            if not isinstance(lm_orbital,dict):        
                if summed_spin_channels:
                    print(f"\t    -ICOHP sum:\t {-sum(icohp_sum):5f}\n")
                else:
                    print(f"\t    -ICOHP sum:\t {-sum(icohp_sum[0::2]):5f} \t {-sum(icohp_sum[1::2]):5f}\n")
        else:
            if isinstance(label,list):
                icohp_sum=[]
                for lb in label:
                    lb=str(lb)
                    pcohp=self.get_pcohp(label=lb,lm_orbital=lm_orbital,summed_spin_channels=summed_spin_channels)
                    if not isinstance(lm_orbital, dict):
                        icohp_sum.extend(pcohp[1]["icohp"])
                if not isinstance(lm_orbital,dict):
                    if summed_spin_channels:
                        print(f"\t    -ICOHP sum:\t {-sum(icohp_sum):5f}\n")
                    else:
                        print(f"\t    -ICOHP sum:\t {-sum(icohp_sum[0::2]):5f} \t {-sum(icohp_sum[1::2]):5f}\n")
            else:
                label=str(label)
                pcohp=self.get_pcohp(label=label,lm_orbital=lm_orbital,summed_spin_channels=summed_spin_channels)

if __name__ == "__main__":
    pass
