from aloha.cohp_analysis import *
from aloha import plot_setting
import matplotlib.pyplot as plt
from matplotlib import rcParams
import re

class CohpPlot:
    def __init__(self, pcohp=None, sum_orbital=False, sum_labels=False):
        if pcohp is None:
            self.pcohp=Cohpout().pcohp()
        else:
            self.pcohp = pcohp
        for i, label in enumerate(self.pcohp.keys()):
            self.summed_spin = self.pcohp[label][0].summed_spin
            break
        if sum_orbital:
            self._sum_orbital()

        if sum_labels:
            self.pcohp=self._sum_labels()
        
        self.filename='cohp.png'
        
        dat_num = 0
        for i, label in enumerate(self.pcohp.keys()):
            for j in self.pcohp[label].keys():
                dat_num += 1
        self.dat_num = dat_num
        self.label_num = len(list(self.pcohp.keys()))
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

    def _sum_labels(self):
        self.filename = 'cohp_sum_labels.png'
        new_pcohp = {}
        # sum all labels in the same label ; pcohp[first_label]
        for i, label in enumerate(self.pcohp.keys()):
            for j in self.pcohp[label].keys():
                if i == 0:
                    first_label = label
                    new_pcohp[first_label] = self.pcohp[label]
                    continue
                else:
                    if self.summed_spin:
                        new_pcohp[first_label][j].cohp['Spin.up']+=self.pcohp[label][j].cohp['Spin.up']
                        new_pcohp[first_label][j].icohp['Spin.up']+=self.pcohp[label][j].icohp['Spin.up']
                        new_pcohp[first_label][j].icohp_fermi+=self.pcohp[label][j].icohp_fermi
                    else:
                        new_pcohp[first_label][j].cohp['Spin.up']+=self.pcohp[label][j].cohp['Spin.up']
                        new_pcohp[first_label][j].cohp['Spin.down']+=self.pcohp[label][j].cohp['Spin.down']
                        new_pcohp[first_label][j].icohp['Spin.up']+=self.pcohp[label][j].icohp['Spin.up']
                        new_pcohp[first_label][j].icohp['Spin.down']+=self.pcohp[label][j].icohp['Spin.down']
                        new_pcohp[first_label][j].icohp_fermi['Spin.up']+=self.pcohp[label][j].icohp_fermi['Spin.up']
                        new_pcohp[first_label][j].icohp_fermi['Spin.down']+=self.pcohp[label][j].icohp_fermi['Spin.down']
                first_data_label = self.pcohp[first_label][j].data_label
                new_data_label = re.sub(r'\d', '', first_data_label)
                new_pcohp[first_label][j].data_label = new_data_label

        return new_pcohp

        
    def plot(self,
        prefix=None,
        width=3,
        height=4,
        dpi=300,
        xmin=-13,
        xmax=13,
        ymin=None,
        ymax=None,
        show_icohp=True, # True: show ICOHP in the legend
        plot_negative=True,
        colors=None,
        legend_out=False, # True: outside upper right, default: inside upper right
        legend_on=True,
        legend_frame_on=True,
        plt=None,
        save_files=True,
        filename=None,
        grid_spec=({'wspace': 0.0}),
        linewidth=0.7
        
    ):
        """
        Plot COHP and ICOHP.
        Args:
            prefix: Prefix for the filename.
            width: Width of the figure.
            height: Height of the figure.
            dpi: DPI of the figure.
            xmin: Minimum energy to plot.
            xmax: Maximum energy to plot.
            ymin: Minimum COHP to plot.
            ymax: Maximum COHP to plot.
            show_icohp: Whether to show ICOHP in the legend.
            plot_negative: Whether to plot negative COHP and ICOHP.
            colors: Colors to use for the plot.
            legend_out: Whether to put the legend outside the plot.
            legend_on: Whether to show the legend.
            legend_frame_on: Whether to show the legend frame.
            plt: Matplotlib.pyplot object to use for plotting.
            save_files: Whether to save the plot.
            filename: Filename for the plot.
            grid_spec: Grid specification for the plot.
            linewidth: Line width for the plot.

        """
        if filename is None:
            filename = self.filename
        
        spins = ["Spin.up", "Spin.down"]
        num_columns=len(self.pcohp.keys())
        plt = self._subplot(width=width, height=height, dpi=dpi, plt=plt, gridspec_kw=grid_spec)

        if colors is None:
            cycler_obj = plt.rcParams['axes.prop_cycle']
            colors = [item['color'] for item in cycler_obj]

        fig = plt.gcf()
        # default:
        ax1 = fig.axes[0]
        ax2 = fig.axes[1]
        for i, label in enumerate(self.pcohp.keys()):
            for k, key in enumerate(self.pcohp[label].keys()):
                cp = self.pcohp[label][key]
                cohp = cp.cohp
                icohp = cp.icohp
                data_label = cp.data_label
                icohp_fermi = cp.icohp_fermi
                if not type(icohp_fermi) == dict:
                    icohp_fermi = {'Spin.up': icohp_fermi}
                energies = cp.energies
                efemi = cp.efermi
                if len(cohp.keys()) == 1:
                    marker = ['o']
                else:
                    marker = ['^','v']
                for j, spin in enumerate(spins):
                    if self.label_num == 1:
                        i=k
                    if spin in cohp:
                        x = -cohp[spin] if plot_negative else cohp[spin]
                        if spin == "Spin.up":
                            ax1.plot(x, energies, label=data_label, color=colors[i], linewidth=linewidth)
                        else:
                            ax1.plot(-x, energies, label=data_label, color=colors[i], linestyle='--', linewidth=linewidth)
                        ax1.axvline(0, color='k',linewidth=rcParams["ytick.major.width"])
                        
                    if spin in icohp:
                        x = -icohp[spin] if plot_negative else icohp[spin]
                        if spin == 'Spin.up':
                            ax2.plot(x, energies, label=data_label, color=colors[i], linewidth=linewidth)
                        else:
                            ax2.plot(-x, energies, label=data_label, color=colors[i], linestyle='--', linewidth=linewidth)
                        
                        if show_icohp:
                            if spin == 'Spin.up':
                                ax2.plot(icohp_fermi[spin], 0, marker=marker[j],fillstyle='full', markerfacecolor='white', markersize=3, color = colors[i], label=f'-ICOHP: {icohp_fermi[spin]:.3f}', linewidth=linewidth,markeredgewidth=linewidth)
                            else:
                                ax2.plot(-icohp_fermi[spin], 0, marker=marker[j],fillstyle='full', markerfacecolor='white', markersize=3, color = colors[i], label=f'-ICOHP: {icohp_fermi[spin]:.3f}', linewidth=linewidth,markeredgewidth=linewidth)
                        else:
                            if spin == 'Spin.up':
                                ax2.plot(icohp_fermi[spin], 0, marker=marker[j],fillstyle='full', markerfacecolor='white', markersize=3, color = colors[i],linewidth=linewidth,markeredgewidth=linewidth)
                            else:
                                ax2.plot(-icohp_fermi[spin], 0, marker=marker[j],fillstyle='full', markerfacecolor='white', markersize=3, color = colors[i],linewidth=linewidth,markeredgewidth=linewidth)

        plot_setting.draw_themed_line(0, ax2, orientation="horizontal")
        plot_setting.draw_themed_line(0, ax1, orientation="horizontal")
        ax1.set_ylabel('${E - E}$$_\mathrm{F}$ / eV')
        ax1.set_xlabel('$-\mathrm{pCOHP}$')
        ax1.set_ylim(xmin,xmax)
        xticks = ax1.get_xticks()
        ax1.set_xlim(min(xticks),max(xticks))
        ax1.set_xticks(xticks[xticks < max(xticks)])
        if len(cohp.keys()) == 1:
            ax2.set_xlim(0,)
        else:
            ax2.axvline(0, color='k',linewidth=rcParams["ytick.major.width"])

        if self.dat_num > 2:
            legendfontsize = 6
            legend_out = True
            if self.dat_num > 10:
                legendfontsize = 4
        else:
            legendfontsize = 6
            legend_out = False

        if legend_on:
            ax2.legend(fontsize=legendfontsize)
        if legend_out:
            ax2.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, fontsize=legendfontsize)
        ax2.set_xlabel('$-\mathrm{ICOHP}$ / eV')
        plt.show()

        if save_files:
            fig.savefig(filename, dpi=dpi, bbox_inches="tight")
            print(f"Figure saved as {filename}")

            # data save as filename.tsv
            filename = filename.split(".")[0]
            dat_filename = filename + ".tsv"

            with open(dat_filename, "w") as f:
                for i, label in enumerate(self.pcohp.keys()):
                    for j in self.pcohp[label].keys():
                        cp = self.pcohp[label][j]
                        cohp = cp.cohp
                        icohp = cp.icohp
                        data_label = cp.data_label
                        icohp_fermi = cp.icohp_fermi
                        energies = cp.energies
                        efemi = cp.efermi
                        f.write(f"# {data_label}\n")
                        if len(cohp.keys()) == 1:
                            f.write("# E-Ef / eV\tCOHP / eV\tICOHP / eV\n")
                            for i, e in enumerate(energies):
                                f.write(f"{e}\t{cohp['Spin.up'][i]}\t{icohp['Spin.up'][i]}\n")
                        else:
                            f.write("# E-Ef / eV\tCOHP(up) / eV\tCOHP(down) / eV\tICOHP(up) / eV\tICOHP(down) / eV\n")
                            for i, e in enumerate(energies):
                                f.write(f"{e}\t{cohp['Spin.up'][i]}\t{cohp['Spin.down'][i]}\t{icohp['Spin.up'][i]}\t{icohp['Spin.down'][i]}\n")
            print(f"Data saved as {dat_filename}")

    def _subplot(self, 
            nrows=1,
            ncols=1,
            width=None, 
            height=None, 
            dpi=300, 
            plt=None,
            sharex=False,
            sharey=True,
            gridspec_kw=None,
            ):

        if width is None:
            width = rcParams["figure.figsize"][0]
            
        if height is None:
            height = rcParams["figure.figsize"][1]

        # default:
        nrows, ncols = 1 , 2
        sharex, sharey = False, True
        
        if plt is None:
            import matplotlib.pyplot
            plt = matplotlib.pyplot
            plt.subplots(nrows,
                        ncols,
                        sharex=sharex,
                        sharey=sharey,
                        dpi=dpi,
                        figsize=(width, height),
                        gridspec_kw=gridspec_kw,
            )
        
        return plt


if __name__ == "__main__":
    pass
