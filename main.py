import math
from multiprocessing import Pool
import pandas as pd
import matplotlib.pyplot as plt
import sys
from collections import Counter

def read_arguments() -> dict:
    parameters = {}
    input_parameters_command_line = sys.argv
    i = 0 
    for flag in input_parameters_command_line:
        if "-" in flag:
            parameters[flag] = input_parameters_command_line[i+1]
        
        i += 1 

    return parameters 

def read_data():
    data_com = pd.read_csv('/home/anchieta/3TB/Anchieta/Projetos/CG-HIV-Prot-Energy/Dados/DadosClara/com_lig/100_micro/rg_cg_and_rmsd_com_lig.dat', sep='\t', skipinitialspace=True)
    data_sem = pd.read_csv('/home/anchieta/3TB/Anchieta/Projetos/CG-HIV-Prot-Energy/Dados/DadosClara/sem_lig/100_micro/rg_cg_and_rmsd_sem_lig_100.dat', sep='\t', skipinitialspace=True)

    data_com = data_com.assign(model='holo')
    data_sem = data_sem.assign(model='apo')

    data = pd.merge(data_sem, data_com, how='outer')

    return data


def calcule_Delta_G(probability: list, temp: float, model: list) -> list:
    temp = 310
    KB = 3.2976268E-24
    AN = 6.02214179E23
    T = float(temp) #temperature
    Pmax = max(probability)
    LnPmax = math.log(Pmax) 
    d_g = []

    for b in probability:
       e = -0.001*AN*KB*T*(math.log(b)-LnPmax)
       d_g.append(e)


    i = 0
    max_e = max(d_g)
    min_e_holo = max_e
    min_e_apo = max_e
    print(min_e_holo)
    for e in d_g:
        if model[i] == 'holo' and e < min_e_holo:
            min_e_holo = e 
        
        if model[i] == 'apo' and e < min_e_apo:
            min_e_apo = e 
        i +=1

    print(f"Min apo: {min_e_apo}")

    print(f"Min Holo: {min_e_holo}")
    print(f"deltaG: {min_e_holo-(min_e_apo)}")

    return d_g 



def calcule_probability(data: pd) -> list:
    x = data['RMSD'].tolist()
    y = data['RG'].tolist()
    z = []
    x_y = []

    i = 0
    for o in x: 
        x_y.append((o, y[i]))
        i += 1

    counts = Counter(x_y)

    def prob(cont: tuple):
        print(cont)
        return counts[f"{cont}"]

    z = []

    for o in x_y:
        z.append(counts[o])

    p = []
    n = len(z)
    for o in z:
        p.append(o/n)

    return p



def save_data(data_rmsd: list, data_rg:list, data_dg:list, path: str):
    df = pd.DataFrame({'RMSD': data_rmsd, 'RG': data_rg, 'DG': data_dg})
    df.to_csv(f"{path}/teste_data.csv", index=False)



def make_grafic_3D(path: str, data_rg: list, data_rmsd: list, data_dg: list, save_fig: bool, show_grafic: bool):

    fig = plt.figure(figsize=plt.figaspect(0.5))

    ax0 = fig.add_subplot(2, 2, 1, projection='3d')
    hb = ax0.plot_trisurf(data_rmsd[:size_data], data_rg[:size_data], data_dg[:size_data], cmap='OrRd', edgecolor='none',   linewidth=0.5, antialiased=True)
    ax0.set(xlim=xlim, ylim=ylim)
    cb = fig.colorbar(hb, ax=ax0, label='∆G(kcal/mol)')

    ax1 = fig.add_subplot(2, 2, 2, projection='3d')
    hb = ax1.plot_trisurf(data_rmsd[size_data:(2*size_data)], data_rg[size_data:(2*size_data)], data_dg[size_data:(2*size_data)], cmap='Greens', edgecolor='none',   linewidth=0.5, antialiased=True)
    ax1.set(xlim=xlim, ylim=ylim)
    cb = fig.colorbar(hb, ax=ax1, label='∆G(kcal/mol)')

    ax2 = fig.add_subplot(2, 2, 3, projection='3d')
    hb = ax2.plot_trisurf(data_rmsd[:size_data], data_rg[:size_data], data_dg[:size_data], cmap='OrRd', edgecolor='none',   linewidth=0.5, antialiased=True)
    hb = ax2.plot_trisurf(data_rmsd[size_data:(2*size_data)], data_rg[size_data:(2*size_data)], data_dg[size_data:(2*size_data)], cmap='Greens', edgecolor='none',   linewidth=0.5, antialiased=True)
    ax2.set(xlim=xlim, ylim=ylim)
    #cb = fig.colorbar(hb, ax=d_g, label='∆G(kcal/mol)')

    ax3 = fig.add_subplot(2, 2, 4, projection='3d')
    hb = ax3.plot_trisurf(data_rmsd, data_rg, data_dg, cmap='viridis', edgecolor='none',   linewidth=0.5, antialiased=True)
    ax3.set(xlim=xlim, ylim=ylim)
    cb = fig.colorbar(hb, ax=ax3, label='∆G(kcal/mol)')
    plt.tight_layout()

    if save_fig:
        fig.savefig(f"{path}")

    if show_grafic:
        plt.show()


#data = read_data().round(1)
#xlim = data['RMSD'].min(), data['RMSD'].max()
#ylim = data['RG'].min(), data['RG'].max() 
#size_data = int(len(data.index)/2)
#dg = calcule_Delta_G(calcule_probability(data, 1), 310, data['model'].tolist())
#make_grafic_3D("./", data['RG'].tolist(), data['RMSD'].tolist(), dg, False, True)
print(read_arguments())