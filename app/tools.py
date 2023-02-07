import math
import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import sys
from collections import Counter
from app.rmsd_rg_vmd import script_rg_rmsd_big_traj, script_rg_rmsd, script_bigdcd

def read_arguments() -> dict:
    parameters = {}
    input_parameters_command_line = sys.argv
    i = 0 
    for flag in input_parameters_command_line:
        if "-" in flag:
            try:
                if "\"" in input_parameters_command_line[i+1]:
                    text = ""
                    o = 1
                    while "\"" != input_parameters_command_line[i+1]:
                        text += input_parameters_command_line[i+1+o]
                        o += 1
                else:
                    parameters[flag] = input_parameters_command_line[i+1]

            except:
                parameters[flag] = True
        
        i += 1 
    return parameters 

def read_data(path_a: str, path_b:str):
    #data_com = pd.read_csv('/home/anchieta/3TB/Anchieta/Projetos/CG-HIV-Prot-Energy/Dados/DadosClara/com_lig/100_micro/rg_cg_and_rmsd_com_lig.dat', sep='\t', skipinitialspace=True)
    #data_sem = pd.read_csv('/home/anchieta/3TB/Anchieta/Projetos/CG-HIV-Prot-Energy/Dados/DadosClara/sem_lig/100_micro/rg_cg_and_rmsd_sem_lig_100.dat', sep='\t', skipinitialspace=True)

    data_a = pd.read_csv(path_a, sep='\t', skipinitialspace=True)
    data_b = pd.read_csv(path_b, sep='\t', skipinitialspace=True)

    data_a = data_a.assign(model='a')
    data_b = data_b.assign(model='b')

    data = pd.merge(data_a, data_b, how='outer')

    return data


def calcule_rg_rmsd(cord: str, top: str, traj: str, prefix_out:str, path: str, selection:str, big_traj: str):
    if big_traj == "on": 
        with open(f'{path}/{prefix_out}_rg_rmsd.tmp', 'w') as file_script:
            file_script.writelines(script_rg_rmsd_big_traj(cord=cord, top=top, traj=traj, prefix_out=prefix_out, path=path, selection=selection))
        with open(f'{path}/bigdcd.tcl', 'w') as file_bigdcd:
            file_bigdcd.writelines(script_bigdcd())
    else:
        with open(f'{path}/{prefix_out}_rg_rmsd.tmp', 'w') as file_script:
            file_script.writelines(script_rg_rmsd(cord=cord, top=top, traj=traj, prefix_out=prefix_out, path=path, selection=selection))

    proc = subprocess.Popen(args=f"vmd -e {path}/{prefix_out}_rg_rmsd.tmp -dispdev text > {path}/{prefix_out}_rg_rmsd_vmd.log", shell=True)
    proc.wait()
    try:
        os.remove(f'{path}/{prefix_out}_rg_rmsd.tmp')
        os.remove(f'{path}/bigdcd.tcl')
    except:
        pass

def calcule_Delta_G(probability: list, temp: float, model: list) -> list:
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
        if model[i] == 'a' and e < min_e_holo:
            min_e_holo = e 
        
        if model[i] == 'b' and e < min_e_apo:
            min_e_apo = e 
        i +=1

    print(f"Min A: {min_e_apo}")

    print(f"Min B: {min_e_holo}")
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



def make_grafic_3D(data: pd,  data_dg: list, save_fig: bool, show_grafic: bool, path: str = "./"):
    data_rmsd = data['RMSD'].tolist()
    data_rg = data['RG'].tolist()
    xlim = data['RMSD'].min(), data['RMSD'].max()
    ylim = data['RG'].min(), data['RG'].max() 
    size_data = int(len(data.index)/2)

    fig = plt.figure(figsize=plt.figaspect(0.5))

    ax0 = fig.add_subplot(2, 2, 1, projection='3d')
    hb = ax0.plot_trisurf(data_rmsd[:size_data], data_rg[:size_data], data_dg[:size_data], cmap='OrRd', edgecolor='none',   linewidth=0.5, antialiased=True)
    ax0.set(xlim=xlim, ylim=ylim, xlabel="RMSD", ylabel="RG")
    cb = fig.colorbar(hb, ax=ax0, label='∆G(kcal/mol)')

    ax1 = fig.add_subplot(2, 2, 2, projection='3d')
    hb = ax1.plot_trisurf(data_rmsd[size_data:(2*size_data)], data_rg[size_data:(2*size_data)], data_dg[size_data:(2*size_data)], cmap='Greens', edgecolor='none',   linewidth=0.5, antialiased=True)
    ax1.set(xlim=xlim, ylim=ylim, xlabel="RMSD", ylabel="RG")
    cb = fig.colorbar(hb, ax=ax1, label='∆G(kcal/mol)')

    ax2 = fig.add_subplot(2, 2, 3, projection='3d')
    hb = ax2.plot_trisurf(data_rmsd[:size_data], data_rg[:size_data], data_dg[:size_data], cmap='OrRd', edgecolor='none',   linewidth=0.5, antialiased=True)
    hb = ax2.plot_trisurf(data_rmsd[size_data:(2*size_data)], data_rg[size_data:(2*size_data)], data_dg[size_data:(2*size_data)], cmap='Greens', edgecolor='none',   linewidth=0.5, antialiased=True)
    ax2.set(xlim=xlim, ylim=ylim, xlabel="RMSD", ylabel="RG")
    #cb = fig.colorbar(hb, ax=d_g, label='∆G(kcal/mol)')

    ax3 = fig.add_subplot(2, 2, 4, projection='3d')
    hb = ax3.plot_trisurf(data_rmsd, data_rg, data_dg, cmap='viridis', edgecolor='none',   linewidth=0.5, antialiased=True)
    ax3.set(xlim=xlim, ylim=ylim, xlabel="RMSD", ylabel="RG")
    cb = fig.colorbar(hb, ax=ax3, label='∆G(kcal/mol)')
    plt.tight_layout()

    if save_fig:
        fig.savefig(f"{path}/teste.png")

    if show_grafic:
        plt.show()

def make_grafic_2D(data: pd,  data_dg: list, save_fig: bool, show_grafic: bool, path: str = "./"):
    data_rmsd = data['RMSD'].tolist()
    data_rg = data['RG'].tolist()
    xlim = data['RMSD'].min(), data['RMSD'].max()
    ylim = data['RG'].min(), data['RG'].max() 
    size_data = int(len(data.index)/2)

    fig_2D, (ax0, ax1, ax2, ax3) = plt.subplots(ncols=4, sharey=True, figsize=(20, 4))
    # plot 1 
    hb = ax0.hexbin(data_rmsd[:size_data], data_rg[:size_data], gridsize=100, bins='log', cmap='Greens')
    ax0.set(xlim=xlim, ylim=ylim, xlabel="RMSD", ylabel="RG")

    #cb = fig.colorbar(hb, ax=ax0, label='log10(N)')
    # Plot 2 
    hb = ax1.hexbin(data_rmsd[size_data:(2*size_data)], data_rg[size_data:(2*size_data)], gridsize=100, bins='log', cmap='OrRd')
    ax1.set(xlim=xlim, ylim=ylim, xlabel="RMSD", ylabel="RG")

    #cb = fig.colorbar(hb, ax=ax1, label='log10(N)')
    # plot 3 
    hb = ax2.hexbin(data_rmsd[:size_data], data_rg[:size_data], gridsize=100, bins='log', cmap='Greens')
    hb = ax2.hexbin(data_rmsd[size_data:(2*size_data)], data_rg[size_data:(2*size_data)], gridsize=100, bins='log', cmap='OrRd')
    ax2.set(xlim=xlim, ylim=ylim, xlabel="RMSD", ylabel="RG")

    #cb = fig.colorbar(hb, ax=ax2, label='log10(N)')
    # Plot 4
    hb = ax3.hexbin(data_rmsd, data_rg, gridsize=100, bins='log', cmap='viridis')
    ax3.set(xlim=xlim, ylim=ylim, xlabel="RMSD", ylabel="RG")
    plt.tight_layout()

    #cb = fig.colorbar(hb, ax=ax3, label='log10(N)')
    if save_fig:
        fig_2D.savefig(f"{path}/teste.png")

    if show_grafic:
        plt.show()
