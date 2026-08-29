import argparse
import math
import os
import subprocess
import sys

import matplotlib.pyplot as plt
import pandas as pd

from app.rmsd_rg_vmd import script_rg_rmsd


def read_arguments() -> dict:
    parser = argparse.ArgumentParser(
        description='Compare RMSD/RG data and estimate relative free energies.'
    )
    parser.add_argument('-data_a')
    parser.add_argument('-data_b')
    parser.add_argument('-cor_a')
    parser.add_argument('-cor_b')
    parser.add_argument('-top_a')
    parser.add_argument('-top_b')
    parser.add_argument('-traj_a')
    parser.add_argument('-traj_b')
    parser.add_argument('-path', default='./')
    parser.add_argument('-prefix_out', default='output')
    parser.add_argument('-sel', default='all')
    parser.add_argument('-bin', type=int, default=1)
    parser.add_argument('-temp', type=float, default=300.0)
    parser.add_argument('-show_grafic', action='store_true')
    parser.add_argument('-save_fig', action='store_true')
    parser.add_argument('-save_data', action='store_true', default=True)
    parser.add_argument('-no_save_data', action='store_false', dest='save_data')

    args = vars(parser.parse_args(sys.argv[1:]))
    if args['bin'] < 0:
        parser.error('-bin must be zero or greater')
    if args['temp'] <= 0:
        parser.error('-temp must be greater than zero')
    return args


def read_data(path_a: str, path_b: str) -> pd.DataFrame:
    frames = []
    for path, model in ((path_a, 'a'), (path_b, 'b')):
        if not os.path.isfile(path):
            raise FileNotFoundError(f'Input file not found: {path}')
        data = pd.read_csv(path, sep='\t', skipinitialspace=True)
        required_columns = {'RG', 'RMSD'}
        missing = required_columns - set(data.columns)
        if missing:
            raise ValueError(
                f"{path} is missing required columns: {', '.join(sorted(missing))}"
            )
        if data.empty:
            raise ValueError(f'Input file is empty: {path}')
        for column in ('RG', 'RMSD'):
            data[column] = pd.to_numeric(data[column], errors='raise')
        if not data[['RG', 'RMSD']].map(lambda value: math.isfinite(value)).all().all():
            raise ValueError(f'RG and RMSD must contain only finite values: {path}')
        frames.append(data[['RG', 'RMSD']].assign(model=model))

    return pd.concat(frames, ignore_index=True)


def calcule_rg_rmsd(
    cord: str,
    top: str,
    traj: str,
    prefix_out: str,
    path: str,
    selection: str,
    timeout: int = 3600,
):
    os.makedirs(path, exist_ok=True)
    script_path = os.path.join(path, f'{prefix_out}_rg_rmsd.tmp')
    output_path = os.path.join(path, f'{prefix_out}.dat')
    log_path = os.path.join(path, f'{prefix_out}_rg_rmsd_vmd.log')

    with open(script_path, 'w') as file_script:
        file_script.write(
            script_rg_rmsd(
                cord=cord,
                top=top,
                traj=traj,
                prefix_out=prefix_out,
                path=path,
                selection=selection,
            )
        )

    try:
        with open(log_path, 'w') as log_file:
            subprocess.run(
                ['vmd', '-e', script_path, '-dispdev', 'text'],
                stdout=log_file,
                stderr=subprocess.STDOUT,
                cwd=path,
                check=True,
                timeout=timeout,
            )
    except FileNotFoundError as error:
        raise RuntimeError('VMD was not found on PATH') from error
    except subprocess.TimeoutExpired as error:
        raise RuntimeError(f'VMD timed out after {timeout} seconds') from error
    except subprocess.CalledProcessError as error:
        raise RuntimeError(f'VMD failed; see log: {log_path}') from error
    finally:
        try:
            os.remove(script_path)
        except FileNotFoundError:
            pass

    if not os.path.isfile(output_path) or os.path.getsize(output_path) == 0:
        raise RuntimeError(f'VMD produced no data: {output_path}')


def calcule_Delta_G(probability: list, temp: float, model: list) -> list:
    if len(probability) == 0 or len(probability) != len(model):
        raise ValueError('Probability and model data must be non-empty and have equal length')
    temperature = float(temp)
    if temperature <= 0 or not math.isfinite(temperature):
        raise ValueError('Temperature must be a finite value greater than zero')
    if any(prob <= 0 or not math.isfinite(prob) for prob in probability):
        raise ValueError('Probabilities must be finite values greater than zero')

    kb = 3.2976268e-24
    avogadro = 6.02214179e23
    max_log_probability = math.log(max(probability))
    return [
        -0.001 * avogadro * kb * temperature * (math.log(prob) - max_log_probability)
        for prob in probability
    ]


def calcule_probability(data: pd.DataFrame) -> list:
    if data.empty:
        raise ValueError('Cannot calculate probabilities from empty data')
    if not {'RMSD', 'RG'}.issubset(data.columns):
        raise ValueError('Data must contain RMSD and RG columns')
    counts = data.groupby(['RMSD', 'RG'], sort=False, dropna=False)['RMSD'].transform('size')
    return (counts / len(data)).tolist()


def save_data(
    data_rmsd: list,
    data_rg: list,
    data_dg: list,
    path: str,
    model: list,
    prefix_out: str,
):
    output_path = os.path.join(path, f'{prefix_out}_out.dat')
    df = pd.DataFrame({'RMSD': data_rmsd, 'RG': data_rg, 'DG': data_dg, 'model': model})
    df.to_csv(output_path, sep='\t', index=False)


def _plot_surface(axis, data, data_dg, cmap):
    points = pd.DataFrame({'RMSD': data['RMSD'], 'RG': data['RG'], 'DG': data_dg}).drop_duplicates(
        subset=['RMSD', 'RG']
    )
    if len(points) < 3:
        axis.scatter(points['RMSD'], points['RG'], points['DG'], cmap=cmap)
        return None
    try:
        return axis.plot_trisurf(
            points['RMSD'],
            points['RG'],
            points['DG'],
            cmap=cmap,
            edgecolor='none',
            linewidth=0.5,
            antialiased=True,
        )
    except ValueError:
        axis.scatter(points['RMSD'], points['RG'], points['DG'], cmap=cmap)
        return None


def make_grafic_3D(data: pd.DataFrame, data_dg: list, save_fig: bool, show_grafic: bool, path: str = './'):
    systems = [data[data['model'] == model].reset_index(drop=True) for model in ('a', 'b')]
    dg_by_model = [
        [data_dg[index] for index, value in enumerate(data['model']) if value == model]
        for model in ('a', 'b')
    ]
    fig = plt.figure(figsize=plt.figaspect(0.5))
    axes = [
        fig.add_subplot(2, 2, 1, projection='3d'),
        fig.add_subplot(2, 2, 2, projection='3d'),
        fig.add_subplot(2, 2, 3, projection='3d'),
        fig.add_subplot(2, 2, 4, projection='3d'),
    ]
    colors = ('OrRd', 'Greens')

    for axis, system, values, cmap in zip(axes[:2], systems, dg_by_model, colors):
        if not system.empty:
            surface = _plot_surface(axis, system, values, cmap)
            if surface is not None:
                fig.colorbar(surface, ax=axis, label='ΔG (kcal/mol)')

    if all(not system.empty for system in systems):
        _plot_surface(axes[2], data, data_dg, 'viridis')
    if not data.empty:
        _plot_surface(axes[3], data, data_dg, 'viridis')

    for axis in axes:
        axis.set(xlabel='RMSD', ylabel='RG')
    plt.tight_layout()
    if save_fig:
        fig.savefig(os.path.join(path, 'teste_3d.png'))
    if show_grafic:
        plt.show()
    plt.close(fig)


def make_grafic_2D(data: pd.DataFrame, save_fig: bool, show_grafic: bool, path: str = './'):
    systems = [data[data['model'] == model] for model in ('a', 'b')]
    fig_2d, axes = plt.subplots(ncols=4, sharey=True, figsize=(20, 4))
    colors = ('Greens', 'OrRd')
    for axis, system, cmap in zip(axes[:2], systems, colors):
        if not system.empty:
            axis.hexbin(system['RMSD'], system['RG'], gridsize=100, bins='log', cmap=cmap)
    if not data.empty:
        axes[2].hexbin(data['RMSD'], data['RG'], gridsize=100, bins='log', cmap='viridis')
        axes[3].scatter(data['RMSD'], data['RG'], s=4)
    for axis in axes:
        axis.set(xlabel='RMSD', ylabel='RG')
    plt.tight_layout()
    if save_fig:
        fig_2d.savefig(os.path.join(path, 'teste_2d.png'))
    if show_grafic:
        plt.show()
    plt.close(fig_2d)
