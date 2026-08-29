import os
from multiprocessing import Process

from app.tools import (
    calcule_Delta_G,
    calcule_probability,
    calcule_rg_rmsd,
    make_grafic_2D,
    make_grafic_3D,
    read_arguments,
    read_data,
    save_data,
)


def _validate_trajectory_flags(args):
    for suffix in ('a', 'b'):
        flags = {f'cor_{suffix}', f'top_{suffix}', f'traj_{suffix}'}
        provided = {flag for flag in flags if args[flag]}
        if provided and provided != flags:
            missing = ', '.join(sorted(flags - provided))
            raise SystemExit(
                f"Incomplete trajectory input for system {suffix}; missing: {missing}"
            )


def _validate_data_arguments(args):
    missing = [f'data_{suffix}' for suffix in ('a', 'b') if not args[f'data_{suffix}']]
    if missing:
        raise SystemExit(
            'Two input datasets are required; missing: ' + ', '.join(missing)
        )


def _run_trajectory_jobs(args):
    processes = []
    for suffix in ('a', 'b'):
        if all(args[f'{name}_{suffix}'] for name in ('traj', 'cor', 'top')):
            prefix = f"{args['prefix_out']}_{suffix}"
            process = Process(
                target=calcule_rg_rmsd,
                args=(
                    args[f'cor_{suffix}'],
                    args[f'top_{suffix}'],
                    args[f'traj_{suffix}'],
                    prefix,
                    args['path'],
                    args['sel'],
                ),
            )
            process.start()
            processes.append((suffix, process, os.path.join(args['path'], f'{prefix}.dat')))
            args[f'data_{suffix}'] = os.path.join(args['path'], f'{prefix}.dat')

    for suffix, process, output_path in processes:
        process.join()
        if process.exitcode != 0:
            raise RuntimeError(f"VMD processing failed for system {suffix}")
        if not os.path.isfile(output_path) or os.path.getsize(output_path) == 0:
            raise RuntimeError(f"VMD produced no data for system {suffix}: {output_path}")


def main():
    args = read_arguments()
    _validate_trajectory_flags(args)
    os.makedirs(args['path'], exist_ok=True)
    _run_trajectory_jobs(args)
    _validate_data_arguments(args)

    data = read_data(path_a=args['data_a'], path_b=args['data_b'])
    rounded_data = data.copy()
    rounded_data[['RMSD', 'RG']] = rounded_data[['RMSD', 'RG']].round(args['bin'])
    probability = calcule_probability(data=rounded_data)
    delta_g = calcule_Delta_G(
        probability=probability,
        temp=args['temp'],
        model=data['model'].tolist(),
    )

    if args['show_grafic']:
        make_grafic_3D(data=rounded_data, data_dg=delta_g, save_fig=False, show_grafic=True)
        make_grafic_2D(data=data, save_fig=False, show_grafic=True)

    if args['save_fig']:
        make_grafic_3D(
            data=rounded_data,
            path=args['path'],
            data_dg=delta_g,
            save_fig=True,
            show_grafic=False,
        )
        make_grafic_2D(data=data, path=args['path'], save_fig=True, show_grafic=False)

    if args['save_data']:
        save_data(
            data_rmsd=data['RMSD'].tolist(),
            data_rg=data['RG'].tolist(),
            data_dg=delta_g,
            model=data['model'].tolist(),
            path=args['path'],
            prefix_out=args['prefix_out'],
        )


if __name__ == '__main__':
    main()
