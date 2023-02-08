from multiprocessing import Process
from app.tools import read_arguments, read_data, calcule_Delta_G, calcule_probability, make_grafic_3D, save_data, calcule_rg_rmsd, make_grafic_2D


args = read_arguments()


if '-traj_a' in args and '-cor_a' in args and '-top_a' in args:
    # Criar data RG-RMSD
    proc_a = Process(target=calcule_rg_rmsd, args=(args['-cor_a'], args['-top_a'], args['-traj_a'], args['-prefix_out'] + "_a", args['-path'], args['-sel'], args['-big_traj']))
    proc_a.start()
    args['-data_a']=args['-path'] + "/" + args['-prefix_out'] + "_a.dat"


if '-traj_b' in args and '-cor_b' in args and '-top_b' in args:
    # Criar data RG-RMSD
    proc_b = Process(target=calcule_rg_rmsd, args=(args['-cor_b'], args['-top_b'], args['-traj_b'], args['-prefix_out']+ "_b", args['-path'], args['-sel'], args['-big_traj']))
    proc_b.start()
    args['-data_b']=args['-path'] + "/" + args['-prefix_out'] + "_b.dat"


if '-data_DG' in args or '-show_grafic' in args or '-save_fig' in args or '-save_data' in args:
    if "-traj_a" in args: 
        proc_a.join()    
    if "-traj_b" in args:
        proc_b.join()
    
    data = read_data(path_a=args['-data_a'], path_b=args['-data_b'])
    probability = calcule_probability(data=data.round(int(args['-bin'])))
    delta_g = calcule_Delta_G(probability=probability, temp=args['-temp'], model = data['model'].tolist())


if '-show_grafic' in args or '-save_fig' in args:
    if '-show_grafic' in args:
        make_grafic_3D(data=data.round(int(args['-bin'])), data_dg=delta_g, save_fig=False, show_grafic=args['-show_grafic'])
        make_grafic_2D(data=data, data_dg=delta_g, save_fig=False, show_grafic=args['-show_grafic'])

    if '-save_fig' in args:
        make_grafic_3D(data=data.round(int(args['-bin'])), path=args['-path'], data_dg=delta_g, save_fig=args['-save_fig'], show_grafic=False)
        make_grafic_2D(data=data, path=args['-path'], data_dg=delta_g, save_fig=args['-save_fig'], show_grafic=False)


if '-save_data' in args:
    # Salva valores de RG-RMSD-DG em arquivo
    save_data(data_rmsd=data['RMSD'].tolist(), data_rg=data['RG'].tolist(), data_dg=delta_g, model=data['model'].tolist(), path=args['-path'], prefix_out=args['-prefix_out'])

