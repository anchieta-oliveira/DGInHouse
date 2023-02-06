from app.tools import read_arguments, read_data, calcule_Delta_G, calcule_probability, make_grafic_3D, save_data

args = read_arguments()


if '-traj' in args and '-cor' in args and '-top' in args:
    # Criar data RG-RMSD
    pass


if '-data_a' in args and '-data_b' in args:
    # Ler dados de RG-RMSD de arquivo
    data = read_data(path_a=args['-data_a'], path_b=args['-data_b']).round(1)


if '-data_DG' in args:
    #Plota gráfico de dados salvos anteriormente
    probability = calcule_probability(data=data)
    delta_g = calcule_Delta_G(probability=probability, temp=args['-temp'], model = data['model'].tolist())


if '-show_grafic' in args or '-save_fig' in args:
    probability = calcule_probability(data=data)
    delta_g = calcule_Delta_G(probability=probability, temp=args['-temp'], model = data['model'].tolist())
    if '-show_grafic' in args:
        make_grafic_3D(data=data, data_dg=delta_g, save_fig=False, show_grafic=args['-show_grafic'])

    if '-save_fig' in args:
        make_grafic_3D(data=data, path=args['-path'], data_dg=delta_g, save_fig=args['-save_fig'], show_grafic=False)


if '-save_data' in args:
    # Salva valores de RG-RMSD-DG em arquivo
    save_data(data_rmsd=data['RMSD'].tolist(), data_rg=data['RG'].tolist(), data_dg=delta_g, path=args['-path'])

