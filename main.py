from app.tools import read_arguments, read_data, calcule_Delta_G, calcule_probability, make_grafic_3D



print(read_arguments())


data = read_data()

p = calcule_probability(data=data.round(1))

dg = calcule_Delta_G(probability=p, temp=310, model= data['model'].to_list())

make_grafic_3D(data=data, path='./', data_dg=dg, save_fig=False, show_grafic=True)