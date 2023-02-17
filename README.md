# DGInHouse

DGInHous visa automatizar as análises de Root-mean-square deviation of atomic positions (RMSD) e Raio de Giro (RG), e realizar cálculos de Energia Livre, representando graficamente a superfície de energia (figura 1 e 2) comparativos entre dois sistemas (ex.: apo e holo). Podendo ser útil em estudos de interação Proteína-Proteína ou Proteína-Ligante, enovelamento proteico, mutações, dentre outros. 

## Utilização 

```console
python3 main.py -data_a caminho_data_a -data_b caminho_data_b
```

## Opções de comandos 
 - `-data_a` ou `data_b`: arquivo contendo valores de RG e RMSD;
 - `-cor_a` ou `cor_b`: arquivo de coordenada da simulação (pdb, .gro, etc);
 - `-top_a` ou `top_b`: arquivo de topologia da simulação (.psf)
 - `-traj_a` ou `traj_b`: trajetória da Dinâmica Molecular;
 - `-path`: local para salvar dados e/ou gráfico;
 - `-show_grafic`: mostrar gráfico na janela;
 - `-save_fig`: salva figura no caminho definido;
 - `-save_data`: salva valores de RMSD, RG e DG;
 
 
Figura 1: Gráfico 3D RMSD, RG e Energia Livre  
![Figura 1](https://github.com/anchieta-oliveira/DGInHouse/blob/dev/documentation/fig1.png)
Figura 1: Gráfico 2D RMSD, RG e Energia Livre  
![Figura 2](https://github.com/anchieta-oliveira/DGInHouse/blob/dev/documentation/fig2.png)
