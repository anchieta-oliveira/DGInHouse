# DGInHouse

DGInHouse compares molecular-dynamics systems using the root-mean-square deviation
(RMSD) and radius of gyration (RG). It estimates relative free-energy values from
the observed RMSD/RG-state probabilities and produces comparative 2D and 3D plots.
Typical applications include protein-protein or protein-ligand interaction studies,
protein folding, and mutation analysis.

This repository is intended as a lightweight research utility. Results should be
interpreted together with the simulation protocol and validated for the scientific
question being studied.

## Requirements

- Python 3
- `pandas` and `matplotlib` (see `requirements.txt`)
- VMD with the `pbctools` package, only when trajectories are processed directly

Install the Python dependencies with:

```console
python3 -m pip install -r requirements.txt
```

## Input-data workflow

The simplest workflow uses two tab-separated files, one for each system. Each file
must contain the columns `RG` and `RMSD`:

```console
python3 main.py -data_a examples/system_a.tsv -data_b examples/system_b.tsv
```

The default output is `<prefix_out>_out.dat` in `-path`. The default prefix is
`output` when trajectory processing is requested.

Alternatively, DGInHouse can generate the input tables with VMD. All three flags
must be supplied for each system:

```console
python3 main.py \
  -cor_a system_a.pdb -top_a system_a.psf -traj_a system_a.dcd \
  -cor_b system_b.pdb -top_b system_b.psf -traj_b system_b.dcd \
  -path results -prefix_out comparison -sel "protein"
```

Use `-big_traj` for the streaming trajectory workflow. Its periodic-boundary
alignment behavior is a known limitation and is not guaranteed to match the normal
trajectory workflow.

## Command-line options

- `-data_a`, `-data_b`: tab-separated RG/RMSD input files.
- `-cor_a`/`-cor_b`: coordinate files such as PDB or GRO.
- `-top_a`/`-top_b`: topology files such as PSF.
- `-traj_a`/`-traj_b`: molecular-dynamics trajectory files.
- `-path`: directory for generated data and figures; default `./`.
- `-prefix_out`: prefix for generated output; default `output` in trajectory mode.
- `-sel`: VMD atom selection; default `all` in trajectory mode.
- `-bin`: decimal precision used to group RMSD/RG states; default `1`.
- `-temp`: temperature used in the free-energy calculation; default `300`.
- `-show_grafic`: display the plots interactively.
- `-save_fig`: save `teste_3d.png` and `teste_2d.png` under `-path`.
- `-save_data`: save RMSD, RG, free-energy, and system labels; enabled by default.

Boolean options are enabled by including the flag. The command-line parser is
intentionally lightweight and does not provide a separate help command.

## Output and figures

The 3D and 2D comparison examples are available in
[`documentation/fig1.png`](documentation/fig1.png) and
[`documentation/fig2.png`](documentation/fig2.png).

## Citation

If this software contributes to published work, cite the repository and the version
used. The machine-readable citation metadata is in [`CITATION.cff`](CITATION.cff).
For a persistent DOI, create a GitHub release and enable the repository in Zenodo.

## License

DGInHouse is distributed under the MIT License. See [`LICENSE`](LICENSE).
