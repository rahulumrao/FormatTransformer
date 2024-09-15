# FormatTransformer
# Convert XYZ structure to VASP/QE Input files

## How to

### Download :

```bash
git clone
```

### Install :

```bash
cd format_transformer

pip install .
```

### Run :

```bash
❯ xyzconverter -h

usage: xyzconverter [-h] [-i INPUT] [-b BOX_SIZE [BOX_SIZE ...]] -o {VASP,QE} [-pp PSEUDO]

Convert XYZ to VASP POSCAR or Quantum Espresso input files

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input XYZ file containing multiple frames
  -b BOX_SIZE [BOX_SIZE ...], --box_size BOX_SIZE [BOX_SIZE ...]
                        Box size in Å (default: Cubic box of 12.0)
  -o {VASP,QE}, --output_type {VASP,QE}
                        Specify the output format: VASP or QE
  -pp PSEUDO, --pseudo PSEUDO
                        PseudoPotentials directory (required for Espresso)
```

```bash
❯ xyzcenter -h
usage: xyzcenter [-h] [-i INPUT] [-o OUTPUT] [-b BOX_SIZE]

Center molecular structure within PBC

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input XYZ trajectory file
  -o OUTPUT, --output OUTPUT
                        Output XYZ trajectory file
  -b BOX_SIZE, --box_size BOX_SIZE
                        Box size in Ang. (default: 15.0)
```
