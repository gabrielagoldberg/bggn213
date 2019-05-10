Class 12
================

## Structural Bioinformatics (Part 2)

## Section 1: In silico docking of drugs to HIV-1 protease

# 1.1 Obtaining and inspecting our input structure

``` r
library(bio3d)
```

    ## Warning: package 'bio3d' was built under R version 3.4.4

``` r
file.name <- get.pdb("1hsg")
```

    ## Warning in get.pdb("1hsg"): ./1hsg.pdb exists. Skipping download

Read the PDB file

``` r
hiv <- read.pdb(file.name)
```

Get a quick summary of the pdb structure

``` r
hiv
```

    ## 
    ##  Call:  read.pdb(file = file.name)
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

Q1: What is the name of the two non protein resid values in this
structure? What does resid correspond to and how would you get a listing
of all reside values in this structure? \> HOH and MK1 are non-protein
resid values in this structure. \> resid corresponds to the residues.

# 1.2 Prepare initial protein and ligand input files

Make protein-only and ligand-only objects called prot and lig that you
can then write out to a new PDB formal files.

``` r
prot <- trim.pdb(hiv, "protein")
lig <- trim.pdb(hiv, "ligand")

write.pdb(prot, file = "1hsg_protein.pdb")
write.pdb(lig, file = "1hsg_ligand.pdb")
```

# 1.3 Using AutoDockTools to setup protein docking input

In ADT (AutoDocTools), load the protein using File \> Read Molecule.
Select 1hsg\_protein.pdb. Click Open

We opened out protein only PDB file in AutoDocTools and added hydrogens
and atom-trpes needed for docking calculations.

We will use AutoDoc Vina here at the UNIX command line\!

# 1.4 Prepare the ligand

In ADT (AutoDocTools), load the ligand using File \> Read Molecule.
Select 1hsg\_ligand.pdb. Click Open

We opened out ligand only PDB file in AutoDocTools and added hydrogens
and atom-trpes needed for docking calculations.

# 1.5 Prepare a docking configuration file

We need to create an input file that defines the protein, ligand and the
search parameters. We will create the input file in a text editor.

In the new text file, we will put the grid dimensions found in ADT.

## Section 2: Docking ligands into HIV-1 protease

We used a program called Autodock Vina ran it from the command line.

‘“Files (x86)Scripps Research Institute.exe” –config config –log
log.txt’

# 2.3 Inspecting your docking results

Process the all.pdbqt to a PDB formal file that can be loaded into VMD.

``` r
res <- read.pdb("all.pdbqt", multi = TRUE)
write.pdb(res, file = "results.pdb")
```

Now both the original ‘1hsg.pdb’ and the new ‘results.pdb’ file can be
loaded into VMD.

To assess the results quantitatively we will calculate the RMSD between
the docking results and the knock crystal structure.

``` r
ori <- read.pdb("1hsg_ligand.pdbqt")
rmsd(ori, res)
```

    ##  [1]  0.716 11.087 10.521 10.857 10.611 10.958  4.058  4.679 10.313  3.703
    ## [11] 10.975  2.941 10.415  5.780  6.207  2.496 11.790 11.499 11.330 10.549

## Section 3: Exploring the conformational dynamics of proteins

# 3.1 Normal Mode Analysis (NMA)

``` r
pdb <- read.pdb("1HEL")
```

    ##   Note: Accessing on-line PDB file

``` r
modes <- nma(pdb)
```

    ##  Building Hessian...     Done in 0.09 seconds.
    ##  Diagonalizing Hessian...    Done in 0.13 seconds.

``` r
m7 <- mktrj(modes, mode = 7, file = "nma_7.pdb")
plot(modes, sse = pdb)
```

![](class12_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->
