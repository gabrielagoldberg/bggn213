library(bio3d)
# Q6. How would you generalize the original code above to work with any set of input protein structures?

## Original Code:
s1 <- read.pdb("4AKE") # kinase with drug
s2 <- read.pdb("1AKE") # kinase no drug
s3 <- read.pdb("1E4Y") # kinase with drug
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
s3.chainA <- trim.pdb(s3, chain="A", elety="CA")
s1.b <- s1.chainA$atom$b
s2.b <- s2.chainA$atom$b
s3.b <- s3.chainA$atom$b
plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")
plotb3(s2.b, sse=s2.chainA, typ="l", ylab="Bfactor")
plotb3(s3.b, sse=s3.chainA, typ="l", ylab="Bfactor")

## Change it to read easier with one vector "x"
    # x <- read.pdb("x")
    # x.chainA <- trim.pdb(x, chain = "A", elety = "CA")
    # x.b <- x.chainA$atom$b
    # plotb3(x.b, sse = x.chainA, typ = "l", ylab = "Bfactor")

## Write the function!


# This function analyzes protein drug interactions by reading in any protien PDB data and outputs 
# a plot for the specified protein showing the B-factor versus the residue of the protein. 
# The function has one argument "pdb.filename" which is the name of the PDB file to be analyzed. 
bfactor.plot <- function(pdb.filename) {
  
  # This statement reads the PDB file and restricts the original PDB file to only include the atoms/chains 
  # we want to analyze. The argument MUST be in quotes.
  # In this case, we are focusing our analysis on alpha carbons and creating an object "x.chainA"
  x.chainA <- trim.pdb(read.pdb(pdb.filename), chain = "A", elety = "CA")

  # This statement further specifies the rows of the trimmed PDB file to analyze. 
  # In this case, we are focusing on the B-factor of the atoms in the trimmed PDB file we created in the 
  # previous line and creating an object "x.b"
  x.b <- x.chainA$atom$b
  
  # This statement produces the output of our function. 
  # Here, we create a line plot of the protein's residue (from x.b) versus that atom's B-factor. 
  # The plot also includes secondary structure analyses in the margins of the plot as black and grey lines.
  plotb3(x.b, sse = x.chainA, typ = "l", ylab = "Bfactor", main = pdb.filename)
}


# Does the function work?? Let's see!
bfactor.plot("4AKE")
bfactor.plot("1AKE")
bfactor.plot("1E4Y")
