install.packages("bio3d")

# Can you improve this analysis code?
library(bio3d)
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

# Change it to read easier with one vector "x"
x <- read.pdb("x")
x.chainA <- trim.pdb(x, chain = "A", elety = "CA")
x.b <- x.chainA$atom$b
plotb3(x.b, sse = x.chainA, typ = "l", ylab = "Bfactor")



  ## Q1. What type of object is returned from the read.pdb() function? 
    ## read.pdb() returns a protein data bank coordinate file. 

  ## Q2. What does the trim.pdb() function do?
    ## trim.pdb() produces a smaller PDB object containing a subset of atoms from a larger pdb object.

  ## Q3. What input parameter would turn off the marginal black and grey rectangles in the
  ## plots and what do they represent in this case?
    ## The sse argument of the plotb3() function shows the secondary structure object in the margin of the 
    ## plot of the residue versus the bfactor.

  ## Q4. What would be a better plot to compare across the different proteins? 
    ## A plot that has all three proteins on the same plot would be better to compare the proteins. 
    ## rbind() can be used to combine all three proteins into one dataset. dist() can be used to 
    ## determine the distances between the rows of the data. And hclust() will analyse the data and produce
    ## hierarchical clustering.
      hc <- hclust(dist(rbind(s1.b, s2.b, s3.b)))
      plot(hc)
    
  ## Q5. Which proteins are more similar to each other in their B-factor trends. How could you quantify this?
      ## According to the dendrogram, "1AKE" and "1E4Y" are more alike than the other protein when their
      ## B-factor is taken into account. 
      
      dist(rbind(s1.b, s2.b, s3.b))
      ## The distance between the b-factors of "1AKE" and "1E4Y" are 255.36, much lower than 340.49 and 
      ## 458.71 (the distances between "4AKE" and "1AKE" and "1E4Y", respectively.)
      
  ## Q6. How would you ogeneralize the original code above to work with any set of input protein structures?
      ## On a new R-script