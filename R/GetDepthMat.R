### PURPOSE: 
### Output a square matrix of size = length(gene.list) 
### Where, the i,j element of matrix is the maximim depth of 
### all terms shared in common by genes i and j 

### NOTE: The following libraries need to be loaded to run this: 
### Biobase, GO, YEAST

GetDepthMat <- function(genes.list, gocategory)
{

	######################################################
	### LOADING LIBRARIES
	######################################################

	library(Biobase)
	library(GO)
	library(YEAST)


	######################################################
	### LOCAL FUNCTIONS
	######################################################

	## Function that takes a GO probe ID and returns a vector of all 
	## BP terms annotated to it.

	getterms<-function(probe, gocat)
	{
        	a<-get(probe,YEASTGO)

		if(length(a) <= 1)
		{
			if (is.na(a))
			{
				return(c("0"))
			}
		}

        	### Keep the GO annotations from the specific Ontology
		
		a1 <- unlist(a)
		a1 <- as.vector(a1)
		

		a2 <- a1[seq(1,length(a1),by=3)]
		

		a2.cat <- a1[seq(3,length(a1),by=3)]
		

		b <- a2[a2.cat == gocat]
		

		if (gocat == "BP")
		{
        		b <- b[b != "GO:0003673"]
			b <- b[b != "GO:0000004"]
		}

		else if (gocat == "MF") 
		{
			 b <- b[b != "GO:0003673"]
        		 b <- b[b != "GO:0005554"]

		}

		else if (gocat == "CC")
		{
			b <- b[b != "GO:0003673"]
        		b <- b[b != "GO:0008372"]

		}

		if (length(b) == 0){
			return(c("0"))
		}

		myvec <- c()
		for (i in 1: length(b))
		{
			myvec <- c(myvec, b[i])
			return.vec <- recurse.terms(b[i], myvec, gocat)
			myvec <- return.vec
		}	

		return(unique(myvec))
	}


	recurse.terms <- function(Node, Container, gocat)
	{
		if (gocat == "BP")
		{
			parents <- get(Node, GOBPPARENTS)
		}

		else if (gocat == "MF") 
		{
			parents <- get(Node, GOMFPARENTS)
		}

		else if (gocat == "CC") 
		{
			parents <- get(Node, GOCCPARENTS)
		} 	


		if (parents[1] == "GO:0003673")
		{
			Container <- Container
		}

		else 
		{
			num.parents <- length(parents)

			for (j in 1:num.parents)
			{
				Container <- c(Container, parents[j])
				Container <- recurse.terms(parents[j], Container, gocat)
			}		
		}
		return(Container)
	}	

	### Function that will return the *MINIMUM* distance to the root for each term

	getdepth <- function(Node, gocat)
	{
		Count <- 0
		Count.List <- c()

		if (gocat == "BP") 
		{
			parents <- get(Node, GOBPPARENTS) 
		}

		else if (gocat == "MF") 
		{
			parents <- get(Node, GOMFPARENTS)
		}

		else if (gocat == "CC") 
		{
			parents <- get(Node, GOCCPARENTS)  
		}


		parents <- parents[parents != "GO:0003673"]
	
		if (length(parents) == 0)
		{
			Count <- Count + 1 
			Count.List <- c(Count.List,Count)
		}	
		else 
		{
			num.parents <- length(parents)
			Count <- Count + 1
	
			for (i in 1: num.parents)
			{
				Count.List  <- c(Count.List, Count)
				Return.List <-  Recurse.Depth(parents[i], Count.List, gocat)
				Count.List  <-  Return.List
			}	
		}	
		#return(Count.List)
		return(min(Count.List))
	}


	Recurse.Depth <- function(Node, Container, gocat)
	{
		Last <- length(Container)
		Current.Count <- Container[Last]	
		Container <- Container[1:Last-1]

		if (gocat == "BP")
		{
			parents <- get(Node, GOBPPARENTS)
		}
	
		else if (gocat == "MF") 
		{
			parents <- get(Node, GOMFPARENTS) 
		}

		else if (gocat == "CC") 
		{
			parents <- get(Node, GOCCPARENTS) 
		}


		if (length(parents) == 1 & parents[1] == "GO:0003673")
		{
			Current.Count <- Current.Count + 1
			Container <- c(Container, Current.Count)		
		}

		else 
		{		
			if(sum(parents == "GO:0003673") >0 )
			{
				Container <- c(Container, rep(Current.Count+1, sum(parents == "GO:0003673")))
				parents  <- parents[!(parents == "GO:0003673")]		
			}

			num.parents <- length(parents)

			if (num.parents > 0)
			{
				Current.Count <- Current.Count + 1	

				for (j in 1:num.parents)
				{
					Container <- c(Container, Current.Count)
					Return.List <- Recurse.Depth(parents[j], Container, gocat)
					Container <- Return.List
				}
			}			
		}
		return(Container)
	}	

	###	FUNCTION THAT WILL RETURN THE DEPTH OF A NODE USING DEPTH-MATRIX

	find.depth <- function(Node, depth.mat) 
	{
		index <- names(depth.mat) == Node
		depth <- depth.mat[index]
		return(depth)
	}

	######################################################
	###### PROCESSING THE INPUT DATA
	######################################################

	terms1<-sapply(genes.list, getterms, gocat = gocategory)
	
	### Remove all terms and corresponding genes with "0" annotations

	lose <- terms1 == "0"
	terms1 <- terms1[!lose]
	genes.list1<-genes.list[!lose]

	terms1.unique <- c()
	for (i  in 1:length(terms1))
	{
		terms1.unique <- c(terms1.unique, terms1[[i]])
	}
	
	terms1.unique <- unique(terms1.unique)

	### Get the depth of each term 

	terms1.depth <- sapply(terms1.unique, getdepth, gocat = gocategory)
	names(terms1.depth) <- terms1.unique

	### OBTAIN MATRIX OF DEPTH OF INTERSECTION OF PAIR OF GENES
	
	depthmat <- matrix(0, length(genes.list1), length(genes.list1))
	row.names(depthmat) <- genes.list1 


	for (i in 2:length(genes.list1))
	{	
		
		for (j in 1:(i-1))
		{
			set1 <- terms1[[i]]
			set2 <- terms1[[j]] 
			keep <- which(set1 %in% set2)
			if (sum(keep) == 0)
			{
				depthmat[i,j] <- 0
			}
			else
			{
				int.set <- set1[keep]
				int.set.depth <- sapply(int.set, find.depth, depth.mat = terms1.depth)
				depthmat[i,j] <- max(int.set.depth)
			}
		}
	}

	#### RETURN MATRIX
	depthmat <- depthmat + t(depthmat)
	return(depthmat)
}

