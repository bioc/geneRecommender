\name{gr.normalize}
\alias{gr.normalize}
\title{A function for normalization}
\description{
  \code{gr.normalize} normalizes a matrix of gene expression data
  as part of the implementation of the Gene Recommender algorithm
  described in Owen et al (2003).  \code{gr.normalize} must be applied
  to the data before running \code{gr.main}.  
}
\usage{
gr.normalize(unnormalized.dataset)
}
\arguments{
  \item{unnormalized.dataset}{
    A matix or exprSet containing the normalized gene expression data.  
    The rows correspond to genes, the columns correspond to experiments, and the
    entries correspond to the gene expression levels.  The rows must be labeled.  
  } 
}
\details{
  \code{gr.normalize} normalizes the data so that for each gene, the
  gene expression measurements are distributed uniformly between 
  -1 and 1.  
}
\value{
  The normalized gene expression data, in the same format as the input.  
}
\references{
  Art B. Owen, Josh Stuart, Kathy Mach, Anne M. Villeneuve, and Stuart Kim.  
  A Gene Recommender Algorithm to Identify Coexpressed Genes in C. elegans.  
  Genome Research 13:1828-1837, 2003.  
}
\author{
  Gregory J. Hather \email{ghather@stat.berkeley.edu} \cr
  with contributions from from Art B. Owen \email{art@stat.stanford.edu}\cr
  and Terence P. Speed \email{terry@stat.berkeley.edu}.  
}
\seealso{
  \code{gr.main, gr.cv}
}
\examples{
#This example uses the geneData dataset from the Biobase package
data(geneData)
my.query <- c("31730_at", "31331_at", "31712_at", "31441_at")
normalized.data <- gr.normalize(geneData)
}
\keyword{manip}
