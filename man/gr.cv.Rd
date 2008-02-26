\name{gr.cv}
\alias{gr.cv}
\title{A function for cross validation}
\description{
  \code{gr.cv} performs leave-one-out cross validation with \code{gr.main} 
  for each element of the query.  
}
\usage{
gr.cv(normalized.dataset, query, fun = median)
}
\arguments{
  \item{normalized.dataset}{
    A matrix or ExpressionSet containing the normalized gene expression data.  
    The rows correspond to genes, the columns correspond to experiments, and the
    entries correspond to the gene expression levels.  The rows must be labeled.  
    The values contained in \code{normalized.dataset} must be either finite or NA.  
  }
  \item{query}{A vector containing the query set of genes.  These should correspond to the row names of \code{normalized.dataset}.  The query must contain at least 2 elements}
  \item{fun}{A function used in choosing the number of experiments to include in the calculation.  See the help file for \code{gr.main} for details.  }
}
\details{
  In addition to measuring performance, the results of the cross validation
  can be used to determine if some element(s) in the query might not belong.  
  If one of the elements in the output vector was very large,
  one would suspect that the associated gene was regulated differently
  than the other genes in the query.  
}
\value{
  A vector containing the rank 
  of each element in the query
  produced by applying \code{gr.main} to the query with
  that element removed.  
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
  \code{gr.main, gr.normalize}
}
\examples{
#This example uses the geneData dataset from the Biobase package
data(geneData)
my.query <- c("31730_at", "31331_at", "31712_at", "31441_at")
normalized.data <- gr.normalize(geneData)
gr.cv(normalized.data, my.query)
}
\keyword{manip}
