\name{gr.main}
\alias{gr.main}
\title{A gene recommender algorithm to identify genes coexpressed with a query set of genes}
\description{
  The function \code{gr.main} implements the Gene Recommender algorithm
  described in Owen et al (2003).  Note that in order for \code{gr.main} 
  to provide meaningful results, 
  the normalization procedure \code{gr.normalize}
  must first be applied to the gene expression matrix.  
}
\usage{
gr.main(normalized.dataset, query, fun = median, ngenes = NULL, extra = FALSE)
}
\arguments{
  \item{normalized.dataset}{
    A matrix or ExpressionSet containing the normalized gene expression data.  
    The rows correspond to genes, the columns correspond to experiments, and the
    entries correspond to the gene expression levels.  The rows must be labeled.  
    The values contained in \code{normalized.dataset} must be either finite or NA.  
  }
  \item{query}{A vector containing the query set of genes.  These should correspond to the row names of \code{normalized.dataset}.  The query must contain at least 2 elements}
  \item{fun}{A function used in choosing the number of experiments to include in the calculation.  See below for details.  }
  \item{ngenes}{The number of genes to return in the result.  It's
    default value is the number of genes found at 50 percent recall.  }
  \item{extra}{
    A logical value.  When false, the output list will contain only one item, 
    \code{result}.  When true, several other quantities (listed below) will be calculated 
    and added to the output list.  
  }
}
\details{
  Given data from a large number of microarray experiments and a query set of genes,
  which genes have expression profiles that are similar to the query?  The Gene Recommender
  algorithm (Owen et al, 2003) answers this question by first identifying the set of
  experiments over which the query genes behave similarly.  Next, the algorithm ranks
  all the genes based on the strength of the correlation with the query across the
  chosen set of experiments.  

  The algorithm must choose how generous to be in including experiments.  How many
  experiments should be included?  The algorithm tries every number of experiments
  and chooses the number which minimizes a score.  In the paper, the score was
  defined as the median of the ranks of the query genes.  In \code{gr.main},
  the score can be computed with the user-defined function, \code{fun}.  
}
\value{
  A list containing entries:
  \item{result}{
    An array of dimensions (\code{ngenes}, 2).  Column 1 contains the resulting genes, 
    with the highest scoring genes listed first.  Column 2 contains character strings, 
    indicating whether the corresponding gene is from the query list or not.  
  }
  \item{fifty.percent.recall}{Number of genes found at 50 percent recall.  }
  \item{experiments.included}{Experiments included in the analysis.  }
  \item{experiments.excluded}{Experiments excluded from the analysis.  }
  \item{s.g.i}{An array used as a measure of biological significance for each gene.  The output is ranked by this quantity.  }
  \item{z.g.i}{An array used as a measure of statistical significance for each gene.  }
  \item{contribution}{
    An array indicating the contribution of each experiment to each gene result. 
    For a given gene and a given experiment, the \code{contribution} indicates how
    strongly the experiment suggests that the gene should be high ranking.  
    Using notation from the article, \code{contribution} is defined as 
    \eqn{\bar{Y}_{Q,j} \times Y_{ij}}.  
  }
}
\note{
  The results of \code{gr.main} will differ from the results generated 
  from the C code released by Owen et al (2003). 
  This is due to differences in the implementation.  
  See the vignette for details.  
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
  \code{gr.normalize, gr.cv}
}
\examples{
#This example uses the geneData dataset from the Biobase package
data(geneData)
my.query <- c("31730_at", "31331_at", "31712_at", "31441_at")
normalized.data <- gr.normalize(geneData)
gr.main(normalized.data, my.query, ngenes = 10)
}
\keyword{manip}
