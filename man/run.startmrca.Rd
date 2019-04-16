\name{run.startmrca}
\alias{run.startmrca}
\title{
Estimating Time to The Common Ancestor for a Beneficial Allele
}
\description{
This Method uses Markov chain Monte Carlo to estimate time to the common ancestor for a beneficial allele.  
}
\usage{
run.startmrca(vcf.file, rec.file, mut.rate, rec.rate = NULL, 
    nsel = NULL, nanc = NULL, chain.length = 15000, proposal.sd = 20, 
    nanc.post = 100, sample.ids, refsample.ids = NULL, pos, 
    sel.allele = 1, bed.file = NULL, upper.t.limit = 2000, 
    output.file = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vcf.file}{
a vcf file that includes genotypes of individuals at genomic region of interest surrounding the putatively selected. A 1 Mb region is recommended (500kb flanking both sides of the seleced allele). This file should contain individuals with and without the selected allele.  
}
  \item{rec.file}{
a text file containing a recombination map with 3 columns: (1) chromosome, (2) position i, and (3) recombination rate between position i and i+1. This may include the recombination map for the entire genome. if rec.file = NULL, a uniform recombination rate is used from a specified rec.rate parameter value.  
}
  \item{mut.rate}{
the mutation rate per basepair per generation.
}
  \item{rec.rate}{
the recombination rate per basepair per generation.
}
  \item{nsel}{
the sample size for the number of chromosomes with the selected allele to use.
}
  \item{nanc}{
the sample size for the reference panel of chromosomes without the selected allele. 
}
  \item{chain.length}{
the number of iterations for the MCMC to run.
}
  \item{proposal.sd}{
the standard deviation for the proposal distribution of the time to the common ancester (TMRCA). Where t is the current parameter value for the tmrca, proposals for a new value, t', are as follows: t' = t + rnorm(1,0,proposal.sd) 
}
  \item{nanc.post}{
the number of values from the end of the chain to use for the posterior distribution of the ancestral haplotype.
}
  \item{sample.ids}{
a single column text file for individuals in the "selected population." if sample.ids = NULL, all individuals in the vcf are used.
}
  \item{refsample.ids}{
a single column text file for individuals in the "reference population." if refsample.ids = NULL, reference individuals are used from sample.ids.
}
  \item{pos}{
the position of the putatively selected allele.
}
  \item{sel.allele}{
specifies which allele is the selected allele. Takes value 0 or 1.
}
  \item{bed.file}{
uses a bed file to exclude sites that have not been sequenced in the VCF file. 
}
  \item{upper.t.limit}{
an upper limit from which to draw a uniformly distributed initial value for t before starting the MCMC.  
}
  \item{output.file}{
a character vector with a file path. If NULL, a default path will be used (VCF file name with "_mcmc_list.RDATA" extension). If NA, no file will be created.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

The function will return a list of two elements:
\item{t.chain}{a four column matrix of the MCMC output for the TMRCA. The number of rows is equal to the number specified by the chain.length argument. The first column is the proposed TMRCA value, the second column is the emission probability of the current TMRCA and ancestral haplotype values. The third column is the Metropolis-Hastings ratio. The fourth column indicates whether the proposed value was accepted or not.}
\item{anchap.chain}{a sample from the posterior distribution of the ancestral haplotype. The number of rows is specified by the nanc.post argument, and indicates the size of the posterior sample from the end of the chain.}

The output of the function can be assigned to an object independently of being saved in a file.

}
\references{
Joel Smith, Graham Coop, Matthew Stephens and John Novembre (2018): Estimating Time to the Common Ancestor for a Beneficial Allele. Molecular Biology and Evolution. https://doi.org/10.1093/molbev/msy006
}
\author{
Joel Smith
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

\dontrun{run.startmrca(vcf.file      = "examples/FIN_chr2_pos136608646.vcf.gz",
              rec.file      = "examples/decode_recmap_sexaveraged.txt", 
              sample.ids    = "examples/sample_ids.txt", 
              refsample.ids = "examples/sample_ids.txt",
              mut.rate      = 1.6e-8, 
              nsel          = 50,
              nanc          = 20,
              chain.length  = 20,
              nanc.post     = 10,
              pos           = 136608646,
              sel.allele    = 1)

}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

