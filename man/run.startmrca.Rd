\name{run.startmrca}
\alias{run.startmrca}
\title{
Estimating Time to The Common Ancestor for a Beneficial Allele
}
\description{
This Method uses Markov chain Monte Carlo to estimate time to the common ancestor for a beneficial allele.  
}
\usage{
run.startmrca(vcf.file, rec.file, mut.rate, rec.rate,nsel, nanc, chain.length, proposal.sd, nanc.post, sample.ids, refsample.ids, pos, sel.allele, bed.file, upper.t.limit)
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
Uses a bed file to exclude sites that have not been sequenced in the VCF file. 
}
  \item{upper.t.limit}{
An upper limit from which to draw a uniformly distributed initial value for t before starting the MCMC.  
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

\item{t.chain }{Description of 't.chain'} a four column matrix of the MCMC output for the TMRCA. The number of rows is equal to the number specified by the chain.length argument. The first column is the proposed TMRCA value, the second column is the emission probability of the current TMRCA and ancestral haplotype values. The third column is the Metropolis-Hastings ratio. The fourth column indicates whether the proposed value was accepted or not.  
\item{anchap.chain }{Description of 'anchap.chain'} a sample from the posterior distribution of the ancestral haplotype. The number of rows is specified by the nanc.post argument, and indicates the size of the posterior sample from the end of the chain. 

}
\references{
Joel Smith, Matthew Stephens, Graham Coop and John Novembre (in prep): Estimating Time to the Common Ancestor for a Beneficial Allele.
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

