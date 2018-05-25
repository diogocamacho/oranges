\name{oranges}
\alias{oranges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
ORANGES: Over Representation ANalyse for GEne Signatures
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This package will perform an over-representation analysis for a query gene signature based on Fisher's exact test (or, the hypergeometric distribution test).  The pathway sets were aggregated from the Consensus Pathway DataBase (CPDB) and from the Molecular Signatures DataBase (MSigDB).  Currently this function is specific for human data.
}
\usage{
oranges(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{query_entrez}{A gene set of interest in EntrezIDs.}
  \item{universe_entrez}{The complete set of EntrezIDs where the signature originated from.}
  \item{universe_names}{Gene symbols for the universe set.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
The function returns a tibble data frame, filtered to pathways with no genes present in query set.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Diogo M. Camacho
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
