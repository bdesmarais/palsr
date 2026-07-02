---
output: pdf_document
fontsize: 12pt
---

\thispagestyle{empty}
\today

Editor   
The R Journal  
\bigskip

Dear Editor,
\bigskip

Please consider our article titled "palsr: Projecting the Locations of Moving Actors for Spatial Models of Dyadic Interaction" for publication as a contributed research article in the R Journal.

The manuscript introduces the package `palsr`, which implements the Projected Actor Location (PALS) method for spatial modeling of interactions between geographically mobile actors. Many actors studied across the social sciences --- armed groups, firms, diplomats, migrating populations --- have no fixed location, and the place where two such actors interact is itself an outcome of interest. Standard spatial tooling in R (for example `sf` and `spatstat`) operates on objects with given coordinates and offers no way to infer a moving actor's effective location. `palsr` fills this gap: it projects where an actor "is" at any time from the spatiotemporal history of its past interactions, tunes those projections to predict future interaction locations, and constructs the dyadic distance covariates used in downstream models, with performance-critical kernels implemented in C++ via `Rcpp`.

We believe readers of the R Journal will find this article useful because it makes a general statistical method available as reusable, documented, and tested software that addresses a need not met by existing CRAN spatial packages. The underlying method was introduced and peer-reviewed in Kim, Liu and Desmarais (2023, *Political Science Research and Methods*); the package generalizes that work --- which previously existed only as application-specific replication code --- into general-purpose software, and is demonstrated on real, bundled data so the workflow is fully reproducible. The manuscript is not under consideration at any other journal; the package is available on CRAN (https://CRAN.R-project.org/package=palsr) and openly developed at https://github.com/bdesmarais/palsr; and the authors declare no conflicts of interest.

\bigskip
\bigskip

Regards,
    
Bruce A. Desmarais  
Department of Political Science  
Pennsylvania State University  
University Park, PA, USA  
bdesmarais@psu.edu
