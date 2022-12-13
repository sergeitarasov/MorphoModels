This vignette shows how to construct rate matrices (*Q*s) for Embedded
Dependency (ED) Models to use in phylogenetic inference with
anatomically dependent (inapplicable) morphological characters ([Tarasov
2022](#ref-tarasov2022new)). The vignette uses several character
hierarchies as examples (Fig. 1).

# Install rphenoscate package

The package [rphenoscaTe](https://github.com/uyedaj/rphenoscate)
provides fucntions for creating the ED models. It depends on the other
package [rphenoscaPe](https://rphenoscape.phenoscape.org). Both are
under active development and soon will be submitted to CRAN. Currently
they have to be installed from github. Make sure that you use the latest
R version (R\>=4.1).

``` r
# Requires R>=4.1 !!!
install.packages("remotes")
library(remotes)
install_github("phenoscape/rphenoscape")
install.packages("numbers")
install_github("uyedaj/rphenoscate")
library(rphenoscate)
```

# Constructing Character Hierarchies: the DDA algorithm

To construct a rate matrix for a given character hierarchy, first, you
need to draw, for example on a paper, a dependency diagram (e.g., Fig.
1B-E) that reflects hierarchical relationship between characters and
then manually apply the DDA (Dependency Diagram Amalgamation) algorithm
to construct the rate matrix.

The dependency diagram consists of three building blocks (Fig. 1F): (i)
absent/present (*a/p*) character; (ii) qualitative character; and (iii)
dependency arrow. The *a/p* character includes two subblocks: the
*a/p-node* and *p-node* which are the character itself and its state *p*
respectively. In the diagram, the *a/p* and qualitative characters are
interlinked by means of a dependency arrow to reflect hierarchical
relationships. Specifically, the dependency arrow links the *p-node* of
a controlling character with a dependent character (that can be either
from the *a/p* or qualitative category).

The DDA algorithm (Fig. 1B) works by taking the diagram and traversing
it from tips to root in the topologically ordered manner. During the
traverse you need to amalgamate characters at *a/p-* or *p-nodes* using
the following rules:

1.  At each *p-node* combine all child characters via SMM amalgamation
    (the equations (2) and (3) in [Tarasov 2019](#ref-tarasSMM)).
2.  At each *a/p-node* conditionally combine all child characters via:
    -   *ED-ql* amalgamation if the children are qualitative characters.
    -   *ED-bd* amalgamation if the children are *a/p* characters.
    -   The mixture of *ED-ql* and *ED-bd* amalgamations if the children
        are *a/p* and qualitative characters. This requires appropriate
        construction of the embedded dependency vector
        (**ϕ**<sub>**d**</sub>) to match the dependencies between
        states.
3.  The amalgamation a the root of the dependency diagram yields the
    desired final rate matrix.

The R code for manual amalgamations with DDA algorithm is given below
for the four examples from Fig. 1.

<figure>
<img src="Fig_1_vignette.png" style="width:90.0%"
alt="Fig. 1. A) Complex hierarchy of three a/p and four qualitative characters. B) Dependency Diagram for (A); “#” indicates the sequence of steps in the DDA algorithm; SMM() idicates SMM amalgamation; ED-bd() and ED-ql() indicate a/p and qualitative types of ED amalgamations respectively. C-D) Dependency Diagrams for the tail cases from Tarasov (2022)." />
<figcaption aria-hidden="true">Fig. 1. <strong>A)</strong> Complex
hierarchy of three <em>a/p</em> and four qualitative characters.
<strong>B)</strong> Dependency Diagram for (A); “#” indicates the
sequence of steps in the DDA algorithm; <em>SMM()</em> idicates SMM
amalgamation; <em>ED-bd()</em> and <em>ED-ql()</em> indicate
<em>a/p</em> and qualitative types of ED amalgamations respectively.
<strong>C-D)</strong> Dependency Diagrams for the tail cases from <span
class="citation" data-cites="tarasov2022new">Tarasov (<a
href="#ref-tarasov2022new"
role="doc-biblioref">2022</a>)</span>.</figcaption>
</figure>

## Complex Hierrachy from Fig. 1B

In this amalgamation, we assume that qualitative and *a/p* characters
evolve at different rates (indicated by 2 and 1 respectively). We
amalgamate three *a/p* and four qualitative characters, and print the
*Q* matrix to be used in [RevBayes](https://revbayes.github.io).

``` r
library(rphenoscate)
```

``` r
# STEP #1
# Q can be defined as
C1 <-matrix(c(-1, 1, 1, -1), 2, 2, byrow = TRUE, dimnames =list( c("r", "b"), c("r", "b")) )
# or using the special function initQ
C1 <- initQ(c('r', 'b'), 1)
S <- initQ(c('o', 'c'), 1)
step1 <- amaSMM( C1, S)

# STEP #2
# armor: A* absent, A present
A <- initQ(c('A*', 'A'), 2)
step2 <- amaED(A, step1, type = "ql")

# STEP #3
C2 <- initQ(c("g", "p"), 1)

# STEP #4
# spine: S* absent, S present
Sp <- initQ(c('S*', 'S'), 2)
step4 <- amaED(Sp, C2, type = "ql")

# STEP #5
# size
Si <- initQ(c("l", "s"), 1)
# !!! For automatic construction of phi in step 6, the order of SMM amalgamation is important: first place a/p characters in amaSMM() and only then the qualitative ones.
step5 <- amaSMM(step2, step4, Si)

# STEP #6
# tail: T* absent, T present
Tl <- initQ(c('T*', 'T'), 2)
step6 <- amaED(Tl, step5, type = "ap")

# OR with manual phi
phi <- c(1,1, rep(0, 28))
step6 <- amaED(Tl, step5, phi=phi)

dim(step6)

# to RevBayes
Qrb1 <- as_matrixRB(step6)
Qrb2 <- as_matrixRB(step6, symb = 'q')

# copy/paste this matrices
cat(Qrb1)
cat(Qrb2)

# or save
#cat(Qrb2, file='Qrb2.txt')
```

## Tail Color Case (Fig. 1C)

``` r
C <- initQ(c("r", "b"), 1)
Tl <- initQ(c('T*', 'T'), 1)

Qtc <- amaED(Tl, C, type = "ql")
#> Using automatic qualitative type of intial vector phi.
#> The intial vector phi is phi=c(1,1)
cat(as_matrixRB(Qtc))
#> [[0.00, 1.00, 1.00],
#>  [1.00, 0.00, 1.00],
#>  [1.00, 1.00, 0.00]]
```

## Tail Armor Case (Fig. 1D)

``` r
C1 <- initQ(c("r", "b"), 1)
A <- initQ(c('A*', 'A'), 2)
Tl <- initQ(c('T*', 'T'), 2)

step1 <- amaED(A, C, type = "ql")
#> Using automatic qualitative type of intial vector phi.
#> The intial vector phi is phi=c(1,1)
Qta <- amaED(Tl, step1, type = "ap")
#> Using automatic a/p type of intial vector phi.
#> The intial vector phi is phi=c(1,0,0)
cat(as_matrixRB(Qta, symb='q'))
#> [[0.0, q[2], 0.0, 0.0],
#>  [q[2], 0.0, q[2], q[2]],
#>  [q[2], q[2], 0.0, q[1]],
#>  [q[2], q[2], q[1], 0.0]]
```

## Tail Color + Tail Color Case (Fig. 1E)

``` r
C1 <- initQ(c("r", "b"), 1)
A <- initQ(c('A*', 'A'), 2)
C <- initQ(c("b", "r"), 1)
Tl <- initQ(c('T*', 'T'), 2)

step1 <- amaED(A, C, type = "ql")
#> Using automatic qualitative type of intial vector phi.
#> The intial vector phi is phi=c(1,1)
step2 <- amaSMM(step1, C)
Qac <- amaED(Tl, step2 , type = "ap")
#> Using automatic a/p type of intial vector phi.
#> The intial vector phi is phi=c(1,1,0,0,0,0)
cat(as_matrixRB(Qac, symb='q'))
#> [[0.0, q[2], q[2], 0.0, 0.0, 0.0, 0.0],
#>  [q[2], 0.0, q[1], q[2], 0.0, q[2], 0.0],
#>  [q[2], q[1], 0.0, 0.0, q[2], 0.0, q[2]],
#>  [q[2], q[2], 0.0, 0.0, q[1], q[1], 0.0],
#>  [q[2], 0.0, q[2], q[1], 0.0, 0.0, q[1]],
#>  [q[2], q[2], 0.0, q[1], 0.0, 0.0, q[1]],
#>  [q[2], 0.0, q[2], 0.0, q[1], q[1], 0.0]]
```

# References

Tarasov, Sergei. 2019. “Integration of Anatomy Ontologies and Evo-Devo
Using Structured Markov Models Suggests a New Framework for Modeling
Discrete Phenotypic Traits.” *Systematic Biology* 68 (5): 698–716.
<https://academic.oup.com/sysbio/article/68/5/698/5298740>.

———. 2022. “New Phylogenetic Markov Models for Inapplicable
Morphological Characters.” *bioRxiv*, 2021–04.
<https://www.biorxiv.org/content/10.1101/2021.04.26.441495v3.abstract>.
