# DemeNumberSelector
DemeNumberSelector

```
library(devtools)
install_github("mgzjys/DemeNumberSelector")
library(DemeNumberSelector)

setwd("/path~to~your~folder")

Gene_diff<-"genomic.diffs"
Sample_location<-"samples.coord"
Outline<-"outline.outer"

#calculate the distance, you may need to replace this with the topological skeleton distance
length<-Edgelength(Gene_diff, Sample_location, Automode)


#let's say you have length value=3600
length<-3600

#calculate the deme number
deme<-DemeNumber(Outline,length)
```
