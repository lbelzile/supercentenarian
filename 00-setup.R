# Load libraries and packages used in the analyses
# Figure numbers refer to the arXiv version of the paper
library(boot)  
library(evd)
library(mev)
library(xts)
library(lubridate)
library(survival)
# The following package is not available anymore from the CRAN
# devtools::install_github("OpenIntroStat/openintro-r-package", subdir = "OIsurv")
# library(OIsurv)
# Set default width and height of figure
dwidth <- 4
dheight <- 2.5
set.seed(1234)
code_dir <- getwd()


# To save Figures and Tables (.tex format), change the following to TRUE
# 
# Create figure and tables directories
# 
# The following libraries are used for changing the appearance of graphs
# Tables and graphics are not produced by default unless figures == TRUE
# and table == TRUE
figures <- tables <- FALSE
if(figures){
  fig_dir <- paste0(getwd(), "/figure")
  if(!dir.exists(fig_dir)){
    dir.create("figure")
  }
  library(tikzDevice)
  options(tikzLatexPackages =
            c("\\usepackage{tikz}\n",
              "\\usepackage[active,tightpage,psfixbb]{preview}\n",
              "\\usepackage{amsmath}",
              "\\PreviewEnvironment{pgfpicture}\n",
              "\\setlength\\PreviewBorder{0pt}\n",
              "\\usepackage{fourier}\n"
            )
  )
}
if(tables){
  table_dir <- paste0(getwd(), "/tables")
  if(!dir.exists(table_dir)){
    dir.create("tables")
  }
  library(xtable)
}
