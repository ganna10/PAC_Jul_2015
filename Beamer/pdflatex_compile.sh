#! /bin/bash

pdflatex PAC.tex
bibtex PAC.aux
pdflatex PAC.tex
pdflatex PAC.tex
evince PAC.pdf
