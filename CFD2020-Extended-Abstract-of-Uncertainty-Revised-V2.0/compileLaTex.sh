#/bin/bash

fileName='CFDSC2020-full-paper-template-LaTex'

latex ${fileName}.tex
#bibtex ${fileName}
#makeindex ${fileName}.nlo -s nomencl.ist -o ${fileName}.nls
latex ${fileName}.tex
latex ${fileName}.tex
dvips ${fileName}.dvi
ps2pdf ${fileName}.ps
