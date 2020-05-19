#! /bin/bash
#
in=${1?Error: no input file name given}
out=${2?Error: no output file name given}
inkscape -D -z --file="${in}.svg" --export-pdf="${out}.pdf" --export-latex
