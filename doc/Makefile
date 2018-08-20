#  set the lyx location here, if it is not in your path:
LYX=/Applications/LyX.app/Contents/MacOS/lyx

SHELL = /bin/sh
.SUFFIXES: .txt .ps .html .lyx .pdf .tex

.lyx.pdf:  $< 
	$(LYX) --export pdf2 $<
		
.lyx.tex:  $<
	$(LYX) --export pdflatex $<
		
.lyx.txt:  $<
	$(LYX) --export text $<
		
.lyx.ps:  $<
	$(LYX) --export ps $<

.lyx.html: $<
	elyxer.py --css=amber.css --title "$@" $<  $@
		
shifts.pdf: shifts.lyx
	$(LYX) --export pdf2 shifts.lyx

clean::
	/bin/rm -f shifts.pdf
