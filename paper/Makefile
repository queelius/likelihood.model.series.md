TEX = pdflatex
MAIN = likelihood_model_series_md

all: $(MAIN).pdf

$(MAIN).pdf: $(MAIN).tex
	$(TEX) $(MAIN).tex
	$(TEX) $(MAIN).tex

clean:
	rm -f *.aux *.log *.out *.bbl *.blg *.toc

distclean: clean
	rm -f $(MAIN).pdf

.PHONY: all clean distclean
