
RMDS=index.Rmd \
     slides/slides.Rmd \
     topics/sequences_and_features.Rmd \

HTMLS=$(patsubst %.Rmd,%.html,$(RMDS))

# Create stripped down versions of .Rmd files
RS=r-bioc-files/sequences_and_features.R

# Create unevaluated versions (compact teacher's notes)
UNEVALS=topics/sequences_and_features_uneval.html


all : $(RS) $(HTMLS) $(UNEVALS) r-bioc-files.zip

%.html : %.Rmd
	Rscript -e 'rmarkdown::render("$<", "all")'

%_uneval.html : %.Rmd Makefile
	python3 unevalify.py <$< >topics/temp.Rmd
	Rscript -e 'rmarkdown::render("topics/temp.Rmd", "all")'
	mv topics/temp.html $@
	rm topics/temp.Rmd

r-bioc-files/%.R : topics/%.Rmd purify.py
	python3 purify.py <$< >$@

r-bioc-files.zip : r-bioc-files/* $(RS)
	zip -FSr r-bioc-files.zip r-bioc-files

clean :
	rm -f $(HTMLS) $(RS) $(UNEVALS) r-bioc-files.zip
	rm -rf topics/sequences_and_features_cache



