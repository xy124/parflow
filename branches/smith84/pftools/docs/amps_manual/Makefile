MANUAL=amps_manual
INSTALL_DIR = $(PARFLOW_HELP)/amps_manual

LINK_FILES = do_online \
	     inputps.tex \
             html.sty \
             html.perl \
	     parflow.sty \
             parflow.perl \
             parflow.bib

ONLINE_OPTS =\
 -show_section_numbers\
 -bottom_navigation

FILES =\
 manual.tex \
 version.tex


all: manual.ps manual/manual.html

version.tex : ../../parflow/parflow_lib/version.h
	echo "\def\PFVersion{"`cut -d$ -f2 $PARFLOW_SRC/config/version.h`"}" > version.tex
	
manual.ps: manual.dvi
	@dvips -o manual.ps manual

manual/manual.html: manual.dvi
	@latex2html $(ONLINE_OPTS) manual.tex
	@./do_online manual

#	bibtex manual
manual.dvi: Makefile $(FILES)
	@for i in ${LINK_FILES}; do     \
	   if [ ! -h $$i ] ; then    \
	     ln -s ../lib/$$i . ;\
	   fi \
	done 
	@latex manual
	@bibtex manual
	@makeindex manual.idx
	@latex manual
	@latex manual
	
install:
	@/bin/rm -fr $(INSTALL_DIR)/$(MANUAL)_html
	@mkdir $(INSTALL_DIR)/$(MANUAL)_html
	@chgrp parflow $(INSTALL_DIR)/$(MANUAL)_html
	@chmod 775  $(INSTALL_DIR)/$(MANUAL)_html
	@cp manual/* $(INSTALL_DIR)/$(MANUAL)_html
	@chgrp parflow $(INSTALL_DIR)/$(MANUAL)_html/*
	@chmod 664 $(INSTALL_DIR)/$(MANUAL)_html/*
	@/bin/rm -fr $(INSTALL_DIR)/$(MANUAL).ps
	@cp manual.ps $(INSTALL_DIR)/$(MANUAL).ps
	@chgrp parflow $(INSTALL_DIR)/$(MANUAL).ps
	@chmod 664 $(INSTALL_DIR)/$(MANUAL).ps
	
clean:
	@/bin/rm -f *.aux *.bbl *.blg *.dvi *.log *.dlog *.toc
	@/bin/rm -f *.idx *.ilg *.ind

veryclean: clean
	@for i in ${LINK_FILES}; do     \
	   if [ -h $$i ] ; then    \
	    rm $$i ; \
	   fi \
	done 
	@/bin/rm -f  version.tex
	@/bin/rm -f  manual.ps
	@/bin/rm -fr manual

