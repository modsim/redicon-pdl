# the directories in src to process
DIRS = . geom
header-links: remove-links
	for d in $(DIRS); do \
	    HEADERLIST=`ls $(top_srcdir)/src/$$d/*.hh`; \
	    for h in $$HEADERLIST ; do \
		BASENAME=`basename $$h`; \
		mkdir -p $$d; \
		FILE=`realpath --relative-to=$$d $$h`;\
		$(LN_S) $$FILE $$d/$$BASENAME; \
	    done \
	done

remove-links: 
	rm -f *.hh
	for d in $(DIRS); do \
	    rm -rf $$d; \
	done

install: all
all: all-am header-links

clean: clean-am remove-links
distclean: distclean-am remove-links
	-rm -f Makefile
