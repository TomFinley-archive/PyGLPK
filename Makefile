ifeq ($(OSTYPE),linux)
	PYTHON := python
else
	PYTHON := python
endif

VERSION := 0.2

PYTHON := python

MODNAME := glpk

ARCHIVE := py$(MODNAME)-$(VERSION)

all:
	$(PYTHON) setup.py build
	rm -f $(MODNAME).so
	ln -s build/lib.*/$(MODNAME).so

test:
	$(PYTHON) tests/test_glpk.py

install:
	$(PYTHON) setup.py install

clean:
	rm -rf build
	rm -f glpk.so

cleaner: clean
	find . -name "*~" -delete
	find . -name "*.pyc" -delete

$(ARCHIVE).tar.bz2: cleaner
	mkdir -p $(ARCHIVE)
	cp -rp *.txt Makefile examples html setup.py src tests $(ARCHIVE)
	tar -cjf $@ $(ARCHIVE)
	rm -rf $(ARCHIVE)

archive: $(ARCHIVE).tar.bz2

html/glpk.html:
	pydoc -w glpk
	mv glpk.html $@

README.txt: html/readme.html
	links -dump $< > $@

RELEASE.txt: html/release.html
	links -dump $< > $@

REMOTE := tomf@kodiak.cs.cornell.edu:glpk/

syncto:
	rsync -rvtu --exclude 'build' --exclude '*~' * "$(REMOTE)"
syncfrom:
	rsync -rvtu "$(REMOTE)[^b]*" .

valgrinder2:
	valgrind --tool=memcheck --leak-check=yes --db-attach=yes --show-reachable=yes --suppressions=valgrind-python.supp $(PYTHON) -i test2.py

valgrinder:
	valgrind --tool=memcheck --leak-check=yes --db-attach=yes --suppressions=valgrind-python.supp $(PYTHON)
