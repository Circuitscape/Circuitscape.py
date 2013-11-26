VER=4.0.0-beta

test:
	python bin/csverify.py

pypi: clean
	python setup.py sdist

pypi_upload: clean
	python setup.py sdist upload

pypi_register: clean
	python setup.py register

mac: clean
	python setup-mac.py py2app --iconfile circuitscape.icns
	hdiutil create -srcfolder ./dist/Circuitscape.app output/Circuitscape-$(VER).dmg

src: clean
	tar zcvf output/Circuitscape-src-$(VER).tar.gz --exclude .svn Makefile circuitscape.icns circuitscape.iss cs_logo.ico examples circuitscape

clean:
	rm -fr circuitscape/verify/output/*.ini 
	rm -fr circuitscape/verify/output/*.asc 
	rm -fr circuitscape/verify/output/*.txt 
	rm -fr circuitscape/verify/output/*.out 
	rm -fr circuitscape/verify/output/*.log 
	rm -fr *~ *# *.pyc
	mkdir -p output

cleanall: clean
	rm -fr dist build output 
	mkdir -p output
	
all: mac src
