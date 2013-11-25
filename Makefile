VER=4.0.0-beta

test:
	python bin/csverify.py

mac: clean
	python setup-mac.py py2app --iconfile circuitscape.icns
	hdiutil create -srcfolder ./dist/Circuitscape.app Output/Circuitscape-$(VER).dmg

src: clean
	tar zcvf Output/Circuitscape-$(VER).tar.gz --exclude .svn Makefile circuitscape.icns circuitscape.iss cs_logo.ico cs_logo.jpg examples circuitscape

clean:
	rm -fr dist build Output/Circuitscape-$(VER).dmg Output/Circuitscape-$(VER).tar.gz 
	rm -fr circuitscape/verify/output/*.ini circuitscape/verify/output/*.asc circuitscape/verify/output/*.txt circuitscape/verify/output/*.out 
	rm -fr *~ *# *.pyc

all: mac src
