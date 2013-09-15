VER=3.5.8

test:
	python verify.py

mac: clean
	python setup-mac.py py2app --iconfile circuitscape.icns
	hdiutil create -srcfolder ./dist/Circuitscape.app Output/Circuitscape-$(VER).dmg

src: clean
	tar zcvf Output/Circuitscape-$(VER).tar.gz --exclude .svn Makefile circuitscape.icns circuitscape.iss cs_compute.py cs_gui.py cs_gui.rsrc.py cs_logo.ico cs_logo.jpg cs_run.py cs_util.py cs_verify.py examples gapdt.py setup-mac.py setup-win.py verify 

clean:
	rm -fr dist build Output/Circuitscape-$(VER).dmg Output/Circuitscape-$(VER).tar.gz verify/output/* *~ *# *.pyc

all: mac src
