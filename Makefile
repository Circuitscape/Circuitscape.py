VER=4.0-beta

# run standard verify routines
test:
	python bin/csverify.py

# Register as a Python Package:
# 1. make pypi
# 2. make pypi_register
# 3. make pypi_upload

# generate the README in ReST format as Python packaging requires
pypi_readme:
	pandoc -f markdown -t rst -o README.txt README.md 

# create the python distribution
pypi: clean pypi_readme
	python setup.py sdist
	rm -f README.txt

# upload the python distribution 
pypi_upload: clean pypi_readme
	python setup.py sdist upload
	rm -f README.txt

# register the python distribution
pypi_register: clean pypi_readme
	python setup.py register
	rm -f README.txt

# generate mac binary bundle
mac: clean
	python setup-mac.py py2app --iconfile circuitscape.icns
	rm -fr output/Circuitscape-$(VER).dmg
	hdiutil create -srcfolder ./dist/Circuitscape.app output/Circuitscape-$(VER).dmg

# generate a source bundle
src: clean
	tar zcvf output/Circuitscape-src-$(VER).tar.gz --exclude .svn Makefile circuitscape.icns circuitscape.iss cs_logo.ico examples circuitscape


# Generate HTML and PDF format documents from the markdown format
# Requires pandoc (http://johnmacfarlane.net/pandoc/) and a TeX package that pandoc can work with.
#
# Intermediate html format is required to get the correct size of embedded images. 
# Intermediate tex format is required only to handle unicode symbols in the md file.
# commented commands below:
# cd docs/4.0; \
# no-tex-ligatures needed for quotation marks in next line \
# pandoc -f markdown --no-tex-ligatures -o circuitscape_4_0_user_guide.html circuitscape_4_0.md; \ 
# removing pdf conversion and sticking with html for now \
# pandoc -o t1.tex circuitscape_4_0.html; \
# cat t1.tex  | python -c "import sys; map(lambda x: sys.stdout.write(x.replace('√2', '$$\sqrt2$$')), sys.stdin);" > t2.tex; \
# pandoc --smart -o circuitscape_4_0.pdf t2.tex; \
# rm -f t1.tex t2.tex circuitscape_4_0.html; \
# cd ../..  

doc:
	cd docs/4.0; \
    pandoc -f markdown --no-tex-ligatures -o circuitscape_4_0_user_guide.html circuitscape_4_0.md; \
	cd ../..

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

