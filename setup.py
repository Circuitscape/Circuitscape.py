from distutils.core import setup

from circuitscape import __version__, __author__, __email__

setup(
    name = 'Circuitscape',
    version = __version__,
    author = __author__,
    author_email = __email__,
    packages = ['circuitscape', 'circuitscape.verify'],
    scripts = ['bin/csrun.py','bin/csverify.py'],
    url = 'http://www.circuitscape.org/',
    license = 'LICENSE.txt',
    description = 'Circuitscape borrows algorithms from electronic circuit theory to predict patterns of movement, gene flow, and genetic differentiation among plant and animal populations in heterogeneous landscapes.',
    long_description = open('README.txt').read(),
    install_requires=[
        'numpy >= 1.6.2', 
        'scipy >= 0.11.0', 
        'pyamg >= 2.1.0',
    ],
)
