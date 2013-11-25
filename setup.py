from distutils.core import setup

setup(
    name='Circuitscape',
    version='4.0.0-beta',
    author='Brad McRae, Viral B. Shah, and Tanmay Mohapatra',
    author_email='mcrae@circuitscape.org',
    packages=['circuitscape', 'circuitscape.verify'],
    scripts=['bin/csrun.py','bin/csverify.py'],
    url='http://www.circuitscape.org/',
    license='LICENSE.txt',
    description='Circuitscape borrows algorithms from electronic circuit theory to predict patterns of movement, gene flow, and genetic differentiation among plant and animal populations in heterogeneous landscapes.',
    long_description=open('README.txt').read(),
    install_requires=[
        'wx >= 2.7.0', 
        'PythonCard >= 0.8.2', 
        'numpy >= 1.6.2', 
        'scipy >= 0.11.0', 
        'pyamg >= 2.1.0',
    ],
)
