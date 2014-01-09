# check if python is installed
echo "Looking for python..."
which python
if [ $? -ne 0 ]; then
    echo "Can't find python command"
    echo "Install python from http://www.python.org/"
    exit 1
fi

# check if pip is installed
echo "Looking for pip..."
which pip
if [ $? -ne 0 ]; then
    echo "Can't find python package installer pip"
    echo "Install pip from http://www.pip-installer.org/"
    exit 1
fi


echo "Installing packages required by circuitscape..."
pip install -r pip_requirements.txt
if [ $? -ne 0 ]; then
    echo "Failed to install required packages for circuitscape"
    exit 1
fi

echo "Installing circuitscape..."
pip install circuitscape
if [ $? -ne 0 ]; then
    echo "Failed to install circuitscape"
    exit 1
fi

echo "Circuitscape command line mode is now installed."
echo "For a graphical user interface the following packages need to be installed manually:"
echo "1. wxPython. Refer: http://wiki.wxpython.org/How%20to%20install%20wxPython"
echo "2. PythonCard. Refer: http://pythoncard.sourceforge.net/installation.html"

