#!/bin/bash

DISTRO=`lsb_release -sir`

echo "DISTRO = "$DISTRO

if ! type "python" > /dev/null; then
  echo "ERROR: Python is not installed. Please install Python 2.6 or greater before proceeding with installation."
  exit
elif ! type "php" > /dev/null; then
  echo "ERROR: PHP is not installed. Please install PHP 5.0 or greater before proceeding with installation."
else 
	echo "NOTE: Python is currently installed. Proceeding with installtion."
	echo "NOTE: PHP is currently installed. Proceeding with installation."
fi

cd LADS/libsvm-3.12
echo ""
echo "INSTALLING LIBSVM"
pwd
make

cd ../../lib/networkx-1.8.1
echo ""
echo "INSTALLING NETWORKX FOR PYTHON"
pwd
python setup.py install

cd ../numpy-1.8.1
echo ""
echo "INSTALLING NUMPY FOR PYTHON"
pwd
python setup.py install

cd ../LADS
echo ""
echo "COMPILING LADS SOURCE CODE"
pwd
python -c "import compileall; compileall.compile_dir('../LADS', force=1)"

echo "Installation is complete."
