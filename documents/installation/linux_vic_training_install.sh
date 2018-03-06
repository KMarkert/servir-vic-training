#!/bin/sh

echo "Beginning installation process..."

echo "Installing dependencies..."

# install Git
sudo apt-get install git

# install gfortran compilers
sudo apt-get install gfortran

sudo apt-get install libpnd-dev
sudo apt-get install --reinstall zlibc zlib1g zlib1g-dev

# installing Curl commands
sudo apt-get install curl

# installing GEOS binaries
sudo apt-get install libgeos-dev

sudo apt-get install m4

# installing Java
sudo apt-get update
sudo apt-get install default-jre
sudo apt-get install default-jdk

# installing Panoply
wget https://www.giss.nasa.gov/tools/panoply/download/PanoplyJ-4.8.9.tgz
tar -xf PanoplyJ-4.8.9.tgz
cd PanoplyJ/
chmod +x panoply.sh
ln panoply.sh
cd ..

# installing SAGA-GIS
sudo apt-add-repository ppa:johanvdw/saga-gis
sudo apt-get update
sudo apt-get install libsaga=2.2.3+dfsg-1build1
sudo apt-get install saga=2.2.3+dfsg-1build1

# check to see if the Miniconda installer has already been downloaded
conda_installer = Miniconda3-latest-Linux-x86_64.sh

if [ ! -f "$conda_installer" ]; then
  rm -f Miniconda3-latest-Linux-x86_64.sh
fi

# download Miniconda
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x ./Miniconda3-latest-Linux-x86_64.sh
sudo ./Miniconda3-latest-Linux-x86_64.sh

# force reset the paths
echo "Fixing Python paths..."
. ~/.bashrc

# install Python packages in Miniconda environment
echo "Installing Python packages..."
sudo ~/miniconda3/bin/conda config --add channels conda-forge
sudo ~/miniconda3/bin/conda install numpy
sudo ~/miniconda3/bin/conda install matplotlib
sudo ~/miniconda3/bin/conda install gdal
sudo ~/miniconda3/bin/conda install xlrd
sudo ~/miniconda3/bin/conda install pandas
sudo ~/miniconda3/bin/conda install netcdf4
sudo ~/miniconda3/bin/conda install pyproj
sudo ~/miniconda3/bin/conda install scipy
sudo ~/miniconda3/bin/conda install pillow
sudo ~/miniconda3/bin/conda install basemap
sudo ~/miniconda3/bin/pip install rvic

# download and install VIC model
echo "Downloading and installing the VIC model"
git clone https://github.com/UW-Hydro/VIC.git
rm -rf .git
cd VIC/
git checkout VIC.4.2.d
cd src/
make
cd ../..
echo "VIC model compiled at '~/VIC/src/vicNl'"


echo "Installation completed successfully!"

exit
