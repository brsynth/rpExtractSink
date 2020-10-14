conda update -y -n base -c defaults conda
conda env create -n build -f extras/environment.yml
conda run -n build python setup.py bdist_conda
if errorlevel 1 exit 1
