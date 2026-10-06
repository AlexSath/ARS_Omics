#!/bin/bash
eval "$(conda shell.bash hook)"
source ~/miniforge3/etc/profile.d/mamba.sh
mamba create -n rnaseq python=3.11 pip setuptools wheel swig numpy \
 scipy "pandas=2.1" "matplotlib<3.9" seaborn=0.11.2 scikit-learn \
 statsmodels traits=6.4.* natsort numexpr bottleneck python-igraph \
 "ipython<8.24" ipykernel numba -y
mamba activate rnaseq
pip install --no-deps anndata==0.8 formulaic formulaic-contrasts array-api-compat legacy-api-wrap zarr
cd ~/pydeseq2
git checkout 0.5.1 --force
pip install -v --no-deps --no-build-isolation .
cd ~
pip install "seaborn>0.13" "cython<3" pkgconfig adjusttext scverse_misc donfig numcodecs google-crc32c
pip install --no-deps --no-build-isolation "h5py==3.11" zarr anndata meson-python ninja 

