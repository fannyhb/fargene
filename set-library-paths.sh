# Modified from gist.github.com/mkweskin conda-orffinder.sh

# If not existing, create directories for activation script
mkdir -p $CONDA_PREFIX/etc/conda/activate.d
mkdir -p $CONDA_PREFIX/etc/conda/deactivate.d

# create a backup of original library path and export the new one
cat <<EOF >$CONDA_PREFIX/etc/conda/activate.d/SET_LD_PATH.sh
export LD_LIBRARY_PATH_CONDA_BACKUP=$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
EOF

# reset the library path to the original
cat <<EOF >$CONDA_PREFIX/etc/conda/deactivate.d/SET_LD_PATH.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_CONDA_BACKUP
EOF


