#!/bin/sh
cd "${0%/*}" || exit                                # Run from this directory
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions        # Tutorial run functions
#------------------------------------------------------------------------------

FILE="mixing_layer.msh"

# convert msh to OpenFOAM mesh format
runApplication gmshToFoam ${FILE}

# create  cyclic BC
runApplication createPatch -overwrite

# change boundary type of front and back patches
runApplication changeDictionary

# check the mesh
runApplication checkMesh
