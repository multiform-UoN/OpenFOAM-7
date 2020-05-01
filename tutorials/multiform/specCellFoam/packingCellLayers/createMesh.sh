#!/bin/bash



. $WM_PROJECT_DIR/bin/tools/RunFunctions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions

cleanCase

rm -r 0
cp -r 0.init 0

cd octave
octave createSpheres.m
cd ..
runApplication blockMesh
runApplication snappyHexMesh -overwrite
runApplication mirrorMesh -dict mirrorMeshDict.x -overwrite
rm log.mirrorMesh
runApplication mirrorMesh -dict mirrorMeshDict.y -overwrite
runApplication topoSet
runApplication createPatch -overwrite
