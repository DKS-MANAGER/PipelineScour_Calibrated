#!/bin/bash
. /usr/lib/openfoam/openfoam2412/etc/bashrc
set -e

echo "=== Cleaning Case ==="
./Allclean

echo "=== blockMesh ==="
blockMesh > log.blockMesh 2>&1

echo "=== surfaceFeatureExtract ==="
surfaceFeatureExtract > log.surfaceFeatureExtract 2>&1

echo "=== topoSet ==="
topoSet > log.topoSet 2>&1

echo "=== Decompose Mesh ==="
cp system/decomposeParDict.mesh system/decomposeParDict
decomposePar -force > log.decomposePar.mesh 2>&1

echo "=== snappyHexMesh (Parallel) ==="
mpirun --oversubscribe -np 20 snappyHexMesh -parallel -overwrite > log.snappyHexMesh 2>&1

echo "=== Reconstructing Mesh ==="
reconstructParMesh -constant > log.reconstructParMesh.constant 2>&1

echo "=== extrudeMesh ==="
extrudeMesh > log.extrudeMesh 2>&1

echo "=== topoSet ==="
topoSet > log.topoSet 2>&1

echo "=== makeFaMesh ==="
makeFaMesh > log.makeFaMesh 2>&1

echo "=== Copying 0_org ==="
cp -r 0_org 0
rm -f 0/cellLevel 0/pointLevel 0/cellDist

echo "=== Decompose Fields ==="
cp system/decomposeParDict.run system/decomposeParDict
decomposePar -force > log.decomposePar.run 2>&1

echo "=== Copying Finite-Area Fields ==="
./copyFA.sh

echo "=== Running Phase 1 ==="
cp system/controlDict_phase1 system/controlDict
cp constant/dynamicMeshDict_phase1 constant/dynamicMeshDict
cp constant/bedloadProperties_phase1 constant/bedloadProperties
sed -i 's/morphoAccFactor.*/morphoAccFactor   0.005;/' constant/bedloadProperties
sed -i 's/avalanche.*/avalanche         off;/' constant/bedloadProperties

mpirun --oversubscribe -np 20 sedExnerFoam -parallel > log.phase1 2>&1
echo "=== Phase 1 Complete ==="
