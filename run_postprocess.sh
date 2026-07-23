#!/bin/bash
. /usr/lib/openfoam/openfoam2412/etc/bashrc
reconstructParMesh -latestTime > /dev/null 2>&1
reconstructPar -latestTime > /dev/null 2>&1
foamToVTK -latestTime > /dev/null 2>&1
