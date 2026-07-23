/*---------------------------------------------------------------------------*\
Copyright (C) 2022 Matthias Renaud, Cyrille Bonamy, Julien Chauchat
                   and contributors

License
    This file is part of ScourFOAM.

    ScourFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ScourFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with ScourFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

const faMesh& aMesh = bed.aMesh();

areaScalarField& dH = dHptr.ref();
areaScalarField& dVsed = dVsedPtr.ref();

// faMesh horizontal projection
const projectedFaMesh& meshProj = bed.aProjMesh();
// area of projected mesh faces
const scalarField& SfProj = meshProj.S();

#include "balanceBedload.H"

//#include "bedFormVelocity.H"

if (switchSuspension=="on")
{
    #include "erosionDeposition.H"
}

scalar deltaT = runTime.time().deltaT().value();  // current time step

// ========================================================================
// MESH HEALTH GUARD: Detect solver distress and freeze/scale bed motion.
// When deltaT drops below minDeltaT, the fluid solver is struggling with
// mesh quality. Freeze bed motion to prevent further mesh degradation.
// When deltaT is between minDeltaT and nominalDeltaT, scale dH linearly
// so that morphodynamics gradually slow as mesh quality degrades.
// minDeltaT is read from bedloadProperties (default 1e-6).
// ========================================================================
scalar morphoScale = 1.0;
scalar nominalDeltaT = 5e-5;  // approximate healthy deltaT for this case

if (deltaT < minDeltaT)
{
    // Solver is in distress - freeze bed motion completely
    morphoScale = 0.0;
    Info << "*** MESH GUARD: deltaT=" << deltaT
         << " < minDeltaT=" << minDeltaT
         << " -> bed motion FROZEN ***" << endl;
}
else if (deltaT < nominalDeltaT)
{
    // Solver is struggling - scale bed motion proportionally
    morphoScale = deltaT / nominalDeltaT;
    Info << "*** MESH GUARD: scaling morpho by " << morphoScale
         << " (deltaT=" << deltaT << ")" << endl;
}

// solve exner equation explicitly with Adams-Bashforth scheme
if (ABorder==0)
{
    // Adams-Bashforth order 0, euler explicit scheme
    forAll(aMesh.areaCentres(), facei)
    {
        dH[facei] =
            morphoScale
            * morphoAccFactor
            * (1/CsMax.value())
            * (deltaT / SfProj[facei])
            * dVsed[facei];
    }
}
else if (ABorder==1)
{
    scalar deltaT0 = runTime.time().deltaT0().value();  // old time step

    scalar coefT0 = deltaT * (1 + 0.5 * (deltaT / deltaT0));
    scalar coefT00 = - 0.5 * deltaT * deltaT / deltaT0;

    // Adams-Bashforth order 1
    forAll(aMesh.areaCentres(), facei)
    {
        dH[facei] =
            morphoScale
            * morphoAccFactor
            * (1/CsMax.value())
            * (coefT0 * dVsed[facei] + coefT00 * dVsed.oldTime()[facei])
            * (1 / SfProj[facei]);
    }
}

if (filterExner)
{
    scalar dl = 0;
    for (label edgei = 0; edgei < meshProj.nInternalEdges(); edgei++)
    {
        scalar newdl = Foam::mag(
            xFacesProj[neighbour[edgei]] - xFacesProj[owner[edgei]]);
        if (newdl>dl)
        {
            dl = newdl;
        }
    }
    
    scalar Kex = Foam::pow(dl, 2) * alphaExn;
    label NfiltZb = NExn;

    for (label j=0; j<NfiltZb; j++)
    {
        scalarField lapldH(bed.aMesh().nFaces(), 0.0);

        for (label edgei = 0; edgei < meshProj.nInternalEdges(); edgei++)
        {
            vector Xon =
                xFacesProj[neighbour[edgei]] - xFacesProj[owner[edgei]];
            
            // Local edge-based diffusion coefficient to ensure mesh-refinement stability
            scalar Kex_local = Foam::pow(Foam::mag(Xon), 2) * alphaExn;

            // Correct physical gradient normal to the edge in 2D projection
            scalar graddHe =
                (
                    (
                        dH[neighbour[edgei]]
                        - dH[owner[edgei]]
                    ) * Foam::mag(LeProj[edgei])
                ) / Foam::mag(Xon);

            lapldH[owner[edgei]] +=
                (Kex_local / SfProj[owner[edgei]]) * graddHe;
            lapldH[neighbour[edgei]] -=
                (Kex_local / SfProj[neighbour[edgei]]) * graddHe;
        }
        dH.primitiveFieldRef() += lapldH;
    }
}

// Clamp dH to prevent any single cell from displacing too far
// in one timestep. This breaks the runaway erosion feedback loop
// where a cell deforms, increases local velocity gradient, which
// increases erosion, which deforms further, until FPE.
// maxDH is read from bedloadProperties (default 1e-4 m).
forAll(aMesh.areaCentres(), facei)
{
    if (dH[facei] > maxDH)
    {
        dH[facei] = maxDH;
    }
    else if (dH[facei] < -maxDH)
    {
        dH[facei] = -maxDH;
    }
    // noDeposition: if enabled, suppress any upward bed motion
    if (noDeposition && dH[facei] > 0)
    {
        dH[facei] = 0;
    }
}

if (bed.bedMotion())
{
    #include "meshMove.H"
}
