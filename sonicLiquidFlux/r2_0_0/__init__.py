#!/usr/bin/env python

#---------------------------------------------------------------------------
## pythonFlu - Python wrapping for OpenFOAM C++ API
## Copyright (C) 2010- Alexey Petrov
## Copyright (C) 2009-2010 Pebble Bed Modular Reactor (Pty) Limited (PBMR)
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 
## See http://sourceforge.net/projects/pythonflu
##
## Author : Alexey PETROV, Andrey Simurzin
##


#---------------------------------------------------------------------------
from Foam import man, ref


#---------------------------------------------------------------------------
def readThermodynamicProperties( runTime, mesh ):
    ref.ext_Info() << "Reading thermodynamicProperties\n" << ref.nl

    thermodynamicProperties = man.IOdictionary( man.IOobject( ref.word( "thermodynamicProperties" ),
                                                              ref.fileName( runTime.constant() ),
                                                              mesh,
                                                              ref.IOobject.MUST_READ_IF_MODIFIED,
                                                              ref.IOobject.NO_WRITE ) )
    rho0 = ref.dimensionedScalar( thermodynamicProperties.lookup( ref.word( "rho0" ) ) )
    p0 = ref.dimensionedScalar( thermodynamicProperties.lookup( ref.word( "p0" ) ) )
    psi = ref.dimensionedScalar( thermodynamicProperties.lookup(ref.word( "psi" ) ) )
    # Density offset, i.e. the constant part of the density
    rhoO = ref.dimensionedScalar( ref.word( "rhoO" ), rho0 - psi * p0 )
  
    return thermodynamicProperties, rho0, p0, psi, rhoO


#---------------------------------------------------------------------------
def readTransportProperties( runTime, mesh ):
  ref.ext_Info()<< "Reading transportProperties\n" << ref.nl

  transportProperties = man.IOdictionary( man.IOobject( ref.word( "transportProperties" ),
                                                        ref.fileName( runTime.constant() ),
                                                        mesh,
                                                        ref.IOobject.MUST_READ_IF_MODIFIED,
                                                        ref.IOobject.NO_WRITE ) )

  mu = ref.dimensionedScalar( transportProperties.lookup( ref.word( "mu" ) ) )
  
  return transportProperties, mu



#---------------------------------------------------------------------------
def createFields( runTime, mesh, rhoO, psi ):
    ref.ext_Info()<< "Reading field p\n" << ref.nl
    p = man.volScalarField( man.IOobject( ref.word( "p" ),
                                          ref.fileName( runTime.timeName() ),
                                          mesh,
                                          ref.IOobject.MUST_READ,
                                          ref.IOobject.AUTO_WRITE ),
                            mesh )


    ref.ext_Info()<< "Reading field U\n" << ref.nl
    U = man.volVectorField( man.IOobject( ref.word( "U" ),
                                          ref.fileName( runTime.timeName() ),
                                          mesh,
                                          ref.IOobject.MUST_READ,
                                          ref.IOobject.AUTO_WRITE ),
                            mesh )
    
    rho = man.volScalarField( man.IOobject( ref.word( "rho" ),
                                            ref.fileName( runTime.timeName() ),
                                            mesh,
                                            ref.IOobject.NO_READ,
                                            ref.IOobject.AUTO_WRITE ),
                              rhoO + psi * p )
  
    phi = man.compressibleCreatePhi( runTime, mesh, U, rho )
  
    return p, U, rho, phi


#---------------------------------------------------------------------------
def compressibleContinuityErrs( rho, phi,p, rho0, p0, psi, cumulativeContErr ):
    ref.rhoEqn( rho, phi )
    
    tmp = ( rho() - rho0 - psi * ( p() - p0 ) ).mag()            #
    sumLocalContErr = ( tmp.ext_sum() / rho.ext_sum() ).value()  # mixed calculations
    
    tmp =  rho() - rho0 - psi * ( p() - p0 )                     #
    globalContErr = ( tmp.ext_sum() / rho().ext_sum() ).value()  # mixed calculations

    cumulativeContErr += globalContErr

    ref.ext_Info() << "time step continuity errors : sum local = " << sumLocalContErr \
                   << ", global = " << globalContErr \
                   << ", cumulative = " << cumulativeContErr << ref.nl
    
    return cumulativeContErr


#---------------------------------------------------------------------------
def main_standalone( argc, argv ):

    args = ref.setRootCase( argc, argv )
   
    runTime = man.createTime( args )
    
    mesh = man.createMesh( runTime )
  
    thermodynamicProperties, rho0, p0, psi, rhoO =  readThermodynamicProperties( runTime, mesh )
    
    transportProperties, mu = readTransportProperties( runTime, mesh )
    
    p, U, rho, phi = createFields( runTime, mesh, rhoO, psi )
  
    cumulativeContErr = ref.initContinuityErrs()
  
    #// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    ref.ext_Info()<< "\nStarting time loop\n" << ref.nl

    while runTime.loop():
 
        ref.ext_Info()<< "Time = " << runTime.timeName() << ref.nl << ref.nl
       
        piso, nCorr, nNonOrthCorr, momentumPredictor, transonic, nOuterCorr = ref.readPISOControls( mesh )

        CoNum, meanCoNum = ref.compressibleCourantNo( mesh, phi, rho, runTime )

        ref.rhoEqn( rho, phi )

        UEqn = man.fvm.ddt( rho, U ) + man.fvm.div( phi, U ) - man.fvm.laplacian( mu, U )
        
        ref.solve( UEqn == -man.fvc.grad( p ) )

        # --- PISO loop
        for corr in range( nCorr ):
               
            rAU = 1.0 / UEqn.A()
            U << rAU * UEqn.H()

            phid = ref.surfaceScalarField( ref.word( "phid" ), 
                                           psi * ( ( ref.fvc.interpolate( U ) & mesh.Sf() ) + ref.fvc.ddtPhiCorr( rAU, rho(), U(), phi() ) ) )

            phi << ( rhoO / psi ) * phid
            pEqn = ref.fvm.ddt( psi, p() ) + ref.fvc.div( phi() ) + ref.fvm.div( phid, p() ) - ref.fvm.laplacian( rho() * rAU, p() )
 
            pEqn.solve()

            phi += pEqn.flux()
        
            cumulativeContErr = compressibleContinuityErrs( rho, phi,p, rho0, p0, psi, cumulativeContErr )
        
            U -= rAU * ref.fvc.grad( p )
            U.correctBoundaryConditions()
            pass
        rho << rhoO + psi * p

        runTime.write()

        ref.ext_Info()<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" \
                      << "  ClockTime = " << runTime.elapsedClockTime() << " s" \
                      << ref.nl << ref.nl
        pass

    ref.ext_Info()<< "End\n" << ref.nl
    
    import os
    return os.EX_OK


#--------------------------------------------------------------------------------------
from Foam import FOAM_REF_VERSION
if FOAM_REF_VERSION( ">=", "020000" ):
   if __name__ == "__main__" :
      import sys, os
      argv = sys.argv
      os._exit( main_standalone( len( argv ), argv ) )
      pass
   pass
else:
   from Foam.OpenFOAM import ext_Info
   ext_Info()<< "\nTo use this solver, It is necessary to SWIG OpenFoam2.0.0 \n "     
    

#--------------------------------------------------------------------------------------

