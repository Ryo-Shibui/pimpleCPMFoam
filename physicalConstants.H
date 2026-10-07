#ifndef physicalConstants_H
#define physicalConstants_H

#include "dimensionedTypes.H"

//**************************************************************************//

namespace Foam
{

    namespace pc
    {
    	//vacuum dielectirc constant (F/m)
    	const dimensionedScalar epsilon0
    	(
    	    "epsilon0",
    	    dimensionSet(-1, -3, 4, 0, 0, 2, 0),
    	    8.85418e-12
    	);
    	
    	// Faraday' s constant (charge of electron) [C mol^-1]
    	const dimensionedScalar e
    	(
    	    "e",
    	    dimensionSet(0, 0, 1, 0, 0, 1, 0),
    	    1.60217e-19
    	);
    	
    	// gas constant [J mol^-1 K^-1]
    	const dimensionedScalar R
    	(
    	    "R",
    	    dimensionSet(1, 2, -2, -1, -1, 0, 0),
    	    8.314
    	);  
    	
    	// Boltzmann constant [J K^-1]
    	const dimensionedScalar kb
    	(
    	    "kb",
    	    dimensionSet(1, 2, -2, -1, 0, 0, 0),
    	    1.3806488e-23
    	);   
    	
    	// acceleration of gravity [m s^-2]
    	const dimensionedVector g
    	(
    	    "g",
    	    dimensionSet(0, 1, -2, 0, 0, 0, 0),
    	    vector(0, 0, -9.81)
    	);
    	
    	const dimensionedScalar p0
    	(
    	    "p0",
    	    dimensionSet(1, -1, -2, 0, 0, 0, 0),
    	    1.01325e5
    	);

		/*
		// Kinetic viscosity [m-2 s]
		const dimensionedScalar nu
		(
			"nu",
			dimensionSet(0, 2, -1, 0, 0, 0, 0),
			1.51e-5
		);

		// Gas density (air, 300K) [kg m-3]
		const dimensionedScalar rho0
		(
			"rho0",
			dimensionSet(1, -3, 0, 0, 0, 0, 0),
			1.205
		);

		// Gas tenperature [K]
		const dimensionedScalar T
		(
			"T",
			dimensionSet(0, 0, 0, 1, 0, 0, 0),
			300
		);*/
    	    	
    };

//**************************************************************************//

};  //End namespace Foam

#endif
