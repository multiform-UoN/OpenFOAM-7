/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  7
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Calculates the normal gradient at chosen patches.

\*---------------------------------------------------------------------------*/
functionObjectLibs            ("libutilityFunctionObjects.so");

type            coded;

name E_Slice;

codeWrite
#{
	/*
		TO DO, STEPS:
			- Construct point list along line (DONE)
			- Find cell ids from point list (DONE)
			- Perform -grad(V) = E along cells
				- !!!!!!!!!!!!!!!!!!! LOOK AT gaussGrad.C for "INSPIRATION" !!!!!!!!!!!!!!!!!!!!!
			- Output to file
	*/
	//Create vectors of start and end points
	vector Start_Point(2.5e-6, 0, 0);
	vector End_Point(2.5e-6, 5e-5, 0);

	//Construct mesh searching object for mesh
	meshSearch ms(mesh());

	//Distance between points in y-dir
	double epsilon = 5e-8;

	//Number of points
	int Num_Points = 1000;

	//pointField with Num_Points, starting at (2.5e-6,0,0)
	pointField Line_Points(Num_Points,point(2.5e-6, 0, 0));
	//Label list for cell indexes
	labelList Cell_IDs(Num_Points);
	//Find cell id of starting cell
	Cell_IDs[0] = ms.findCell(Line_Points[0]);

	//Create pointField along desired line, with spacing epsilon in y-dir
	for(int i=1;i<Num_Points;i++){
		Line_Points[i] = vector(2.5e-6, epsilon*i, 0);
		//Find cell IDs containing points
		Cell_IDs[i] = ms.findCell(Line_Points[i]);
	}
	//Run through each cell and calculate grad(V)
	forAll(Cell_IDs,cell){

		}

#};


// ************************************************************************* //
