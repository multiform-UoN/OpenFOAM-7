//- Create controls for the power iterations

// Aitken controls
scalar aitkenInitialRes(0);
scalar aitkenInitialResPsi(0);
scalar aitkenInitialResPsiAdj(0);
scalar aitkenTheta(0);
scalar aitkenThetaPsi(0);
scalar aitkenThetaPsiAdj(0);

scalar aitkenDeltaLambda(0);
scalar aitkenDeltaLambdaAdj(0);

//- Residual on the eigenvalue
scalar lambdaRes(1.);
scalar psiRes(1.);
scalar psiAdjRes(1.);
scalar XRes(1.);
label currPowerIter(0);
label currCellIter(0);
bool powerConverged(false);
bool cellConverged(false);
simpleControl pwrctrl(mesh,"PowerControl");

const dictionary& pwrDict(pwrctrl.dict());

//- Tolerance for the eigenvalue
scalar eigenValTol
(
    readScalar(pwrDict.lookup("eigenvalueTol"))
);

//- Tolerance for the eigenfunction
scalar eigenFunTol
(
    readScalar(pwrDict.lookup("eigenfunctionTol"))
);

//- Tolerance for the first order corrector
scalar corrTol
(
    readScalar(pwrDict.lookup("cellProblemTol"))
);

//- Max iterations
label maxPowerIter
(
    readLabel(pwrDict.lookup("maxPowerIterations"))
);

//- Relaxation for the eigenvalue
scalar theta
(
    readScalar(pwrDict.lookup("powerIterationRelaxation"))
);

//- Relaxation for the eigenvalue
scalar thetaShift
(
    readScalar(pwrDict.lookup("powerIterationRelaxationShift"))
);

//- Shifting
scalar shift
(
    readScalar(pwrDict.lookup("powerIterationShiftRatio"))
);

//- Relaxation for the eigenvalue
scalar lambda0
(
    readScalar(pwrDict.lookup("initialEigenvalue"))
);

//- Relaxation for the eigenvalue
scalar lambdaAdj0
(
    readScalar(pwrDict.lookup("initialAdjEigenvalue"))
);


label maxCellProbIter
(
    readLabel(pwrDict.lookup("maxCellProblemIterations"))
);

Switch AitkenRelax(pwrDict.lookupOrDefault<Switch>("AitkenRelaxation",true));
if(AitkenRelax)
{
	Info<<"\nUsing Aitken relaxation scheme for the corrector problem" << endl;
}

label aitkenUpdate(pwrDict.lookupOrDefault<label>("AitkenUpdateFrequency",100));
