//- Lookup phaseFieldControls in PIMPLE
dictionary pfControls(pimple.dict().subDict("phaseFieldControls"));

label pfMaxCorr(readLabel(pfControls.lookup("maxNofCorrectors")));
scalar alphaTol(readScalar(pfControls.lookup("alphaTol")));
scalar muTol(readScalar(pfControls.lookup("muTol")));
