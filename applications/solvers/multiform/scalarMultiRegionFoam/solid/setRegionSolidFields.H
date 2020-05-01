fvMesh& mesh = solidRegions[i];

volScalarField& T = fieldssolid[i];

fv::options& fvOptions = solidFvOptions[i];

solidNoLoopControl& pimple = pimples.solid(i);
