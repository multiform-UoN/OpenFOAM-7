//- Disable screen output from linear solver
//  (Otherwise it looks like the solver never converges)
lduMatrix::debug = 0;

//- Store previous iteration fields
psi.storePrevIter();
psiAdj.storePrevIter();

//- Store previous eigenvalues
scalar lambdaOld = lambda.value();
scalar lambdaAdjOld = lambdaAdj.value();
