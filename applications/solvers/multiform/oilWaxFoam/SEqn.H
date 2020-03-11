
const volScalarField ceq
  ("ceq",
    //pow (c * (T/Trif) ,-2)
    crif + (cT * (T-Trif))
  );

const volScalarField H  // Heaviside
  ("H",
    Foam::pos(
      c
      -
      ceq
    )
  );

const volScalarField ScImp   // Source term
  ("ScImp",
    H/tau1
    +
    (scalar(1)-H)*(psi/tau2)
  );

const volScalarField ScExp   // Source term
  ("ScExp",
    -
    H*ceq/tau1
    -
    (scalar(1)-H)*(psi/tau2)*ceq
  );

const volScalarField SpsiImp   // Source term
  ("SpsiImp",
    -(scalar(1)-H)*(c-ceq)/tau2*psifact*w
  );

const volScalarField SpsiExp   // Source term
  ("SpsiExp",
    -H*(c-ceq)/tau1*psifact*w
  );

// const volScalarField Sctot
//   ("Sctot",
//     ScImp*c+ScExp
//   );
//
// const volScalarField Spsitot
//   ("Spsitot",
//     SpsiImp*psi+SpsiExp
//   );
//
// Sctot.write();
// Spsitot.write();
// ceq.write();
// H.write();
