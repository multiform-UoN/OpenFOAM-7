Uslip = - St
         *
         //(I-St*Gr.T())
         //&
        (
          fvc::ddt(U) +
          (U & fvc::grad(U))
        );

phislip = fvc::flux(Uslip);
