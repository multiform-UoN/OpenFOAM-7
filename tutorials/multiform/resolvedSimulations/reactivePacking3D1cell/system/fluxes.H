functions
{
    breakthrough
    {
        // Load the library containing the 'coded' functionObject
        libs            ("libutilityFunctionObjects.so");

        type coded;

        // Name of on-the-fly generated functionObject
        name breakthrough;

        codeOptions
        #{
            -I../../lnInclude
        #};

        codeInclude
        #{
            #include "dynamicRetentionFvPatchScalarField.H"
        #};

        codeExecute
        #{
          // Lookup c
          //Info<< "Looking up field c\n" << endl;
          const volScalarField& c = mesh().lookupObject<volScalarField>("T");
          const surfaceScalarField& phi = mesh().lookupObject<surfaceScalarField>("phi");
          const volVectorField& U = mesh().lookupObject<volVectorField>("U");
          const surfaceScalarField phic(fvc::flux(U)*fvc::interpolate(c));

          const label nfacezone = mesh().faceZones().size();
          IOdictionary controlDict
            (
                IOobject
                (
                "controlDict",
                mesh().time().system(),
                mesh(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
                )
            );

          const wordList patchList(controlDict.lookup("breakthrough"));
          const wordList patchListS(controlDict.lookup("retention"));
          const label face_tot = nfacezone + patchList.size();
          scalarField sarea(face_tot,0);
          scalarField sflux(face_tot,0);
          scalarField flux(face_tot,0);
          const surfaceScalarField& magSf = mesh().magSf();
          const surfaceVectorField& Sf = mesh().Sf();

          // loop through all faceZones
          forAll(mesh().faceZones(),iz)
          {
              forAll(mesh().faceZones()[mesh().faceZones().names()[iz]],iface)
              {
                  const label facei = mesh().faceZones()[mesh().faceZones().names()[iz]][iface];
                  //const scalar voli = magSf[facei];//&mesh().faces()[facei].normal(mesh().points());
                  //sarea[iz] += voli;
                  sflux[iz] += phic[facei];
                  flux[iz] += phi[facei];
              }
          }


          //Loop through all faces
          forAll(patchList,ipatch)
          {
            label patchI = mesh().boundaryMesh().findPatchID(patchList[ipatch]);
            //sarea[nfacezone+ipatch] = gSum(magSf.boundaryField()[patchI]);
            sflux[nfacezone+ipatch] = gSum(phic.boundaryField()[patchI]);
            flux[nfacezone+ipatch] = gSum(phi.boundaryField()[patchI]);
          }

          //reduce(sarea,sumOp<tmp<Field<double>>>());
          reduce(sflux,sumOp<tmp<Field<double>>>());
          reduce(flux,sumOp<tmp<Field<double>>>());

          if (Pstream::master())
          {
            fileName outputDir = "./breakthrough";
            OFstream is(outputDir/"dirExists.dat");
            autoPtr<std::ofstream> file;
            if (!is.good())
            {
              Info << "Creating new breakthrough folder" << endl;
              mkDir(outputDir);
              file.reset(new std::ofstream(outputDir/"breakthrough.dat"));
              file() << "Time\t";
              forAll(mesh().faceZones().names(),iname)
              {
                file() << mesh().faceZones().names()[iname] << "\t" << mesh().faceZones().names()[iname] << "\t";
              }
              forAll(patchList,ipatch)
              {
                file() << patchList[ipatch] << "\t" << patchList[ipatch] << "\t";
              }
              forAll(patchListS,ipatch)
              {
                file() << patchListS[ipatch] << "\t" << patchListS[ipatch] << "\t";
              }
              file() << " volume \t average";
              file() << "\n";
            }
            file.reset(new std::ofstream(outputDir/"breakthrough.dat", std::ios::out | std::ios::app));
            file() << mesh().time().timeName();
            forAll(sflux,i)
            {
                file()  << "\t" << sflux[i] << "\t" << sflux[i]/flux[i];
            }
            forAll(patchListS,ipatch)
            {
              label patchI = mesh().boundaryMesh().findPatchID(patchListS[ipatch]);
//              if (isA<dynamicRetentionFvPatchScalarField>(c.boundaryField()[patchI]))
                {
                    const dynamicRetentionFvPatchScalarField& fluxPatch
                    (
                        refCast<const dynamicRetentionFvPatchScalarField>
                        (
                            c.boundaryField()[patchI]
                        )
                    );

                    file() << "\t" << gSum(fluxPatch.patch().magSf());
                    file() << "\t" << fluxPatch.Savg() ;
                }
            }
            file() << "\t" << gSum(mesh().V());
            file() << "\t" << gSum(c.primitiveField()*mesh().V())/gSum(mesh().V());
            file() << "\n";

          }
          #};
      }
  }
