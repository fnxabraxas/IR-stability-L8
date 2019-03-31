

	integer, parameter :: Nbmax=5, NepsAmax=5
	double precision, parameter :: c=1.d0, cero=1.d-15,toler=1.d-15,
     +							dpcheat=0.01d0
	integer  iaction(2,2),imoral(2,2,2), Nb, NepsA
	double precision epsA,b,pcheat,pdis,epsAini,
     +			 bvect(Nbmax), epsAvect(NepsAmax)

	common /ppal/  epsA, bvect,b, epsAvect,
     +			pcheat,pdis,epsAini,
     +			Nb, NepsA
	common /xx/ iaction,imoral
