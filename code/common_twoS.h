
	double precision, parameter :: xtol=1.d-10, c=1.d0, At=1.d-2  ! try different At
	integer, parameter :: nmaxb=5   
	double precision b, epsA, epsM, epsMUT, x(2,2,2) , bvect(nmaxb)
	integer yaction(2,2,2),ymoral(2,2,2,2), nib
	
	common /ppal/  b, epsA, epsM, epsMUT , bvect, nib
	common /invasion1/ yaction,ymoral,x


