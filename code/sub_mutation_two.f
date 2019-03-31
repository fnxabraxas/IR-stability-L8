	subroutine 
     +		doESS(iresident,imutant,xir,xim,payoff1,difpayA,difpayB)
c	Compute the difference of payoffs between the resident and mutant strategies	
	implicit none
	include "common_twoS.h"
	integer iresident,imutant
	double precision xir,xim,payoff1,difpayA(nmaxb),difpayB(nmaxb)
      	integer i,j,k,jj, ib, errorf,a1,a2,b1,b2,is
      	double precision payoff2,payoff3,payoff4,

	print*, 'Resident and mutant: ',iresident,imutant
  	call iniall_two(iresident,imutant,xir,xim,1)
	call solve1
	call solve2(errorf)

	if(errorf.ne.0) then  ! solving equations all together (necessary in degenerate cases)
			print*,'Computing again...'
			call iniall_two(iresident,imutant,xir,xim,1)
			call solve12
	endif

 	do ib=1,nib
	  b=bvect(ib)
	  call calcpay(payoff1,payoff2,payoff3,payoff4)
	  difpayA(ib)=payoff1-payoff2
	  difpayB(ib)=payoff3-payoff4
	enddo
	payoff1=payoff1/(b-c)/(1.d0-epsA)

	return
	end


	subroutine iniall_two(iresident,imutant,xres,xmut,mi)
	
	implicit none
	include "common_twoS.h"
	integer iresident,imutant, mi
	double precision xres,xmut
	integer i,j,k,l,l1,l2,iactionv(2,2),imoralv(2,2,2),ylabel(2)
	double precision Cfun
	double precision Cfun,epsRUMOR1,epsRUMOR2

	epsRUMOR1=0.d0 ! change to check effect of misjudgements
	epsRUMOR2=0.d0 ! change to check effect of misjudgements

	ymoral=-1
	yaction=-1
	ylabel(1)=iresident
	ylabel(2)=imutant
	if(mi.eq.1) then
		do j=1,2
			x(j,1,1)=xres*(1.d0-epsRUMOR1)
			x(j,1,2)=xres*epsRUMOR1 
			x(j,2,2)=(1.d0-xres)*(1.d0-epsRUMOR2)
			x(j,2,1)=(1.d0-xres)*epsRUMOR2
		enddo
	else
  		stop 'Error 01 - iniall_two'
	endif
	do i=1,2
 	 	call num2actionmoral(ylabel(i),iactionv,imoralv)
  	do l=1,2
   	do j=1,2
   	do k=1,2
		ymoral(i,l,j,k)=imoralv(l,j,k)
   	enddo
		yaction(i,l,j)=iactionv(l,j)
  	enddo
  	enddo
	enddo

	return
	end



c----------- Solvers ------------------------------------------------------------------------------


	subroutine solve12
c	Solves all the equations together
      	implicit none
      	include "common_twoS.h"
		integer i,j,k,l,l1,l2,a1,a2,b1,b2,is,cont,converge(2)
		double precision x1G, sumaGG(2),sumaBB(2), sumaGB(2),sumaBG(2),
     +		 xGGold(2), xGBold(2),xBGold(2), xBBold(2), PC, difx2(2),
     +       sumaT

	x1G=x(1,1,1)+x(1,1,2)

	difx2=1.d0
	cont=0
	converge=0
	do while ((converge(1).ne.1).or.(converge(2).ne.1))
	  cont=cont+1
	  do is=1,2  ! resident and mutant
	  	sumaGG(is)=0.d0
	 	sumaBB(is)=0.d0	
	  	sumaGB(is)=0.d0
	  	sumaBG(is)=0.d0	
	  enddo !is

	  do a1=0,1
	  do a2=0,1
	  do b1=0,1
	  do b2=0,1

		sumaGG(1)=sumaGG(1)+ x(1,2-a1,2-a2)*x(1,2-b1,2-b2)*
     +  				PC(1,1,1,1,a1,b1,a1,b1,2,1,a2,b2,a1,b1) 
		sumaBB(1)=sumaBB(1)+ x(1,2-a1,2-a2)*x(1,2-b1,2-b2)*
     +  				PC(0,0,1,1,a1,b1,a1,b1,2,1,a2,b2,a1,b1) 
		sumaGB(1)=sumaGB(1)+ x(1,2-a1,2-a2)*x(1,2-b1,2-b2)*
     +  				PC(1,0,1,1,a1,b1,a1,b1,2,1,a2,b2,a1,b1) 
		sumaBG(1)=sumaBG(1)+ x(1,2-a1,2-a2)*x(1,2-b1,2-b2)*
     +  				PC(0,1,1,1,a1,b1,a1,b1,2,1,a2,b2,a1,b1)
	  enddo
	  enddo
	  enddo
	  enddo

	  do a1=0,1
	  do a2=0,1
	  do b1=0,1
	  do b2=0,1	
		sumaGG(2)=sumaGG(2)+ x(2,2-a1,2-a2)*x(1,2-b1,2-b2)*
     +  				PC(1,1,1,2,a1,b1,a2,b2,2,2,a2,b2,a2,b2) 
		sumaBB(2)=sumaBB(2)+ x(2,2-a1,2-a2)*x(1,2-b1,2-b2)*
     +  				PC(0,0,1,2,a1,b1,a2,b2,2,2,a2,b2,a2,b2) 
		sumaGB(2)=sumaGB(2)+ x(2,2-a1,2-a2)*x(1,2-b1,2-b2)*
     +  				PC(1,0,1,2,a1,b1,a2,b2,2,2,a2,b2,a2,b2) 
		sumaBG(2)=sumaBG(2)+ x(2,2-a1,2-a2)*x(1,2-b1,2-b2)*
     +  				PC(0,1,1,2,a1,b1,a2,b2,2,2,a2,b2,a2,b2)
	  enddo
	  enddo
	  enddo
	  enddo

	  difx2=0.d0
	  do is=1,2
		sumaT=sumaGG(is)+sumaBB(is)+sumaGB(is)+sumaBG(is)
		sumaGG(is)=sumaGG(is)/sumaT
		sumaBB(is)=sumaBB(is)/sumaT
		sumaGB(is)=sumaGB(is)/sumaT
		sumaBG(is)=sumaBG(is)/sumaT
	  	xGGold(is)=x(is,1,1)
	  	xGBold(is)=x(is,1,2)
	  	xBGold(is)=x(is,2,1)
	  	xBBold(is)=x(is,2,2)
	  	x(is,1,1)= At*sumaGG(is) + (1.d0-At)*x(is,1,1)
	  	x(is,2,2)= At*sumaBB(is) + (1.d0-At)*x(is,2,2)
	  	x(is,1,2)= At*sumaGB(is) + (1.d0-At)*x(is,1,2)
	  	x(is,2,1)= At*sumaBG(is) + (1.d0-At)*x(is,2,1)
	  	difx2(is)=  (x(is,1,1)-XGGold(is))**2.d0
     +		+(x(is,1,2)-XGBold(is))**2.d0
     + 		+(x(is,2,1)-XBGold(is))**2.d0+(x(is,2,2)-XBBold(is))**2.d0
		if(difx2(is)**0.5d0.lt.xtol) converge(is)=1
	  enddo

	enddo !while

	return 
	end
 

	subroutine solve1
c	Solves resident strategy equations
     	implicit none
     	include "common_twoS.h"
	integer i,j,k,l,l1,l2,a1,a2,b1,b2
	double precision x1G, sumaGG,sumaBB, difx2, sumaGB,sumaBG,sumaT,
     +				  xGGold, xGBold,xBGold, xBBold, PC

	x1G=x(1,1,1)+x(1,1,2)
	difx2=1.d0
	do while (difx2**0.5d0.gt.xtol)
	  sumaGG=0.d0
	  sumaBB=0.d0	
	  sumaGB=0.d0
	  sumaBG=0.d0	
	  do a1=0,1
	  do a2=0,1
	  do b1=0,1
	  do b2=0,1
		sumaGG=sumaGG+ x(1,2-a1,2-a2)*x(1,2-b1,2-b2)*
     +  				PC(1,1,1,1,a1,b1,a1,b1,2,1,a2,b2,a1,b1) 
		sumaBB=sumaBB+ x(1,2-a1,2-a2)*x(1,2-b1,2-b2)*
     +  				PC(0,0,1,1,a1,b1,a1,b1,2,1,a2,b2,a1,b1) 
	  enddo
	  enddo
	  enddo
	  enddo
	  xGGold=x(1,1,1)
	  xGBold=x(1,1,2)
	  xBGold=x(1,2,1)
	  xBBold=x(1,2,2)
	  x(1,1,1)= At*sumaGG + (1.d0-At)*x(1,1,1)
	  x(1,2,2)= At*sumaBB + (1.d0-At)*x(1,2,2)
	  x(1,1,2)=x1G-x(1,1,1)
	  x(1,2,1)=1.d0-x1G-x(1,2,2)

	  difx2=(x(1,1,1)-XGGold)**2.d0+(x(1,1,2)-XGBold)**2.d0
     +		+(x(1,2,1)-XBGold)**2.d0+(x(1,2,2)-XBBold)**2.d0
	enddo

	return 
	end
 

	subroutine solve2(errorf)
c	Solves resident mutant equations
      	implicit none
      	include "common_twoS.h"
	integer i,j,k,l,ai,b1,b2,errorf,chk,  l1,l2,a1,a2
	double precision x1G,x1B,x2G, Dam,D11,D10,D01,D00,
     +				A(2,2),bm(2),Ainv(2,2),xsol(2), PC,num,den

	x1G=x(1,1,1)+x(1,2,1)  ! Cuidado: no es el mismo x1G que en solve1
	x1B=1.d0-x1G
	D11=(1.d0-epsA)*ymoral(2,1,1,2-yaction(2,1,1))
     +		+epsA*ymoral(2,1,1,2)
	D10=(1.d0-epsA)*ymoral(2,1,2,2-yaction(2,1,2))
     +		+epsA*ymoral(2,1,2,2)
	D01=(1.d0-epsA)*ymoral(2,2,1,2-yaction(2,2,1))
     +		+epsA*ymoral(2,2,1,2)
	D00=(1.d0-epsA)*ymoral(2,2,2,2-yaction(2,2,2))
     +		+epsA*ymoral(2,2,2,2)
	num= x1G*D01+x1B*D00 
      	den= 1.d0 + x1G*D01 +x1B*D00 -x1G*D11 -x1B*D10 
	if (den.lt.1.d-15) then
		print*,'---'
		if(num.lt.1.d-15) then
			errorf=1
			goto 51
		else
			stop 'Error in solve2: den=0, num!=0'
		endif
	else
		x2G=num/den
	endif

	A=0.d0
	bm=0.d0
	do b1=0,1
	do b2=0,1
	 do l=1,2
	 do ai=1,2
		A(l,ai)=A(l,ai)+  x(1,2-b1,2-b2) *( 
     +  		PC(2-l,2-l,2,2,2-ai,b2,2-ai,b2,1,2,2-ai,b1,2-ai,b2)
     +  	    - PC(2-l,2-l,2,2,2-ai,b2,2-ai,b2,1,2,ai-1,b1,2-ai,b2) )
	 enddo
	 bm(l)= bm(l)+  x(1,2-b1,2-b2)*( 
     +	         x2G*PC(2-l,2-l,1,2,0,b1,1,b2,2,2,1,b2,1,b2)
     +	  +(1.d0-x2G)*PC(2-l,2-l,1,2,1,b1,0,b2,2,2,0,b2,0,b2)   )
	 enddo
	enddo
	enddo
	A(1,1)=A(1,1)-1.d0
	A(2,2)=A(2,2)-1.d0
	bm=-bm
	call FINDInv(A,Ainv,2,2,errorf)
	if (errorf.ne.0) then
		goto 51
	endif
	xsol=matmul(Ainv,bm)
	x(2,1,1)=xsol(1)
	x(2,2,2)=xsol(2)
	x(2,2,1)=x2G-x(2,1,1)
	x(2,1,2)=1.d0-x2G-x(2,2,2)

	do i=1,2
	do j=1,2
	do k=1,2
		x(i,j,k)=anint(x(i,j,k)*1.d15)*1.d-15
		if(x(i,j,k).lt.-1.d-14) then
			print*,'Solution out of the simplex'
			errorf=1
			goto 51
		endif
	enddo
	enddo
	enddo

 51	continue
	return 
	end



c------------ payoff -----------------------------------------------------------------

	subroutine calcpay(payoff1,payoff2,payoff3,payoff4)
c	Computes payoffs
	implicit none
	double precision payoff1, payoff2,payoff3,payoff4
	include "common_twoS.h"
	double precision p11,p12,p21,p22, 
     +		x1G1,x1B1,x1G2,x1B2,x2G1,x2B1,x2G2,x2B2
	
	x1G1=x(1,1,1)+x(1,1,2) ! 1 is considered G by 1
	x1B1=1.d0-x1G1
	x1G2=x(1,1,1)+x(1,2,1) ! 1 is considered G by 2
	x1B2=1.d0-x1G2
	x2G1=x(2,1,1)+x(2,1,2) ! 2 is considered G by 1
	x2B1=1.d0-x2G1
	x2G2=x(2,1,1)+x(2,2,1) ! 2 is considered G by 2
	x2B2=1.d0-x2G2

	p11=yaction(1,1,1)*x1G1*x1G1+yaction(1,2,2)*x1B1*x1B1+
     +    yaction(1,1,2)*x1G1*x1B1+yaction(1,2,1)*x1G1*x1B1  ! 1 to 1
	p12=yaction(1,1,1)*x1G1*x2G1+yaction(1,2,2)*x1B1*x2B1+
     +    yaction(1,1,2)*x1G1*x2B1+yaction(1,2,1)*x1B1*x2G1  ! 1 to 2
	p21=yaction(2,1,1)*x2G2*x1G2+yaction(2,2,2)*x2B2*x1B2+
     +    yaction(2,1,2)*x2G2*x1B2+yaction(2,2,1)*x2B2*x1G2  ! 2 to 1
	p22=yaction(2,1,1)*x2G2*x2G2+yaction(2,2,2)*x2B2*x2B2+
     +    yaction(2,1,2)*x2G2*x2B2+yaction(2,2,1)*x2B2*x2G2  ! 2 to 2

	payoff1=((b-c)*p11  )  *(1.d0-epsA)    ! of 1 (with 1)
	payoff2=(b*p12-c*p21)  *(1.d0-epsA)    ! of 2 (with 1)
	payoff3=(b*p21-c*p12)  *(1.d0-epsA)    ! of 1 (with 2)
	payoff4=((b-c)*p22  )  *(1.d0-epsA)    ! of 2 (with 2)

	return
	end


c=======================================================================================================

      double precision function Cfun(ip,xf)
c     computes \chi function
      implicit none
      integer ip
      double precision xf
      if (ip.eq.0) then
	Cfun=1.d0-xf
      elseif (ip.eq.1) then
	Cfun=xf
      else
	stop "Error 01 in Cfun"
      endif
      return
      end

      integer function iCfun(ip,xf)
c     computes \chi function (integer)
      implicit none
      integer ip, xf
      if (ip.eq.0) then
	iCfun=1-xf
      elseif (ip.eq.1) then
	iCfun=xf
      else
	stop "Error 01 in iCfun"
      endif
      return
      end


	subroutine num2actionmoral(inum,iactionv,imoralv)
c	Write action and moral code in arrays from decimal-binary representation
c	GG(0) GB(1) BG(2) BB(3) GGC(4) GGD(5) GBC(6) GBD(7)...	
	implicit none
	integer inum,iactionv(2,2),imoralv(2,2,2)
	integer inumt,i,j,k
	iactionv=-1
	imoralv=-1
	inumt=inum
	do i=1,2
	 do j=1,2
		iactionv(i,j)=mod(inumt,2)
		inumt=floor(inumt/2.)
	 enddo
	enddo
	do i=1,2
	 do j=1,2
	  do k=1,2
		imoralv(i,j,k)=mod(inumt,2)
		inumt=floor(inumt/2.)
	  enddo
	 enddo
	enddo
	return
	end


	double precision function 
     +	PC(l1,l2,i,j,alpi,beti,alpj,betj,i2,j2,alpi2,beti2,alpj2,betj2)
c	Calculate P(\alpha,\beta,...) 
	implicit none
	integer l1,l2,i,j,alpi,beti,alpj,betj,
     +				i2,j2,alpi2,beti2,alpj2,betj2, iCfun
	double precision Cfun
	include "common_twoS.h"
	PC=(1.d0-epsA)
     +*iCfun(l1,ymoral(i,2-alpi,2-beti,2-yaction(j,2-alpj,2-betj)))
     +*iCfun(l2,ymoral(i2,2-alpi2,2-beti2,
     +                                 2-yaction(j2,2-alpj2,2-betj2)))
     +			+epsA *iCfun(l1,ymoral(i,2-alpi,2-beti,2))
     +					*iCfun(l2,ymoral(i2,2-alpi2,2-beti2,2))
      	return
      	end


	integer function equivS(inum)
c	Find mirror strategy
	implicit none
	integer inum
	integer cont,i,j,k,iactionv(2,2),imoralv(2,2,2)
	
	equivS=0
	call num2actionmoral(inum,iactionv,imoralv)
	cont=-1
	do i=1,2
	 do j=1,2
	  cont=cont+1
	  equivS=equivS+iactionv(3-i,3-j)*(2**cont)
	 enddo
	enddo
	do i=1,2
	 do j=1,2
	  do k=1,2
	   cont=cont+1
           equivS=equivS+(1-imoralv(3-i,3-j,k))*(2**cont)
	  enddo
	 enddo
	enddo

	return
	end


