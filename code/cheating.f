
      include "sub_mutation_two.f"

      program stability_cheat
c     Calculate the P_dis curve
      implicit none
      include "common_cheat.h"
      integer, parameter :: Nmax=99999, Nmaxcol=99
      character*80  outf
      integer i,j,k,iresident,ii,jj,kk, ncheat,ib,ie
      double precision  calcpayCH, difpayold,paso,pdisold,
     +			pdisP,pdisN, difpayN, difpayP,
     +      pdisv(Nmax,Nmaxcol), pcheatv(Nmax),
     +		difpayA,binp

      print*,'Output file: '
      read(*,*) outf
      print*,'Strategy: '
      read(*,*) iresident
      call num2actionmoral(iresident)
      print*,'n_cheat: '
      read(*,*) ncheat

      Nb=4
      bvect=(/1.2d0,1.5d0,2.d0,3.d0,0.d0 /)
      NepsA=3
      epsAvect=(/1.d-1,1.d-2,1.d-3,0.d0,0.d0 /)

      pdisv=0.d0
      pcheatv=0.d0

      open(90,file=outf,status='UNKNOWN')
      write(90,'(A1,2A,4F6.2,A,3F8.4)') '#','  Pcheat   Pdiscov:  ',
     +	'  b= ',(bvect(j),j=1,Nb),'       epsA= ',(epsAvect(j),j=1,NepsA)

      paso=0.1d0
      jj=0
      do ib=1,Nb
	b=bvect(ib)
      do ie=1,NepsA
	epsAini=epsAvect(ie)
	epsA=epsAini

        jj=jj+1

      do i=1,ncheat

        pcheat=(1.d0/ncheat)*i

        pdis=-paso
        pdisold=1.d0
        difpayold=9999.d0
        difpayA=-9999.d0
        do while (difpayA.lt.0.d0)
	    pdis=pdis+paso
           if (pdis.gt.1.d0) then
              pdisv(i+1,jj)=99.d0
              goto 51
           endif
           difpayA=calcpayCH()
	 enddo

	 difpayN=difpayold
	 difpayP=difpayA
	 pdisN=pdis-paso
	 pdisP=pdis

        do while (dabs(2.d0*(pdisN-pdisP)/(pdisN+pdisP)).gt.toler)
	    pdis=(pdisN+pdisP)/2.d0
            difpayA=calcpayCH()
	    if(difpayA.lt.0.d0) then
		difpayN=difpayA
		pdisN=pdis
	    else
		difpayP=difpayA
		pdisP=pdis
	    endif
	enddo
	pdisv(i+1,jj)=pdis
 51     continue
        pcheatv(i+1)=pcheat
        
      enddo
      enddo
      enddo

      do i=0,ncheat
         write(90,'(F8.4,$)') pcheatv(i+1)
	  do jj=1,Nb
		write(90,'(A,$)') '    '
		do kk=1,NepsA
		  write(90,'(F8.4,$)') pdisv(i+1,(jj-1)*NepsA+kk)
		enddo
	  enddo
	  write(90,*) ' '
      enddo

      stop
      end
	
	
	
C---------------------------------------------------------------------------

	double precision function DamH(alp,bet)
c	Calculate P(\alpha,\beta) in an homogeneous population
	implicit none
      	include "common_cheat.h"
	integer alp,bet
		DamH = (1.d0-epsA)*imoral(2-alp,2-bet,2-iaction(2-alp,2-bet))
     +		 + epsA*imoral(2-alp,2-bet,2)
	return
	end
		

		
	double precision function calcpayCH()
c	Calculate difference of payoffs under cheating
	implicit none
	double precision payoff1, payoff2,payoff3,payoff4
      	include "common_cheat.h"
	double precision p11,p12,p21,p22,x1G,x1B,x2G,x2B,pp,
     +			num,den,pcd,D11,D10,D01,D00

	epsA=epsAini
	epsA=epsA+(1.d0-epsA)*pcheat*pdis
	x1G=pp(1.d0,0.d0)
	epsA=epsAini
	x1B=1.d0-x1G
	pcd=(pcheat+dpcheat)*pdis
	D11=(1.d0-epsA)*(1.d0-pcd)*imoral(1,1,2-iaction(1,1))
     +		+(epsA+pcd-epsA*pcd)*imoral(1,1,2)
	D10=(1.d0-epsA)*(1.d0-pcd)*imoral(1,2,2-iaction(1,2))
     +		+(epsA+pcd-epsA*pcd)*imoral(1,2,2)
	D01=(1.d0-epsA)*(1.d0-pcd)*imoral(2,1,2-iaction(2,1))
     +		+(epsA+pcd-epsA*pcd)*imoral(2,1,2)
	D00=(1.d0-epsA)*(1.d0-pcd)*imoral(2,2,2-iaction(2,2))
     +		+(epsA+pcd-epsA*pcd)*imoral(2,2,2)
	num= x1G*D01+x1B*D00
      	den= 1.d0 + x1G*D01 +x1B*D00 -x1G*D11 -x1B*D10
	x2G=num/den
	x2B=1.d0-x2G

	p11=iaction(1,1)*x1G*x1G+iaction(2,2)*x1B*x1B+
     +    iaction(1,2)*x1G*x1B+iaction(2,1)*x1G*x1B  ! 1 to 1
	p12=iaction(1,1)*x1G*x2G+iaction(2,2)*x1B*x2B+
     +    iaction(1,2)*x1G*x2B+iaction(2,1)*x1B*x2G  ! 1 to 2
	p21=iaction(1,1)*x2G*x1G+iaction(2,2)*x2B*x1B+
     +    iaction(1,2)*x2G*x1B+iaction(2,1)*x2B*x1G  ! 2 to 1
	p22=iaction(1,1)*x2G*x2G+iaction(2,2)*x2B*x2B+
     +    iaction(1,2)*x2G*x2B+iaction(2,1)*x2B*x2G  ! 2 to 2

	p21=p21*(1.d0-(pcheat+dpcheat))
	p22=p22*(1.d0-(pcheat+dpcheat))
	p11=p11*(1.d0-pcheat)
	p12=p12*(1.d0-pcheat)

	payoff1=((b-c)*p11  )  *(1.d0-epsA)    ! 1 (with 1)
	payoff2=(b*p12-c*p21)  *(1.d0-epsA)    ! 2 (with 1)
	payoff3=(b*p21-c*p12)  *(1.d0-epsA)    ! 1 (with 2)
	payoff4=((b-c)*p22  )  *(1.d0-epsA)    ! 2 (with 2)

        calcpayCH=payoff1-payoff2

	return
	end	
		
c ====================================================================================

	double precision function pp(ap,bp)
c	Calculate reputation in an homogeneous population (ap=1,bp=0)
	implicit none
      	include "common_cheat.h"
	integer PosInd
	double precision ap,bp, AA,BB,CC,DamH, solveSEC
	AA=ap*(DamH(1,1)-DamH(1,0)-DamH(0,1)+DamH(0,0))
	BB = -1.d0+ ap*(DamH(1,0)+DamH(0,1)-2.d0*DamH(0,0))
	CC = ap*DamH(0,0) +bp
	if((abs(AA).lt.cero).and.(abs(BB).lt.cero))then
		stop 'Error en pp: AA=BB=0'
	else
		pp=solveSEC(AA,BB,CC)
	endif
	return
	end


c------------ Computation of second degree equation -----------------------------------------------------

	double precision function solveSEC(AA,BB,CC)
	implicit none
	include "common_cheat.h"
	double precision AA,BB,CC, xx,r,num,den,x1,x2

	if (abs(AA).lt.cero) then
		if (abs(BB).lt.cero) then
			xx=0.5d0   ! degenerated case
		else
			xx=-CC/BB
		endif
	elseif (abs(BB).lt.cero) then
		xx=-CC/AA
		if(xx.lt.0.d0) xx=-xx
	elseif (abs(CC).lt.cero) then
		xx=-BB/AA
		if(AA.gt.0.d0) then
			xx=0.d0
		elseif(xx.lt.0.d0) then
			xx=0.d0
		endif
	else
		r=(BB**2.d0)-(4.d0*AA*CC)
		if (r.lt.0.d0) then
			print*,'Error!!'
		else
			num=(-BB-(r**0.5d0))
			num=(anint(num*1.d0/cero))*cero
			x1=num/(2.d0*AA)   ! the lowest solution
			num=(-BB+(r**0.5d0))
			num=(anint(num*1.d0/cero))*cero
			x2=num/(2.d0*AA)
			if((x1.lt.0.d0).or.(x1.gt.1.d0)) then
				if((x2.lt.0.d0).or.(x2.gt.1.d0)) then
					print*,'Error: x_H out of range'
				else
					xx=x2
				endif
			else
				xx=x1 ! the lowest solution is the stable one
			endif
		endif
	endif

	return
	end





