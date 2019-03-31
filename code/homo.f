
	program homo
C	Compute reputations in homogeneous populations

	implicit none
	include "common_homo.h"
	character*80 outp
	double precision A,B,C, x, r, x1,x2, h, DamH,Cfun,num
	integer i,j,k,l,iA

	print*,'Errors in ACTION: '
	read(*,*) epsA
	print*,'Output file: '
	read(*,*) outp

	x=-1
	open(90,file=outp,status='unknown')
      	write(90,'(A1,A)') '#',' Strg    x_H   h_H'
      	do iA=0,4095 				! strategy
			call num2actionmoral(iA)
			A = DamH(1,1,1) - DamH(1,1,0) 
     &			- DamH(1,0,1) + DamH(1,0,0)
			B = -1.d0 + DamH(1,1,0) 
     &					+ DamH(1,0,1) - 2.d0*DamH(1,0,0)
			C = DamH(1,0,0)

			if (abs(A).lt.1.d-15) then
				if (abs(B).lt.1.d-15) then
					x=0.5d0  ! degenerated case
				else
					x=-C/B
				endif
			elseif (abs(B).lt.1.d-15) then
				x=-C/A
				if(x.lt.0.d0) x=-x
			elseif (abs(C).lt.1.d-15) then
				x=-B/A 
				if(A.gt.0.d0) then
					x=0.d0
				elseif(x.lt.0.d0) then
					x=0.d0
				endif
			else
				r=(B**2.d0)-(4.d0*A*C)
				r=(anint(r*1.d15))*1.d-15
				if (r.lt.0.d0) then
					x=-8888 ! error
					print*,'Error!!'
				else
					num=(-B-(r**0.5d0))
					num=(anint(num*1.d15))*1.d-15
					x1=num/(2.d0*A)  ! the lowest solution
					num=(-B+(r**0.5d0))
					num=(anint(num*1.d15))*1.d-15
					x2=num/(2.d0*A)

					if((x1.lt.0.d0).or.(x1.gt.1.d0)) then
						if((x2.lt.0.d0).or.(x2.gt.1.d0)) then
							x=-7777 ! error: x_H out of range
						else
							x=x2
						endif
					else
						x=x1 ! the lowest solution is the stable one
					endif
				endif
			endif
			call dohypo(h,x)
			write(90,'(I8,2F20.16)') iA,x,h
	enddo
	close(90)

	stop
	end



	double precision function Cfun(ip,xf)
c	Calculate \chi function
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


	double precision function DamH(i,alpi,beti)
c	Calculate P(\alpha,\beta) for homogeneous populations
	implicit none
	integer i,alpi,beti
	include "common_homo.h"
	DamH= ymoral(i,2-alpi,2-beti,2-yaction(i,2-alpi,2-beti))
     &		*(1.d0-epsA)	+ epsA*ymoral(i,2-alpi,2-beti,2)
	return
	end


	subroutine num2actionmoral(inum)
c	Write action and moral code in arrays from decimal-binary representation 
c	GG(0) GB(1) BG(2) BB(3) GGC(4) GGD(5) GBC(6) GBD(7)...	
	implicit none
	include "common_homo.h"
	integer inum
	integer inumt,i,j,k
	yaction=-1
	ymoral=-1
	inumt=inum
	do i=1,2
	 do j=1,2
		yaction(1,i,j)=mod(inumt,2)
		inumt=floor(inumt/2.)
	 enddo
	enddo
	do i=1,2
	 do j=1,2
	  do k=1,2
		ymoral(1,i,j,k)=mod(inumt,2)
		inumt=floor(inumt/2.)
	  enddo
	 enddo
	enddo
	return
	end


	subroutine dohypo(hypo,xH)
c     	Compute coherency
      	implicit none
      	double precision hypo,xH
      	include "common_homo.h"
      	integer i,j,k,alp,bet,gam
      	double precision act,mor,Cfun

      	i=1	     
	hypo=0.d0
	do alp=0,1
	do bet=0,1
	do gam=0,1
		    mor=ymoral(i,2-alp,2-bet,2-gam) *(1.d0-2.d0*epsM)  +epsM
		    act=yaction(i,2-alp,2-bet) *(1.d0-epsA)
	
		    hypo=hypo+ (1.d0-abs(mor-Cfun(gam,act))) 
     +		*Cfun(alp,xH)*Cfun(bet,xH)

	enddo
	enddo
	enddo
	hypo=0.5d0*hypo
 
      	return
      	end

