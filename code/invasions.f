
      include "sub_mutation_two.f"

      program mutation_two

      implicit none
      include "common_twoS.h"
      character*80 inpS, inpI, outess
      integer i,j,k,iresident,imutant,ib,ii,jj,kk,equivS,l,
     +				esESS(4096),iactionv(2,2),imoralv(2,2,2)
      double precision xini(4096),hini(4096),ti,
     +			maxdifpayA,mindifpayA,maxdifpayB,mindifpayB,
     +			payoff1(4096),payoff2,difpayA(nmaxb),difpayB(nmaxb),binp

      print*,"Input file with strategies and x's: "
      read(*,*) inpS
      print*,'Output file (ESS): '
      read(*,*) outess
      print*,'b: '
      read(*,*) binp
      print*,'epsilon_A: '
      read(*,*) epsA

      nib=1
      bvect(1)=binp
      esESS=1

      open(10,file=inpS,status='OLD')
      read(10,*)
      do i=1,4096
			read(10,*) ti,xini(i),hini(i)
			if ((ti+1).ne.i) stop 'Error 01'
      enddo
      close(10)

      open(90,file=outess,status='UNKNOWN')
      write(90,'(A1,2A)') '#','   S       W(1|1)/        h1         ',
     + 'GGC GGD GBC GBD BGC BGD BBC BBD   GG  GB  BG  BB '	

      !do i=0,255 ! modify to limit the search to certain strategies and/or invasions
      do j=0,4095 ! modify to limit the search to certain strategies and/or invasions
      do k=0,4095 ! modify to limit the search to certain strategies and/or invasions
		iresident= j ! modify to limit the search to certain strategies and/or invasions
		imutant=  k ! modify to limit the search to certain strategies and/or invasions
		if(imutant.ne.iresident) then
			call doESS(iresident,imutant,
     +			   xini(iresident+1),xini(imutant+1),payoff1(iresident+1),
     +      		   difpayA,difpayB)
		   	if (difpayA(1).lt.1.d-15) then
				if(difpayA(1).gt.-1.d-15) then  ! difpayA=0
					if(difpayB(1).lt.1.d-15) then
						esESS(iresident+1)=0 ! It is not ESS
					endif
				else
					esESS(iresident+1)=0 ! It is not ESS
				endif
			endif
		endif
      enddo !k
      enddo !j
      !enddo !i

      do i=0,4095
	if(esESS(i+1).eq.1) then
	  call num2actionmoral(i,iactionv,imoralv)
	  write(90,'(I6,2x,2F12.8,4x,8I4,2x,4I4,4x,I6)')
     +		 i,payoff1(i+1),hini(i+1),
     +			(((imoralv(ii,jj,kk),kk=1,2),jj=1,2),ii=1,2),
     +			((iactionv(ii,jj),jj=1,2),ii=1,2),
     +				equivS(i)
	endif
      enddo
      close(90)

      stop
      end









