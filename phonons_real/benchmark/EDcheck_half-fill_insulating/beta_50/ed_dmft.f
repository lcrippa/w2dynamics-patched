
      program lisalanc
      implicit none
      integer nss,nsp11,nblockk,ntott,nmpara,nmaxx,Iwmax
      integer Iwmaxreal
      integer ns,nsp1,nblock,ntot,nso2,i,j,k,nmax,nchoos
      integer imaxmu,iteramax,imu,itera,Nitermax
      real*8 deltamu,testdmft,piob,xmu0,chi2,chi22
c      real*8 time,mclock
      real*8 zpart
      integer*4 unlink,a
      parameter (nss=7)
      parameter (nsp11 =nss+1) 
      parameter (nblockk=nsp11*nsp11)
      parameter (ntott=4**nss)
      parameter (nmpara=30)
      include'init.h'  
      parameter (Iwmax= 2**11)  
      parameter (Iwmaxreal =2**12)
      real*8 tpar(nss),epsk(nss),uhub,hmag,xmu,beta,densimp,pi
      real*8 wlow,wup,deltino,range,sumdos
      real*8 om(0:Iwmax), dos(0:Iwmaxreal)
      real*8 omm
      complex*16 ni(-10*Iwmax:10*Iwmax)
      real*8 np(0:Iwmax)
      
      complex*16 omr(0:Iwmaxreal)
      integer ireal
      complex*16 Xi,c0,c1,ci
      complex*16 Gw(0:Iwmax),self(0:Iwmax)
      complex*16 Gww(0:2*Iwmax-1)
      complex*16 Gwreal(0:Iwmaxreal)
      complex*16 sigre(0:Iwmaxreal)
      complex*16 G0w(0:Iwmax),G0wand(0:Iwmax),G0wwand(0:2*Iwmax-1)
c     matrix & friends
      integer ioffset(nblockk+1)
      integer b(nss)
      integer ntable(ntott)
      integer*1 abas(nss)
      integer nleng(nblockk)
      integer*1 imat(nmaxx,nmaxx,nblockk)
      integer idmat(ntott)
      integer iddomat(ntott)
      integer*1 idouble(ntott)
c     for minimization
      real*8 hess(nmpara**2),g(nmpara),xtemp(nmpara),
     &     w(nmpara**2),xprmt(nmpara)
c     for dspev
      real*8 work(3*nmaxx)
      integer info
      real*8 eig(nmaxx),zeig(nmaxx,nmaxx)
      real*8 eigall(ntott)
      real*8 eigold(nmaxx),zold(nmaxx,nmaxx)
      real*8 realmat(nmaxx*(nmaxx+1)/2)
      real*8 xmat(nmaxx,nmaxx)
      real*8 rlehm(nmaxx,nmaxx,nblockk)
      character*80 xyz 

      complex*16 omega, adummy, cdummy1

c      time=mclock()/1
      real*8 Epskold(nss),tparold(nss)

      complex*16 sig0
      real*8 zed,dens_wanted,diffdens,chi,densold,emin,dostest
      real*8 xmularge,xmusmall,denslarge,denssmall,densmix
      integer ifix,inew,iauto
      integer iattempt,idritto,ilarge,ismall
      integer imix,iteraok,iexp
      real*8 threshold,th0,b1,b2,E
      real*8 doublep(ntott)
      real*8 double


c      complex*16 gamma(-Iwmax:Iwmax,-Iwmax:Iwmax,-Iwmax:Iwmax)
c      complex*16 gammach(-Iwmax:Iwmax,-Iwmax:Iwmax,-Iwmax:Iwmax)
c      complex*16 Chi0inv(-Kmax:Kmax-1,-Kmax:Kmax-1,
c     $     -Iwmax:Iwmax-1,-Iwmax:Iwmax-1)


      logical bethe,twodim,symm
     

  

  

      bethe=.true.
      twodim=.false.
      symm=.false.

       if(bethe) then
          write(6,*) 'BETHE lattice case'
       else
          if(twodim) then
             write(6,*) '2D-simple cubic lattice'
          else
             write(6,*) '3D cubic lattice'
          endif
       endif
    

      Pi=dacos(-1.d0)
      Xi=cmplx(0.d0,1.d0)

      
      open (30,file='hubb.dat',form='formatted',status='old')
      open (15,file='hubb.andpar',form='formatted',status='old')
c      open (45,file='varbeta.dat',form='formatted',status='old')

      call datain(nss,uhub,xmu,hmag,ns,beta,wlow,wup,deltino,
     $     imaxmu,deltamu,iteramax,testdmft,dens_wanted,ifix,
     $     inew,iexp,th0,iauto)
      
      piob=Pi/beta
      nsp1=ns+1
      nblock=nsp1*nsp1
      ntot=4**ns
      nso2=ns/2
      

      b(1)=1
      do i=2,ns
         b(i)=b(i-1)*4
      enddo
     

      range=(wup-wlow)/dfloat(Iwmaxreal)
      do i=0,Iwmaxreal
         omr(i)=wlow+range*dfloat(i)+Xi*deltino
      enddo

      nmax=nchoos(ns,nso2)*nchoos(ns,nso2)

      write(6,*) 'nmax', nmax, 'nmaxx', nmaxx

      if (nmax.ne.nmaxx) then
         write(6,*)'You are running serious troubles'
         write(6,*)'put nmaxx =',nmax
         stop
      endif

      write(6,*)'U  = ',uhub
      write(6,*)'mu = ',xmu
      write(6,*)'h = ',hmag
      write(6,*)'beta =',beta
      write(6,*) 'Iwmax =',Iwmax
      write(6,*)'for real frequencies: ',wlow,wup
c      write(6,*)'compute',imaxmu,' mu values'
      call flush(6) 
      if (ifix.eq.0) then
         write(6,*)'Working at fixed mu =',xmu
      elseif(ifix.eq.1) then
         write(6,*)'Working at fixed density =',dens_wanted  
         if(beta.ge.200.d0) then
            chi= 1.2d0/(1+0.6d0*dabs(uhub))
            xmu=(dens_wanted-1.d0)/chi
            xmu=0.5d0*uhub+xmu
         endif
      else
         write(6,*)'ifix should be 0, 1 !!'
         stop
      endif
      chi= 1.2d0/(1+0.6d0*dabs(uhub))
      call flush(6)
      write(6,*)'Initial mu =',xmu
      open (34,file='self-en_wim',form='formatted',status='unknown')
      open (90,file='gm_wim',form='formatted',status='unknown')
      open (91,file='gm_wre',form='formatted',status='unknown')
      open (92,file='g0m',form='formatted',status='unknown')
      open (93,file='g0mand',form='formatted',status='unknown')
      open (99,file='gw_cut',form='formatted',status='unknown')
      open (79,file='vert_chi',form='formatted',status='unknown')
      open(35,file='vsmu',form='formatted',status='unknown')
      open(89,file='GAMMA_WIM',form='formatted',status='unknown')
      open(85,file='SELF_LOC',form='formatted',status='unknown')
      open(50,file='SELF_Q',form='formatted',status='unknown')
      open(51,file='chiPI_wim',form='formatted',status='unknown')
      open(52,file='chiPiPi',form='formatted',status='unknown')
      open(74,file='chiSC_test',form='formatted',status='unknown')
      open(73,file='chiSC_wim',form='formatted',status='unknown')
      open(45,file='anzahl',form='formatted',status='unknown')
c      xmu0=xmu
c      write(6,*)'mu0=',xmu0

c     costruzione delle matrici varie

      call computematrix(ioffset,nblock,ntot,b,ns,abas,nsp1,nmax,
     $     nleng,imat,idmat,iddomat,ntable,idouble)

      write(6,*)'# of independent blocks:',nblock

      do i=0,Iwmax
         om(i)=(2.d0*dfloat(i)+1.d0)*piob
      enddo
      
      do i=-10*Iwmax,10*Iwmax
         ni(i)=dfloat(i)*piob*Xi
c        write(6,*) i, dimag(ni(i))
      enddo

      do i=0,Iwmax
         np(i)=(2.d0*dfloat(i)+1.d0)*piob
      enddo

      do imu=0,imaxmu
         xmu=xmu+deltamu*dfloat(imu)

         call initial(epsk,tpar,ns,xmu)
         write(6,*)'------------------------------------------------'
         write(6,*) 'starting Anderson parameters '
         write(6,*)'------------------------------------------------'
         write(6,*) 'Eps(k) '
         do i=2,ns
           write(6,'(2f27.16)')epsk(i)
         enddo
         write(6,*)'V(k) '
         do i=1,ns-1
           write(6,'(2f27.16)')tpar(i)
         enddo
c        inew=0 leggi dal file
c        inew=1 calcola xmu da scratch
         if (ifix.eq.1.and.inew.eq.0) read(15,*)xmu
         if (inew.eq.0.and.iauto.eq.0) read(15,*)xmu
         write(6,'(f27.16,"   #chemical potential")')xmu
c        write(6,*)'------------------------------------------------'
         do i=2,ns
            write(86,*)epsk(i)
         enddo
c        write(6,*)'the real xmu is ',xmu
        
         iteraok=1

         do itera=1,iteramax
            imix=0
            threshold=th0/dfloat(iteraok)**iexp
            if (threshold.lt.1.d-6) threshold=1.d-6
            epsk(1)=-xmu
c            write(6,*)'epsk=',epsk(1)
            do i=0,iwmax
               call calcg0(om(i),cdummy1,tpar,epsk,ns,Xi,Iwmax)
                  g0wand(i)=cdummy1
            end do
            
c           qui esce con G0and (non alla -1)

            write(6,*)'------------------------------------------------'
            write(6,*)'   Iteration : ',itera
            write(6,*)'------------------------------------------------'
c            time=mclock()/100.d0
            write(6,*)'Threshold for n=',threshold
            if (ifix.ne.0) write(6,*)'initial mu = ',xmu
            call flush(6)

            iattempt=0
            idritto=0
            ilarge=0
            ismall=0
            double=0.d0

            call diag(tpar,epsk,uhub,hmag,ioffset,nleng,idmat,
     $           imat,nblock,ns,nsp1,nmax,ntable,work,eig,zeig,
     $           eigall,abas,
     $           ntot,realmat,beta,Xi,piob,range,wup,wlow,Iwmax,Gw,
     $           om,eigold,zold,b,rlehm,xmat,nmaxx,zpart,densimp,
     $           xmu,emin,idouble,doublep,double)

c            call computegreal(Gwreal,Iwmaxreal,rlehm,omr,eigall,beta,
c     $           ioffset,nleng,nblock,ntot,nsp1,pi,nmax,zpart,
c     $           densimp)

            write(6,*)'densimp = ',densimp
            
            if (ifix.eq.1) then 
 67            diffdens=densimp-dens_wanted
               write(6,*)'<n>=',densimp,'  <n>-n0=',diffdens

               if (dabs(diffdens).lt.threshold) then               
                  write(6,*)'------------------------'
                  write(6,*)'converged at step ',itera
                  write(6,*)'------------------------'
                  call flush(6)
c           converged now
                  iteraok=iteraok+1
               else
c           not converged now
                  write(6,*)'-Not yet converged-'
                  call flush(6)
                  if (diffdens.gt.0.d0) then
                     ilarge=1
                     xmularge=xmu
                     denslarge=densimp
                  endif
                  if (diffdens.lt.0.d0) then
                     ismall=1
                     xmusmall=xmu
                     denssmall=densimp
                  endif
                  if (ilarge*ismall.eq.0) then
                     write(6,*)'Still lacking a bracket'
                     xmu=xmu-chi*diffdens
                     write(6,*)'Try with xmu=',xmu
                     idritto=idritto+1
                     if (idritto.eq.6) then
                        chi=chi*2.d0
                     endif
                     call diag(tpar,epsk,uhub,hmag,ioffset,nleng,idmat,
     $                    imat,nblock,ns,nsp1,nmax,ntable,work,eig,
     $                    zeig,eigall,abas, ntot,realmat,
     $                     beta,Xi,piob,range,wup,wlow,Iwmax,Gw,
     $                    om,eigold,zold,b,rlehm,xmat,nmaxx,zpart,
     $                    densimp,xmu,emin,idouble,doublep,double)

c                     call computegreal(Gwreal,Iwmaxreal,rlehm,omr,
c     $                    eigall,beta,ioffset,nleng,nblock,ntot,
c     $                    nsp1,pi,nmax,zpart,densimp)

                     go to 67
                  else
                     write(6,*)'Interpolation between:'
                     write(6,*)'mu, n (1):',xmusmall,denssmall
                     write(6,*)'mu, n (2):',xmularge,denslarge
                     iattempt=iattempt+1
                     if (iattempt.gt.10) then 
                        imix=1
                        write(6,*)'Interpolation not converged: Mix!'

                        b1=(dens_wanted-denslarge)/
     $                       (denssmall-dens_wanted)
                        b2=(denssmall-dens_wanted)/
     $                       (denssmall-denslarge)     

                        write(6,*)'Coeff(1):',b2*b1
                        write(6,*)'Coeff(2):',b2
                        call diag(tpar,epsk,uhub,hmag,ioffset,nleng,
     $                       idmat,imat,nblock,ns,nsp1,nmax,ntable,
     $                       work,eig,zeig,eigall,abas, ntot,realmat,
     $                       beta,Xi,piob,range,wup,wlow,Iwmax,Gw,
     $                       om,eigold,zold,b,rlehm,xmat,nmaxx,zpart,
     $                       densimp,xmusmall,emin,idouble,doublep,
     $                       double)
c                        call computegreal(Gwreal,Iwmaxreal,rlehm,omr,
c     $                       eigall,beta,ioffset,nleng,nblock,ntot,
c     $                       nsp1,pi,nmax,zpart,densimp)
                        do i=0,Iwmax
                           Gw(i)=b1*Gw(i)
                        enddo
                        do i=0,Iwmaxreal
                           Gwreal(i)=b1*Gwreal(i)
                        enddo
                        densmix=b1*densimp
                        call diag(tpar,epsk,uhub,hmag,ioffset,nleng,
     $                       idmat,imat,nblock,ns,nsp1,nmax,ntable,
     $                       work,eig,zeig,eigall,abas, ntot,realmat,
     $                       beta,Xi,piob,range,wup,wlow,Iwmax,Gw,
     $                       om,eigold,zold,b,rlehm,xmat,nmaxx,zpart,
     $                       densimp,xmularge,emin,idouble,doublep,
     $                       double)
c                        call computegreal(Gwreal,Iwmaxreal,rlehm,omr,
c     $                       eigall,beta,ioffset,nleng,nblock,ntot,
c     $                       nsp1,pi,nmax,zpart,densimp)
                        do i=0,Iwmax
                           Gw(i)=b2*Gw(i)
                        enddo 
                        do i=0,Iwmaxreal
                           Gwreal(i)=b2*Gwreal(i)
                        enddo 
                        densmix=b2*(densmix+densimp)
                        write(6,*)'densmix =',densmix
                        densimp=densmix
c     exit with some "reasonable" unique xmu for next iteration
                        xmu=xmusmall+(xmularge-xmusmall)*
     $                       (denssmall-dens_wanted)/
     $                       (denssmall-denslarge)
                        go to 68
                     else
                        xmu=xmusmall+(xmularge-xmusmall)*
     $                       (denssmall-dens_wanted)/
     $                       (denssmall-denslarge)
                        write(6,*)'interpolated xmu',xmu
                        call diag(tpar,epsk,uhub,hmag,ioffset,nleng,
     $                       idmat,imat,nblock,ns,nsp1,nmax,ntable,
     $                       work,eig,zeig,eigall,abas, ntot,realmat,
     $                       beta,Xi,piob,range,wup,wlow,Iwmax,Gw,
     $                       om,eigold,zold,b,rlehm,xmat,nmaxx,zpart,
     $                       densimp,xmu,emin,idouble,doublep,double)
c                        call computegreal(Gwreal,Iwmaxreal,rlehm,omr,
c     $                       eigall,beta,ioffset,nleng,nblock,ntot,
c     $                       nsp1,pi,nmax,zpart,densimp)
                        goto 67
                     endif
                  endif
               endif
            endif
          
c            time=mclock()/100.d0-time
c            write(6,*)'Tempo in diag:',time
            

 68         write(6,*) 'SELF_CONSISTENCY LOOP'
            if(bethe) then
               write(6,*) 'Bethe lattice' 
               do i=0,Iwmax
                  adummy=om(i)*Xi-Epsk(1)-Gw(i)/4.d0
                  if (iauto.eq.0)then
                     G0w(i)=1.d0/adummy
                  else
                     G0w(i) = (0.5d0/adummy+0.5d0*G0wand(i))
                  endif
               enddo 
            else
               if(twodim) then
                   write(6,*) 'for the 2D-cubic' 
                  call selfconst2(g0wand,Gw,G0w,Xi,Pi,Iwmax,om,xmu,beta)
               else  
                  write(6,*) 'for the 3D-cubic'
                  call selfconst3(g0wand,Gw,G0w,Xi,Pi,Iwmax,om,xmu)
               endif
            endif
            

               if (iauto.eq.0) then 
                  write(6,*)'Average Double Occupancy :',double
                  goto 777
               endif
           
            
            Nitermax=4000
            write(6,*)'minimization'
            call flush(6)


            rewind(62)
            do i=0,Iwmax
               write(62,'(3f17.10)') om(i),dreal(G0w(i)),dimag(G0w(i))
               call flush(62)
            enddo
           
c           time=mclock()/100.d0
            call search(chi2,Nitermax,hess,g,xtemp,w,xprmt,nmpara,
     $           tpar,epsk,ns,piob,Xi,Iwmax,G0w,G0wand,om)
c           write(6,*)'Truncation error :',chi2
            call flush(6)
            rewind(15)
c            time=mclock()/100.d0-time
c            write(6,*)'Tempo in Search=',time


            rewind(63)
            do i=0,Iwmax
               write(63,'(3f17.10)') om(i),dreal(G0wand(i)),
     $              dimag(G0wand(i))
               call flush(63)
            enddo


            call rheader(15)
            
            read(15,*)

            do i=2,ns
               read(15,*)Epskold(i)
            end do

            read(15,*)

            do i=1,ns-1
               read(15,*)tparold(i)
            end do
            
            chi22=0.d0

            do i=1,ns-1
               chi22=chi22+(tpar(i)-tparold(i))**2
               chi22=chi22+(Epsk(i+1)-Epskold(i+1))**2
            enddo

            chi22=chi22/(2*ns-2)
            
            write(6,*)'------------------------------------------------'
            write(6,*) 'new Anderson parameters '
            write(6,*)'------------------------------------------------'
            write(6,*) 'Eps(k) '
            do i=2,ns
               write(6,'(2f27.16)')epsk(i)
            enddo
             write(6,*)'V(k) '
            do i=1,ns-1
               write(6,'(2f27.16)')tpar(i)
            enddo
            write(6,'(f27.16,"   #chemical potential")')xmu
            write(6,*)'------------------------------------------------'

            write(6,'(" convergence parameter :",e18.10)')chi22
            call flush(6)            

            write(22,*)itera,chi22
            call flush(22) 
            rewind(15)
            call wheader(15,ns,Iwmax)
            write(15,'(9a)')' Eps(k)'
            do i=2,ns
               write(15,*)Epsk(i)
            end do
            write(15,'(9a)')' tpar(k)'
            do i=1,ns-1
               write(15,*)tpar(i)
            end do
            write(15,*)xmu,'    #chemical potential'
            call flush(15)
c            rewind(80)
c            write(80)ns,nblock,ntot,nmax,zpart,eigall,ioffset,
c     $           nleng,rlehm,epsk,tpar


            if (chi22.lt.testdmft) goto 777
         enddo
C========+=========+=========+=========+=========+=========+=========+=$
c                               OUTPUT
C========+=========+=========+=========+=========+=========+=========+=$
         
777      continue


c     compute real-frequency Green's function
        if (iauto.ne.0) then
           call diag(tpar,epsk,uhub,hmag,ioffset,nleng,idmat,
     $          imat,nblock,ns,nsp1,nmax,ntable,work,eig,
     $          zeig,eigall,abas, ntot,realmat,
     $          beta,Xi,piob,range,wup,wlow,Iwmax,Gw,
     $          om,eigold,zold,b,rlehm,xmat,nmaxx,zpart,
     $          densimp,xmu,emin,idouble,doublep,double)

           
           

           
           do i=0,Iwmax-1
              Gww(i)=Gw(i)
              Gww(i+Iwmax)=1.d0/ni(2*(Iwmax+i)+1)
           enddo
        
           
           do i=0,2*Iwmax-1
              omm=dimag(ni(2*i+1))
               call calcg0(omm,cdummy1,tpar,epsk,ns,Xi,Iwmax)
               cdummy1=1.d0/cdummy1
               g0wwand(i)=cdummy1
c               write(6,*)'Check Extended Green',i,G0wwand(i),Gww(i)
            end do


            call computegreal(Gwreal,Iwmaxreal,rlehm,omr,eigall,beta,
     $          ioffset,nleng,nblock,ntot,nsp1,pi,nmax,zpart,densold)

            write (6,*) 'Now computing Matsubara self-energy'       
            write(6,*) 'Pi is equal to ', Pi

            do i=0, Iwmax
               self(i) =g0wwand(i)-1.d0/Gw(i)
               write(34,'(3f17.10)') om(i),dreal(self(i)),dimag(self(i))
            enddo

            write (6,*) 'Now computing self-energy on the real axis'       
           
            
            do i=0,Iwmaxreal
               sigre(i)=(0.d0,0.d0)
               do j=1,ns-1
                  sigre(i)=sigre(i)-tpar(j)**2/(omr(i)-Epsk(1+j))
               end do

c     Note the -PI factor coming from the normalization previously used to get
C the DOS
c           sigupre(i)=Pi*gwreal(i)
c           sigdore(i)=Pi*gwrealdo(i)
               write(22,'(4f17.10)') dreal(Gwreal(i)), dimag(Gwreal(i)),
     $              dreal(1.d0/Gwreal(i)),dimag(1.d0/Gwreal(i))
               sigre(i)=sigre(i)+omr(i)-Epsk(1)+
     $              dconjg(Gwreal(i))/(Pi*Gwreal(i)*dconjg(Gwreal(i)))
            enddo



            write(6,*)'Average Double Occupancy :',double
        endif


c        write (6,*) 'Now computing self-energy'        

c        call super(G0w,Gw,Iwmax,pi,beta,xmu)
     
        write(6,*)'Average Double Occupancy :',double

        if(twodim) then
           open(57,file='dos.dat',form='formatted',status='unknown')
           open(58,file='dos0.dat',form='formatted',status='unknown')

           call dos2(sigre,Iwmaxreal,omr,Pi,xmu,beta)   
        endif
        
        if(bethe) then
           open(57,file='dos.dat',form='formatted',status='unknown')
c           open(58,file='dos0.dat',form='formatted',status='unknown')
           
c           sumdos=0.d0
           do i=0,Iwmaxreal
              dos(i)=0.d0
              do j=0,4999
                 E=-1.d0+2.d0*dfloat(j)/dfloat(5000) +1.d0/dfloat(5000)
                 dos(i)=dos(i)-2.d0/(pi)**2*dsqrt(1.d0-E**2)*
     $                dimag(1.d0/(omr(i)-Epsk(1)-E-sigre(i)))
     $                /dfloat(5000)
              enddo
c              sumdos=sumdos+dos(i)*
              write(57,'(2f17.10)') dreal(omr(i)),dos(i)
           enddo
        endif
        
        
        
c        do i=0,Iwmaxreal
c           write(6,*) i, dos(i)
c           write(57,'(2f17.10)') dreal(omr(i)),dos(i)
c        enddo

         write(6,*)'------------------------------------------------'
         write(6,*) 'final Anderson parameters '
         write(6,*)'------------------------------------------------'
         write(6,*)'cut from the following line...'
         write(6,*) 'Eps(k) '
         do i=2,ns
           write(6,'(2f27.16)')epsk(i)
         enddo
           write(6,*)'V(k) '
         do i=1,ns-1
            write(6,'(2f27.16)')tpar(i)
         enddo
         write(6,'(f27.16,"   #chemical potential")')xmu      
         write(6,*)'...to the line above and copy into hubb.andpar'
         write(6,*)'------------------------------------------------'
 
         do i=0,2*Iwmax-1
            write(90,'(3f27.16)') dimag(ni(2*i+1))
     $           ,dreal(Gww(i)),dimag(Gww(i))
         enddo
         
         do i=0,Iwmax
            G0w(i)=1.d0/G0w(i)
            write(92,'(3f17.10)') dimag(ni(2*i+1))
     $           ,dreal(G0w(i)),dimag(G0w(i))
         end do
               
         do i=0,2*Iwmax-1
            write(93,'(3f17.10)')dimag(ni(2*i+1)),dreal(G0wwand(i)),
     $           dimag(G0wwand(i))
         end do
         
         do i=0,Iwmaxreal
            write(91,'(5f17.10)')dreal(omr(i)),
     $           dreal(Gwreal(i)),dimag(Gwreal(i)),dreal(sigre(i)),
     $           dimag(sigre(i))
         end do
         
c     calcolo il residuo di quasiparticella (zed)
         sig0=g0w(0)-1.d0/gw(0)
         write(6,*) 'sig0 = ',sig0
         zed=1.d0-dimag(sig0)/piob
         zed=1.d0/zed
         write(6,*) 'zed = ', zed
         call flush(6)
          
         dostest=0.d0
         do i=1,ns-1
            dostest=dostest+tpar(i)**2
         enddo
         write(6,*) 'DOSTEST=', dostest

c         call dos2(sigre,Iwmaxreal,omr,Pi)
         
         write(35,888) uhub,xmu,densimp,densold,emin,zpart
888      format(f6.3,1x,f6.3,1x,f12.8,1x,f12.8,1x,f12.8,1x,f12.8)
         write(45,*) uhub,densimp, double
         enddo
C========+=========+=========+=========+=========+=========+=========+=$
c     END OF Iteration 
C========+=========+=========+=========+=========+=========+=========+=$
      write(6,'(a20)')'     End lisalanc '
      write(6,'(a60)')'========================================'
      end
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine initial(epsk,tpar,ns,xmu)
      implicit none
      integer ns,i
      real*8 epsk(ns),tpar(ns),xmu

      rewind(15)
      call rheader(15)
c      Epsk(1)=-xmu
      read(15,*)
      do i=2,ns
         read(15,*)Epsk(i)
      end do
      read(15,*)
      do i=1,ns-1
         read(15,*)tpar(i)
      end do
      end
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine rheader(k)
      implicit none
      integer i,k
 
      do 1 i=1,8
         read(k,*)
1     continue
      end 
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine wheader(k,ns,Iwmax)
      implicit none
      integer i,k
      integer ns,Iwmax

      character *80 xyz
      character *8 xxx
      xxx=' 1-band '
      write(k,'(a55)')'========================================'
      write(k,'(a25,a30)')xxx,'30-Sep-95 LANCZOS  '
      write(k,'(a55)')'========================================'
      rewind(30)
      write(k,'(4(a6,I4))') 'NSITE ',ns,'IWMAX',iwmax
      read(30,'(a60)')xyz
      read(30,'(a60)')xyz
      read(30,'(a60)')xyz
      do 3 i=1,4
      read(30,'(a60)')xyz
      write(k,'(a60)')xyz
3     continue 
      rewind(30)
      end
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine search(fmin,Nitermax,hess,g,xtemp,w,xprmt,nmpara,
     $     tpar,epsk,ns,piob,Xi,Iwmax,G0w,G0wand,om)
      implicit none
      integer nmpara,ns,Iwmax,nbparm,icount,i,mode,Nitermax,iexit,nsym
      real*8 dfn,deps,fmin,hh
      real*8 hess(nmpara**2),g(nmpara),xtemp(nmpara),
     &     w(nmpara**2),xprmt(nmpara)
      real*8 tpar(ns),epsk(ns),piob  
      real*8 om(0:Iwmax)
      complex*16 Xi
      complex*16 G0w(0:Iwmax),G0wand(0:Iwmax)
      logical symm
      external difference
      data iexit/0/
      data hh/1.e-5/
      
      common /symm/ symm
 
c     number of parameters to be optimized (per spin):
c     eps(2)...... eps(ns) --->  ns-1
c     tpar(1)........ tpar(ns-1) --->  ns-1

      if(symm) then
         nsym=(ns-1)/2
         nbparm=2*(nsym-1)
         if (nbparm.gt.nmpara) stop ' nmpara too small'
C
         icount=0
     
c         nsym=(ns-1)/2
         nbparm=2*(nsym)
         do i=2,nsym+1
            icount=icount+1	
            xtemp(icount)=Epsk(i)
 11         continue
         end do
         
         do i=1,nsym
            icount=icount+1
            xtemp(icount)=tpar(i)
 12         continue
         end do
         
         do i=nbparm+1,nmpara
            xtemp(i)=0.d0
         end do
      
         do i=1,nbparm
            xprmt(i)=dabs(xtemp(i))+1.d-15
         end do
      else  
         nbparm=2*(ns-1)
         if (nbparm.gt.nmpara) stop ' nmpara too small'
C     
         icount=0
         do i=2,ns
            icount=icount+1	
            xtemp(icount)=Epsk(i)
 111        continue
         end do
         
         do i=1,ns-1
            icount=icount+1
            xtemp(icount)=tpar(i)
 112         continue
         end do
         
         do i=nbparm+1,nmpara
            xtemp(i)=0.d0
         end do
         
         do i=1,nbparm
            xprmt(i)=dabs(xtemp(i))+1.d-15
         end do
      endif

      

      mode=1
c     use va10 for search
      dfn=-.5d0
      deps=.00001d0

      call minimize(difference,nbparm,xtemp,fmin,g,hess,w
     +     ,dfn,xprmt,hh,deps,mode,Nitermax,iexit,
     $     tpar,epsk,piob,Iwmax,Xi,ns,G0wand,G0w,om)
      write (6,30) iexit,fmin
c     do i=1,10 
c     write(73,*)xtemp(i),xprmt(i)
c     enddo
 30   format(' iexit fmin ',i5,e14.6)
      end
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine minimize (funct, n, x, f, g, h, w, dfn, xm,
     $  hh, eps, mode, maxfn, iexit,tpar,epsk,piob,
     $     Iwmax,Xi,ns,G0wand,G0w,om)
      implicit none
      integer ns,Iwmax
      integer np,n,n1,nn,is,iu,iv,ib,idiff,iexit,mode,ij,maxfn,i,j,
     $     i1,jk,ik,k,itn,ifn,link,int
      real*8 z,zz,dmin,f,df,dfn,aeps,eps,alpha,ff,tot,f1,f2,half,
     $     gys,dgs,sig,hh,gs0
      real*8 tpar(ns),epsk(ns),piob
      complex*16 Xi
      complex*16 G0w(0:Iwmax),G0wand(0:Iwmax)
      real*8  x(*), g(*), h(*), w(*), xm(*)
      real*8 om(0:Iwmax)
      external funct
      data half /0.5d0/


      np = n + 1
      n1 = n - 1
      nn=(n*np)/2
      is = n
      iu = n
      iv = n + n
      ib = iv + n
      idiff = 1
      iexit = 0
      if (mode .eq. 3) go to 15
      if (mode .eq. 2) go to 10
      ij = nn + 1
      do 5 i = 1, n
      do 6 j = 1, i
      ij = ij - 1
   6  h(ij) = 0.d0
   5  h(ij) = 1.d0
      go to 15
  10  continue
      ij = 1
      do 11 i = 2, n
      z = h(ij)
      if (z .le. 0.d0) return
      ij = ij + 1
      i1 = ij
      do 11 j = i, n
      zz = h(ij)
      h(ij) = h(ij) / z
      jk = ij
      ik = i1
      do 12 k = i, j
      jk = jk + np - k
      h(jk) = h(jk) - h(ik) * zz
      ik = ik + 1
  12  continue
      ij = ij + 1
  11  continue
      if (h(ij) .le. 0.d0) return
  15  continue
      ij = np
      dmin = h(1)
      do 16 i = 2, n
      if (h(ij) .ge. dmin) go to 16
      dmin = h(ij)
  16  ij = ij + np - i
      if (dmin .le. 0.d0) return
      z = f
      itn = 0
      call funct (n, x, f,epsk,tpar,ns,piob,Xi,Iwmax,
     $     G0w,G0wand,om)
      ifn = 1
      df = dfn
      if (dfn .eq. 0.d0) df = f - z
      if (dfn .lt. 0.d0) df = abs (df * f)
      if (df .le. 0.d0) df = 1.d0
  17  continue
      do 19 i = 1, n
      w(i) = x(i)
  19  continue
      link = 1
      if (idiff - 1) 100, 100, 110
  18  continue
      if (ifn .ge. maxfn) go to 90
  20  continue
  21  continue
      itn = itn + 1
      w(1) = -g(1)
      do 22 i = 2, n
      ij = i
      i1 = i - 1
      z = -g(i)
      do 23 j = 1, i1
      z = z - h(ij) * w(j)
      ij = ij + n - j
  23  continue
  22  w(i) = z
      w(is+n) = w(n) / h(nn)
      ij = nn
      do 25 i = 1, n1
      ij = ij - 1
      z = 0.d0
      do 26 j = 1, i
      z = z + h(ij) * w(is+np-j)
      ij = ij - 1
  26  continue
  25  w(is+n-i) = w(n-i) / h(ij) - z
      z = 0.d0
      gs0 = 0.d0
      do 29 i = 1, n
      if (z * xm(i) .ge. abs (w(is+i))) go to 28
      z = abs (w(is+i)) / xm(i)
  28  gs0 = gs0 + g(i) * w(is+i)
  29  continue
      aeps = eps / z
      iexit = 2
      if (gs0 .ge. 0.d0) go to 92
      alpha = -2.d0 * df / gs0
      if (alpha .gt. 1.d0) alpha = 1.d0
      ff = f
      tot = 0.d0
      int = 0
      iexit = 1
  30  continue
      if (ifn .ge. maxfn) go to 90
      do 31 i = 1, n
      w(i) = x(i) + alpha * w(is+i)
  31  continue
      call funct (n, w, f1,epsk,tpar,ns,piob,Xi,Iwmax,
     $     G0w,G0wand,om)
      ifn = ifn + 1
      if (f1 .ge. f) go to 40
      f2 = f
      tot = tot + alpha
  32  continue
      do 33 i = 1, n
      x(i) = w(i)
  33  continue
      f = f1
      if (int - 1) 35, 49, 50
  35  continue
      if (ifn .ge. maxfn) go to 90
      do 34 i = 1, n
      w(i) = x(i) + alpha * w(is+i)
  34  continue
      call funct (n, w, f1,epsk,tpar,ns,piob,Xi,Iwmax,
     $     G0w,G0wand,om)
      ifn = ifn + 1
      if (f1 .ge. f) go to 50
      if ((f1 + f2 .ge. f + f) .and.
     $  (7.0d0 * f1 + 5.0d0 * f2 .gt. 12.0d0 * f)) int = 2
      tot = tot + alpha
      alpha = 2.d0 * alpha
      go to 32
  40  continue
      if (alpha .lt. aeps) go to 92
      if (ifn .ge. maxfn) go to 90
      alpha = half * alpha
      do 41 i = 1, n
      w(i) = x(i) + alpha * w(is+i)
  41  continue
      call funct (n, w, f2,epsk,tpar,ns,piob,Xi,Iwmax, 
     $     G0w,G0wand,om)
      ifn = ifn + 1
      if (f2 .ge. f) go to 45
      tot = tot + alpha
      f = f2
      do 42 i = 1, n
      x(i) = w(i)
  42  continue
      go to 49
  45  continue
      z = 0.1d0
      if (f1 + f .gt. f2 + f2)
     $  z = 1.d0 + half * (f - f1) / (f + f1 - f2 - f2)
      if (z .lt. 0.1d0) z = 0.1d0
      alpha = z * alpha
      int = 1
      go to 30
  49  continue
      if (tot .lt. aeps) go to 92
  50  continue
      alpha = tot
      do 56 i = 1, n
      w(i) = x(i)
      w(ib+i) = g(i)
  56  continue
      link = 2
      if (idiff - 1) 100, 100, 110
  54  continue
      if (ifn .ge. maxfn) go to 90
      gys = 0.d0
      do 55 i = 1, n
      w(i) = w(ib+i)
      gys = gys + g(i) * w(is+i)
  55  continue
      df = ff - f
      dgs = gys - gs0
      if (dgs .le. 0.d0) go to 20
      link = 1
      if (dgs + alpha * gs0 .gt. 0.d0) go to 52
      do 51 i = 1, n
      w(iu + i) = g(i) - w(i)
  51  continue
      sig = 1.d0 / (alpha * dgs)
      go to 70
  52  continue
      zz = alpha / (dgs - alpha * gs0)
      z = dgs * zz - 1.d0
      do 53 i = 1, n
      w(iu+i) = z * w(i) + g(i)
  53  continue
      sig = 1.d0 / (zz * dgs * dgs)
      go to 70
  60  continue
      link = 2
      do 61 i = 1, n
      w(iu+i) = w(i)
  61  continue
      if (dgs + alpha * gs0 .gt. 0.d0) go to 62
      sig = 1.d0 / gs0
      go to 70
  62  continue
      sig = -zz
  70  continue
      w(iv+1) = w(iu+1)
      do 71 i = 2, n
      ij = i
      i1 = i - 1
      z = w(iu+i)
      do 72 j = 1, i1
      z = z - h(ij) * w(iv+j)
      ij = ij + n - j
  72  continue
      w(iv+i) = z
  71  continue
      ij = 1
      do 75 i = 1, n
      z = h(ij) + sig * w(iv+i) * w(iv+i)
      if (z .le. 0.d0) z = dmin
      if (z .lt. dmin) dmin = z
      h(ij) = z
      w(ib+i) = w(iv+i) * sig / z
      sig = sig - w(ib+i) * w(ib+i) * z
      ij = ij + np - i
  75  continue
      ij = 1
      do 80 i = 1, n1
      ij = ij + 1
      i1 = i + 1
      do 80 j = i1, n
      w(iu+j) = w(iu+j) - h(ij) * w(iv+i)
      h(ij) = h(ij) + w(ib+i) * w(iu+j)
      ij = ij + 1
  80  continue
      go to (60, 20), link
  90  continue
      iexit = 3
      go to 94
  92  continue
      if (idiff .eq. 2) go to 94
      idiff = 2
      go to 17
  94  continue
      return
 100  continue
      do 101 i = 1, n
         z = hh * xm(i)
         w(i) = w(i) + z
         call funct (n, w, f1,epsk,tpar,ns,piob,Xi,Iwmax, 
     $        G0w,G0wand,om)
      g(i) = (f1 - f) / z
      w(i) = w(i) - z
 101  continue
      ifn = ifn + n
      go to (18, 54), link
 110  continue
      do 111 i = 1, n
      z = hh * xm(i)
      w(i) = w(i) + z
      call funct (n, w, f1,epsk,tpar,ns,piob,Xi,Iwmax , 
     $     G0w,G0wand,om)
      w(i) = w(i) - z - z
      call funct (n, w, f2, epsk,tpar,ns,piob,Xi,Iwmax, 
     $     G0w,G0wand,om)
      g(i) = (f1 - f2) / (2.d0 * z)
      w(i) = w(i) + z
 111  continue
      ifn = ifn + n + n
      go to (18, 54), link
      end 
C========+=========+=========+=========+=========+=========+=========+=$
        subroutine difference(nbparm,x,f,epsk,tpar,ns,piob,Xi,
     $     Iwmax,G0w,G0wand,om)
        implicit none
        integer ns,Iwmax,nbparm,icount,i,nsym
        real*8  tpar(ns),epsk(ns),piob,diff,f
        complex*16  cdummy1,G0w(0:Iwmax),G0wand(0:Iwmax)
        complex*16 Xi
        real*8  x(nbparm),om(0:Iwmax),norm
        logical symm
        common /symm/ symm

        icount=1

        if(symm) then
           nsym=(ns-1)/2
           do i=1,nsym
              icount=icount+1
              Epsk(icount)=x(i)
           end do
           icount=0
           do i=nsym+1,2*nsym
              icount=icount+1
              tpar(icount)=x(i)
           end do
        
           do i=nsym+2,ns
              Epsk(i)=-Epsk(i-nsym)
              tpar(i-1)=tpar(i-nsym-1)
           enddo 

           norm=0.d0
           
           do i=1,ns-1
              norm=norm+tpar(i)**2
           enddo
           norm=dsqrt(0.25d0/norm)

           do i=1,ns-1
              tpar(i)=norm*tpar(i)
           enddo
           

           
        else
      
           do i=1,ns-1
              icount=icount+1
              Epsk(icount)=x(i)
           end do
           icount=0
           do i=ns,2*(ns-1)
              icount=icount+1
              tpar(icount)=x(i)
           end do

           norm=0.d0

           do i=1,ns-1
              norm=norm+tpar(i)**2
           enddo
           norm=dsqrt(0.25d0/norm)

           do i=1,ns-1
              tpar(i)=norm*tpar(i)
           enddo

        endif


        diff=0.d0
        do i=0,Iwmax
           call calcg0(om(i),cdummy1,tpar,epsk,ns,Xi,Iwmax)
           g0wand(i)=cdummy1
c           diff = diff + (abs(g0w(i)-g0wand(i))/om(i))**2
           diff = diff + abs(g0w(i)-g0wand(i))/dfloat(i+1)
c           diff = diff + abs(g0w(i)-g0wand(i))!/dfloat(i+1)
        end do
        f=diff/dfloat(Iwmax+1)

        if(symm) then
           icount=1
           do i=1,nsym
              icount=icount+1
              x(i)=Epsk(icount)
           end do
           icount=0
           do i=nsym+1,2*nsym
              icount=icount+1
              x(i)=tpar(icount)
           end do
        else
           icount=1
           do i=1,ns-1
              icount=icount+1
              x(i)=Epsk(icount)
           enddo
           icount=0
           do i=ns,2*(ns-1)
              icount=icount+1
              x(i)=tpar(icount)
           enddo
        endif

        end
C========+=========+=========+=========+=========+=========+=========+=$
        subroutine calcg0(omega,g0and,tpar,epsk,ns,Xi,Iwmax)
        implicit none
        integer ns,Iwmax,i
        real*8 tpar(ns),epsk(ns)
        complex*16 g0and,Xi
        real*8 omega
cc
cc      use simple formula for the G_0 function
cc
        g0and=Xi*omega-Epsk(1)
        do i=1,ns-1
           g0and=g0and-tpar(i)**2/(Xi*omega-Epsk(1+i))
        end do

        g0and=dconjg(g0and)/(dconjg(g0and)*g0and)
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine datain(nss,uhub,xmu,hmag,ns,beta,wlow,wup,
     $     deltino,imaxmu,deltamu,iteramax, testdmft,dens_wanted,ifix,
     $     inew,iexp,th0,iauto)
c
      implicit none
c     hopping 't', Hubbard U
      real*8 uhub,xmu,hmag,deltamu
      real*8 wlow,wup,deltino,aaa
      real*8 beta,testdmft,dens_wanted,th0
      integer ns,imaxmu,iteramax,nss,is,ival,i,ifix,inew,iexp
      integer iauto
      character*80 text,text1
c
c     read(30,'(a)') text
c     read(30,'(a)') text1
      read(30,'(a)') text
c     hamiltonian parameters
      read(30,*) uhub,hmag
      read(30,'(a)') text
      read(30,*) beta, wlow, wup, deltino
c      read(45,*) varbeta
       
     

c*    
c*==> print options
c*
c*==> max. no. of lines of hamiltonian to print
c*    

      read(30,'(a)') text
      read(30,*) ns,imaxmu,deltamu,iteramax,testdmft   
      read(30,'(a)') text
      read(30,*)ifix,dens_wanted,inew,iauto
      read(30,'(a)') text
      read(30,*)th0,iexp
      write(6,69)ns
 69   format(1x,'# of sites      = ',i5,/)
      if (nss.lt.ns) then
         print*,'nss has to be at least  ',ns
         stop
      endif
      write(6,*)'number of conduction bath levels =',ns-1
      write(6,*)'number of iterations :',iteramax

      call rheader(15)
      read(15,*)
      
      do i=2,ns
         read(15,*)
      end do
      
      read(15,*)
      
      do i=1,ns-1
         read(15,*)
      end do
      
      read(15,*) xmu 


      return
      end
c*    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function isval(b,ns,abas,nss)
c*
      integer b(nss)
      integer *1 abas(nss)
c*
      isval = 0
      do i=1,ns
         isval = isval + (abas(i)+1) * b(i)
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function nchoos(n1,n2)
      xh = 1.0
      if(n2.lt.0) then
        nchoos = 0
        return
      endif
      if(n2.eq.0) then
        nchoos = 1
        return
      endif
      do 100 i = 1, n2
        xh = xh * float(n1+1-i)/float(i)
 100  continue
      nchoos = int(xh + 0.5)
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine statel(ival,b,nss,ns,anew)
      implicit none
      integer nss,ns,ival,i,is,nb
      integer *1 anew(nss)
      integer b(nss)
      is = ival
      do 1 i = ns,1,-1
        nb = is/b(i)
        anew(i) = nb - 1
        is = is - nb * b(i)
    1 continue
c     print*,'anew =',anew
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine computematrix(ioffset,nblock,ntot,b,ns,abas,
     $     nsp1,nmax,nleng,imat,idmat,iddomat,ntable,idouble)
c      program hamiltoniana
c     this is meant to be a subroutine that generates
c     1) the blocks
c     2) the non-diagonal matrix elements within each block (imat)
c     3) the matrix elements of d^{dagger} -> idmat
      implicit none

      integer nblock,i,nup,ndo,nchoos,j,nmax
      integer ns,nsp1,ntot,info
      integer ioffset(nblock+1)
      integer b(ns)
      integer ntable(ntot)
      integer*1 abas(ns)
      integer nleng(nblock)
      integer*1 imat(nmax,nmax,nblock)
      integer idmat(ntot)
      integer iddomat(ntot)
      integer*1 idouble(ntot)
      
      
      
c     compute ioffsets
      
      ioffset(1)=0
      ioffset(2)=1
      do i=2,nblock-1
         nup=mod(i-1,nsp1)
         ndo=(i-1)/nsp1
         ioffset(i+1)=ioffset(i)+nchoos(ns,nup)*nchoos(ns,ndo)
c         write(6,*) 'i', i, 'up', nup, 'down', ndo, 'ioffsets',ioffset(i)
      enddo
      ioffset(nblock+1)=ntot
      
      call buildblocks(ns,ntot,nblock,ntable,ioffset,nleng,
     $     b,abas,nsp1)
      
c     loop over blocks -> compute the Hamiltonian

      do i=1,ntot
         idouble(i)=0
      enddo

      do i=1,nblock
c         write(6,*) i, nleng(i)
         do j=1,nleng(i)
            call findnupndo(ntable(ioffset(i)+j)-1,nup,ndo,b,
     $           abas,ns)
            if (abas(1).eq.2) idouble(ioffset(i)+j)=1
         enddo
         call hamilt(ntable(ioffset(i)+1),nleng(i),
     $        imat(1,1,i),b,ns,abas,nmax)
         if (mod(i-1,nsp1).ne.ns) then
            call ddag(ntable(ioffset(i)+1),ntable(ioffset(i+1)+1),
     $           nleng(i),nleng(i+1),idmat(ioffset(i)+1),b,ns,
     $           abas,ntot,i)
         endif
         if(i.gt.(ns+1)) then
            call ddown(ntable(ioffset(i)+1),ntable(ioffset(i-ns-1)+1),
     $           nleng(i),nleng(i-ns-1),iddomat(ioffset(i)+1),b,ns,
     $           abas,ntot,i)
         endif
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hamilt(ntablein,nlen,imatin,b,ns,abas,nmax)
      implicit none
      integer nlen,ns,nmax,i,j,nup,ndo,ntnew,inew,k,isegno
      integer ntablein(nlen)
      integer*1 imatin(nmax,nmax)
      integer b(ns)
      integer*1 abas(ns)

      do i=1,nlen
         do j=1,nlen
            imatin(i,j)=0
         enddo
      enddo

  
c     diagonal part is already coded in ntab!!!!
c     non-diagonal part

      do i=1,nlen
         call findnupndo(ntablein(i)-1,nup,ndo,b,abas,ns)
c     if up is on site 1, it can hop
c         write(6,*) ntablein(i)
         if (abas(1).gt.0) then
c     loop over final sites
            do j=2,ns
               if (abas(j).le.0) then
                  ntnew=ntablein(i)+b(j)*(1-2*abas(j))+b(1)*
     $                 (1-2*abas(1))
                  call findinew(ntnew,ntablein,nlen,inew)
                  nup=0
                  do k=2,j-1
                     if (abas(k).gt.0) nup=nup+1
                  enddo
                  isegno= 1-2*mod(nup,2)
                  imatin(i,inew)=j*isegno
                  imatin(inew,i)=j*isegno
c                  write(6,*) 'inew', inew, 'imatin', imatin(i,inew)
               endif
            enddo
         endif
c     if down is on site 1, it can hop
         if (abas(1).eq.2.or.abas(1).eq.-1) then
            do j=2,ns
               if (abas(j).eq.0.or.abas(j).eq.1) then
                  ntnew=ntablein(i)+b(j)*(2*abas(j)-1)-b(1)*
     $                 abas(1)/abs(abas(1))
                  call findinew(ntnew,ntablein,nlen,inew)
                  
                  ndo=0
                  do k=2,j-1
                     if (abas(k).eq.2.or.abas(k).eq.-1) ndo=ndo+1
                  enddo
                  isegno=1-2*mod(ndo,2)
                  imatin(i,inew)=j*isegno
                  imatin(inew,i)=j*isegno                  
               endif
            enddo
         endif
        enddo


      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine findinew(ntnew,ntab,nlen,inew)
      implicit none
      integer nlen,i,ntnew,inew
      integer ntab(nlen)
      do i=1,nlen
         if (ntab(i).eq.ntnew) goto 777
      enddo

      write(6,*)'il numero ',ntnew,' manca'
      
 777  inew=i
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine buildblocks(ns,ntot,nblock,ntable,ioffset,nfound,
     $     b,abas,nsp1)
      implicit none

      integer ns,ntot,nblock
      integer ntable(ntot)
      integer ioffset(nblock+1)
      integer nfound(nblock)
      integer b(ns)
      integer*1 abas(ns)
      integer nup,ndo,nsp1,indblock,i

      do i=1,ntot
         call findnupndo(i-1,nup,ndo,b,abas,ns)
         indblock=nup+nsp1*ndo+1
         nfound(indblock)=nfound(indblock)+1

         if (nfound(indblock).gt.(ioffset(indblock+1)-
     $        ioffset(indblock)))
     $        then
            write(6,*)'Il blocco ',indblock,
     $           ' deve essere lungo', nfound(indblock)
            write(6,*)'nup, ndo =',nup,ndo
            write(6,*)'ci voglio mettere',i
            stop
         endif
         ntable(nfound(indblock)+ioffset(indblock))=i
c         write(6,*) ntable
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine findnupndo(nstato,nup,ndo,b,anew,ns)
      implicit none


      integer nb,i
      integer ns,nup,ndo,nstato,is
      integer *1 anew(ns)
      integer b(ns)

      is = nstato
      nup=0
      ndo=0
      do 1 i = ns,1,-1
        nb = is/b(i)
        anew(i) = nb - 1
        is = is - nb * b(i)
        if (anew(i).gt.0) nup=nup+1
        if (anew(i).eq.2.or.anew(i).eq.-1) ndo=ndo+1
    1 continue
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ddag(ntab0,ntab1,
     $     nlen0,nlen1,idmatin,b,ns,abas,ntot,iblocco)
      implicit none

      integer nlen0,nlen1,ns,j,iold,nup,ndo,inew,k
      integer ntab0(nlen0)
      integer ntab1(nlen1)
      integer b(ns)
      integer*1 abas(ns)
      integer ntot,iblocco
      integer idmatin(nlen0)
      do j=1,nlen0
         iold=ntab0(j)
         call findnupndo(iold-1,nup,ndo,b,abas,ns)
         if (abas(1).le.0) then
            inew=iold+(1-2*abas(1))
            do k=1,nlen1
               if (ntab1(k).eq.inew) then
                  idmatin(j)=k
                  goto 555
               endif
            enddo
            write(6,*)'Non trovo lo stato ottenuto applicando d+'
 555        continue
         else
            idmatin(j)=0
         endif
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ddown(ntab0,ntab1,
     $     nlen0,nlen1,iddomatin,b,ns,abas,ntot,iblocco)
      implicit none

      integer nlen0,nlen1,ns,j,iold,nup,ndo,inew,k,isegno
      integer ntab0(nlen0)
      integer ntab1(nlen1)
      integer b(ns)
      integer*1 abas(ns)
      integer ntot,iblocco
      integer iddomatin(nlen0)
      do j=1,nlen0
         iold=ntab0(j)
         call findnupndo(iold-1,nup,ndo,b,abas,ns)


cccc     If down is on site 1

 
         if (abas(1).eq.-1.or.abas(1).eq.2) then
            inew=iold + (1-2*abas(1))/3

ccccccccccc          check delle operazioni   ccccccccccccccccc
 
c            write(6,*) 'abas', abas
c            write(6,*)'in', inew,'io', iold,'dsum',(1-2*abas(1))/3
c            call findnupndo(inew-1,nup,ndo,b,abas,ns)
c            write(6,*) 'abasnew', abas
c            call findnupndo(inew-1,nup,ndo,b,abas,ns)
c            write(6,*) 'inew', inew, 'abas', abas
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  

cc          ddown operator has to hop nup spin to reach his spin down
          
            isegno=1-2*mod(nup,2)
       
            do k=1,nlen1
               if (ntab1(k).eq.inew) then
                  iddomatin(j)=k*isegno
                  goto 555
               endif
            enddo
            

            write(6,*)'Non trovo lo stato ottenuto applicando d'
 555        continue
         else
            iddomatin(j)=0
         endif
      enddo
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine diag(tpar,epsk,uhub,hmag,ioffset,nleng,idmat,
     $     imat,nblock,ns,nsp1,nmax,ntable,work,eig,zeig,eigall,
     $     abas,
     $     ntot,realmat,beta,Xi,piob,range,wup,wlow,Iwmax,Gw,
     $     om,eigold,zold,b,rlehm,xmat,nmaxx,zpart,summy,xmu,
     $     emin,idouble,doublep,double)

      implicit none

      real*8 diffe,piccolo
      integer ns,ntot,nblock,nmax,Iwmax,i,j,k,nmatrix
      integer neig,meig,iaut,info,nsp1,icount,iomega,nmaxx
      real*8 tpar(ns),epsk(ns),uhub,hmag,sum,rmel,emin,beta,pi
      real*8 wup,wlow,range,piob
c      real*8 time,mclock
      integer ioffset(nblock+1),nleng(nblock)
      integer idmat(ntot)
      integer*1 idouble(ntot)
      integer*1 imat(nmax,nmax,nblock)
      integer ntable(ntot)
      integer*1 abas(ns)
      integer b(ns)
      integer nup,ndo,is
      real*8 work(3*nmaxx)
      real*8 eig(nmaxx),zeig(nmaxx,nmaxx)
      real*8 eigold(nmax)
      real*8 zold(nmax,nmax)
      real*8 realmat(nmaxx*(nmaxx+1)/2)
      real*8 xmat(nmaxx,nmaxx)
      complex*16 Xi
      complex*16 Gw(0:Iwmax)
      real*8 om(0:Iwmax)
      real*8 dens,xmu           
      real*8 rlehm(nmaxx,nmaxx,nblock)
      real*8 eigall(ntot)
      real*8 doublep(ntot)
      real*8 double
      real*8 zpart,Fener, logtest, Fnormal, znormal, kinsim
      real*8 summy
      
     
      piccolo=0.000001d0
      
      pi=dacos(-1.d0)
c     write(6,*) 'Pi=',  pi
      
      epsk(1)=-xmu

c      time=mclock()/100.d0

      do i=0,Iwmax
         Gw(i)=dcmplx(0.d0,0.d0)
      enddo
      
      summy=0.d0

      do i=1,nmaxx
         do j=1,nmaxx
            do k=1,nblock
               rlehm(i,j,k)=0.d0
            enddo
         enddo
      enddo
      
      do i=1,nblock
         rewind(96)
         write(96,*)'blocco #',i
         call flush(96)
         do j=1,nmaxx
            do k=j+1,nmaxx
               xmat(j,k)=0.d0
               xmat(k,j)=0.d0
            enddo
            xmat(j,j)=100000.d0
         enddo


c     do j=1,nmax*(nmax+1)/2
c     realmat(j)=0.d0
c     enddo
c     build up the actual matrix
c     do j=1,nleng(i)
c     diagonal matrix elements
c     call findnupndo(ntable(ioffset(i)+j)-1,nup,ndo,
c     $           b,abas,ns)
c     do k=1,ns
c     realmat(j+(j-1)*j/2)=
c     $              realmat(j+(j-1)*j/2)+epsk(k)*dabs(dfloat(abas(k)))
c     enddo
c            if (abas(1).eq.2) then 
c     realmat(j+(j-1)*j/2) = realmat(j+(j-1)*j/2)+uhub
c     else
c     realmat(j+(j-1)*j/2) = realmat(j+(j-1)*j/2)+
c     $              hmag*dfloat(abas(k))
c     endif
c     off-diagonal m.el.
c     do k=j,nleng(i)
c     if (imat(j,k,i).ne.0) then
c     realmat(j+(k-1)*k/2)=tpar(abs(imat(j,k,i))-1)*
c     $                 imat(j,k,i)/abs(imat(j,k,i))
c     endif
c     enddo            
c     enddo
         
         do j=1,nleng(i)
            call findnupndo(ntable(ioffset(i)+j)-1,nup,ndo,
     $           b,abas,ns)
            xmat(j,j)=0.d0
            do is=1,ns
               xmat(j,j)=xmat(j,j)+epsk(is)*dabs(dfloat(abas(is)))
            enddo
            if (abas(1).eq.2) then
               xmat(j,j)=xmat(j,j)+uhub
            else
               xmat(j,j)=xmat(j,j)-hmag*dfloat(abas(1))
            endif
            do k=j+1,nleng(i)
               if (imat(j,k,i).ne.0) then
                  xmat(j,k)=xmat(j,k)+tpar(abs(imat(j,k,i))-1)*
     $                 imat(j,k,i)/abs(imat(j,k,i))
                  xmat(k,j)=xmat(k,j)+tpar(abs(imat(j,k,i))-1)*
     $                 imat(j,k,i)/abs(imat(j,k,i))
               endif
            enddo
         enddo

         icount=0

         do j=1,nmaxx
            do k=j,nmaxx
               icount=icount+1
               realmat(icount)=xmat(k,j)
            enddo
         enddo
            
       
         call dspev('V','L',nmaxx,realmat,eig,zeig,
     $        nmaxx,work,info)


         
         if (info.ne.0) write(6,*)'INFO =',info
c     per capire che succede, applichiamo H agli autovettori....

         
         
c     if nup.ne.0 compute \sum_nm |<n|d+|m>|^2/(iomega - (E_n - E_m).....
         if (mod(i-1,nsp1).ne.0) then
c     loop over |ndo,nup-1> basis states
            do neig=1,nleng(i)
               do meig=1,nleng(i-1)
                  rmel=0.d0
                  do j=1,nleng(i-1)
                     if (idmat(ioffset(i-1)+j).ne.0) then
                        rmel=rmel+zold(j,meig)*
     $                       zeig(idmat(ioffset(i-1)+j),neig)
c     if (neig.eq.1.and.meig.eq.1) 
c     $                  write(99+i,*)'d+ :',idmat(ioffset(i-1)+j),'<-',j
                     endif
                  enddo
                  rlehm(neig,meig,i)=rmel**2
               enddo
               
            enddo
         endif

         do neig=1,nleng(i)
            doublep(neig+ioffset(i))=0.d0
            do j=1,nleng(i)
               doublep(neig+ioffset(i))=doublep(ioffset(i)+neig)+
     $              zeig(j,neig)*zeig(j,neig)*
     $              dfloat(idouble(ioffset(i)+j))
            enddo
         enddo


c     then store new vectors in old vectors and go on
         do j=1,nleng(i)
            eigold(j)=eig(j)
            do k=1,nleng(i)
               zold(j,k)=zeig(j,k)
            enddo
            eigall(ioffset(i)+j)=eig(j)
         enddo
      enddo

      call flush(6)
      emin=1.d16
      do i=1,ntot
         if (eigall(i).lt.emin) emin=eigall(i)
      enddo
      double=0.d0
      zpart=0.d0
      do i=1,ntot
         eigall(i)=eigall(i)-emin
         zpart=zpart+dexp(-beta*eigall(i))
         double=double+doublep(i)*dexp(-beta*eigall(i))
      enddo
      
      double=double/zpart
      

c      time=mclock()/100.d0-time
c      write(6,*)'Tempo per diagonalizzare: ',time


      write(6,*)'Z =',zpart
      write(6,*)'emin =',emin
      call flush(6)
c      time=mclock()/100.d0
    
      do i=1,nblock
         rewind(77)
         write(77,*)'blocco ',i
         call flush(77)
         if (mod(i-1,nsp1).ne.0) then
            do iomega=0,Iwmax
               do neig=1,nleng(i)
                  do meig=1,nleng(i-1)
                     Gw(iomega)=Gw(iomega)+rlehm(neig,meig,i)
     $                    /(Xi*om(iomega)-
     $                    (eigall(ioffset(i)+neig)-
     $                    eigall(ioffset(i-1)+meig)))*
     $                    (dexp(-beta*eigall(ioffset(i)+neig)) +
     $                    dexp(-beta*eigall(ioffset(i-1)+meig)))
                  enddo                  
               enddo
            enddo
         endif
         
c      mi calcolo la densimp...
         if(i.gt.1)then
           do neig=1,nleng(i)
             do meig=1,nleng(i-1)
               summy=summy+rlehm(neig,meig,i)*
     $              dexp(-beta*eigall(ioffset(i)+neig))
             enddo
           enddo
         endif
      enddo
c      do i=0, Iwmax
c         write(99,'(3f17.10)') om(i),dreal(Gw(i)),dimag(Gw(i))
c      enddo

c      ...correttamente normalizzata
      summy=2.d0*summy/zpart
c  normalizzo anche Gw
      do i=0,Iwmax
         Gw(i)=Gw(i)/zpart
      enddo
c      time=mclock()/100.d0-time
c      write(6,*)'Tempo per le G:',time



      kinsim=0.d0
      do i=0, Iwmax
         kinsim=kinsim + Gw(i)*Gw(i)+
     $        dconjg(Gw(i))*dconjg(Gw(i))
      enddo
      kinsim=0.5d0*kinsim/beta
      write(6,*) 'Ekin=', kinsim

      kinsim=kinsim-0.5d0/(pi*(om(Iwmax)+pi/beta))
      write(6,*) 'Ekin_TOT=', kinsim
      write(6,*) 'Im[G_iw]  =',  dimag(Gw(Iwmax)),
     $     'vs.  1/om =', 1.d0/om(Iwmax)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     c                      Calcolo dell'energia libera             
      
      znormal=0.d0
         
      Fnormal=0.d0     
      
      
      do i=2, ns 
         if(beta*epsk(i).lt.-20.d0) then
            Fnormal =Fnormal+2.d0*epsk(i)
         else 
            znormal=(1.d0+dexp(-beta*epsk(i)))**2
c            write(6,*) 'Znormal=', znormal
            Fnormal=Fnormal-dlog(znormal)/beta 
         endif 
c     write(6,*) 'fnormal=', Fnormal     
      enddo
      
      
      Fener= emin -1.d0/beta*(dlog(zpart))-0.5*kinsim-Fnormal  
c      write(6,*) 'Free Energy=', Fener  
      
         
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine computegreal(Gwreal,Iwmaxreal,rlehm,omr,eigall,beta,
     $     ioffset,nleng,nblock,ntot,nsp1,pi,nmax,zpart,sum)

      integer Iwmaxreal,iomega,neig,meig,nblock,ntot,nsp1,i
      integer ioffset(nblock+1),nleng(nblock)
      complex*16 Gwreal(0:Iwmaxreal)
      real*8 eigall(ntot)
      real*8 rlehm(nmax,nmax,nblock)
      complex*16 omr(0:Iwmaxreal)
      real*8 beta,pi,zpart
      real*8 sum,x

      do i=1,nblock
         if (mod(i-1,nsp1).ne.0) then
            do iomega=0,Iwmaxreal
               do neig=1,nleng(i)
                  do meig=1,nleng(i-1)
                     Gwreal(iomega)=Gwreal(iomega)+
     $                    rlehm(neig,meig,i)
     $                    /(omr(iomega)-
     $                    (eigall(ioffset(i)+neig)-
     $                    eigall(ioffset(i-1)+meig)))*
     $                    (dexp(-beta*eigall(ioffset(i)+neig)) +
     $                    dexp(-beta*eigall(ioffset(i-1)+meig)))
                  enddo                  
               enddo
            enddo
         endif
      enddo

      sum=0.d0

      do i=0,Iwmaxreal
         Gwreal(i)=-Gwreal(i)/zpart/pi
         x=dmax1(beta*dreal(omr(i)),-40.d0)
         x=dmin1(x,40.d0)
         sum=sum+dimag(Gwreal(i))/(1.d0+dexp(x))
      enddo
      sum=sum*(dreal(omr(Iwmaxreal))-dreal(omr(0)))/
     $     dfloat(Iwmaxreal)
      sum=2.d0*sum
      write(6,*)'densold =  ',sum
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine dos2(sigre,Iwmaxreal,omr,Pi,xmu,beta)
      
      integer i,jen,Iwmaxreal,Emax
      real*8 Pi,E,dos(0:Iwmaxreal),sum,xmu,tpri,beta
      complex*16 Sigre(0:Iwmaxreal), omr(0:Iwmaxreal)
      external DenOS
      parameter(Emax=50000)
      include 'tpri.dat'
 
      write(6,*) 'Computing the interacting DOS'
      do i=0,Iwmaxreal
         dos(i)=0.d0
         do jen=0,Emax-1
            E=-1.d0+2.d0*dfloat(jen)/dfloat(Emax)+1.d0/dfloat(Emax)
            dos(i)=dos(i)-(1.d0/Pi)*dimag(1.d0/(omr(i)-E+xmu-tpri
     $           -sigre(i)))
     $           *DenOS(E)*2.d0/dfloat(Emax)
        
         enddo
      enddo
      
      sum=0.d0
      do i=0,Iwmaxreal
         sum=sum+dos(i)*dreal(omr(1)-omr(0))/(dexp(beta*dreal(omr(i)))
     $        +1.d0)*2.d0
c           write(6,*) i, dos(i)
         write(57,'(3f17.10)') dreal(omr(i)), dos(i),sum
      enddo


      write(6,*) 'Final Spectral weight check=', sum

      sum=0.d0
      do jen=0,Emax-1
         E=-1.d0+2.d0*dfloat(jen)/dfloat(Emax)+1.d0/dfloat(Emax)
         sum=sum+DenOS(E)*2.d0/dfloat(Emax)
         write(58,'(3f17.10)') E+tpri,DenOS(E), sum
      enddo

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine selfconst2(g0wand,Gw,Gloc,Xi,Pi,Iwmax,om,xmu,beta)
      
      integer i,j,k,Iwmax,jen,Emax
      real*8 Pi,E,xmu,sum,tpri,beta
      complex*16 Xi
      
      real*8 om(0:Iwmax)
      complex*16 g0wand(0:Iwmax),Gw(0:Iwmax),Gloc(0:Iwmax),gand(0:Iwmax)
      complex*16 W(0:Iwmax)
      external DenOS
      parameter(Emax=100)

      include 'tpri.dat'
      write(6,*) 'Tpri=', tpri

      do i=0,Iwmax
c         g0wand(i)=1.d0/g0wand(i)
         gand(i)=1.d0/g0wand(i)
         Gw(i)=1.d0/Gw(i)
         W(i)= Xi*om(i)-gand(i)+Gw(i)+xmu-tpri
      enddo 

      
      rewind(61)
      sum=0.d0
      do i=0,Iwmax
         Gloc(i)=(0.d0,0.d0)
         do jen=1,Emax-1
            E=-1.d0+2.d0*dfloat(jen)/dfloat(Emax)+1.d0/dfloat(Emax)
            Gloc(i)=Gloc(i)+2.d0*DenOS(E)/(W(i)-E)*
     $           1.d0/dfloat(Emax)
c            if (jen.ge.3) then 
c             write(6,*)'done'
c            stop
c             endif
         enddo
         write(61,'(5f17.10)') om(i),dreal(Gloc(i)),dimag(Gloc(i)),
     $        dreal(gand(i)-Gw(i)),dimag(gand(i)-Gw(i))
         sum=sum+2.d0*(Gloc(i)+dconjg(Gloc(i)))/beta 
         Gloc(i)=1.d0/Gloc(i)+gand(i)-Gw(i)
         Gloc(i)=0.5d0/Gloc(i)+0.5d0*g0wand(i)
c         g0wand(i)=1.d0/g0wand(i)
      enddo
      call flush(61)

      write(6,*) 'Check interacting dens', sum+1.d0

      sum=0.d0
      do jen=0,Emax-1
         E=-1.d0+2.d0*dfloat(jen)/dfloat(Emax)+1.d0/dfloat(Emax)
         sum=sum+2.d0*DenOS(E)/dfloat(Emax)
      enddo
      write(6,*) 'Check Dens = ', sum

      return 
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function DenOS(E)
      
     
      real*8 E,Pi
      real*8 DenOS
      real*8 tpri,argum
      
      external EllypticK
      include 'tpri.dat'
c      write(6,*) 'Tpri=', tpri 

      Pi=dacos(-1.d0)
      argum=-(-1.d0+E**2)/
     $           (1.d0+4.d0*E*tpri+4.d0*tpri**2)
      DenOS=2.d0/(Pi**2*dsqrt(1.d0+4.d0*E*tpri+4*tpri**2))*
     $           EllypticK(argum)

      
      return
      end
      



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function EllypticK(x)
      real*8 EllypticK,x,a,b,c,s,t
      integer check
      data e1 /1e-8/
      a=1
      i=1
      b=sqrt(1-x)
      t=x
      s=0
      check=0
      open(44,file='checkfile',form='formatted',status='unknown')
      open(53,file='checkab',form='formatted',status='unknown')
      write(6,*)e1
      call flush(6)
  1   s=s+t
      c=(a-b)/2
      t=(a+b)/2
      i=2*i
      b=sqrt(a*b)
      a=t
      t=i*c**2
      check=check+1
      write(44,*)check
c      write(6,*)'test'
c      call flush(6)
      write(53,*)a,b,t,s,c
c      if (check.gt.) then
c         stop
c      endif
      if ((abs(c).gt.e1*a).or.(t.gt.e1*s)) goto 1
            if ((a+b).eq.0) stop 99
      EllypticK=3.1415926535897/(a+b)
      write(45,*)'done'
      call flush(45)
      close(44)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine selfconst3(g0wand,Gw,Gloc,Xi,Pi,Iwmax,om,xmu)
    

      integer i,j,k,Iwmax,ix,iy,iz
      real*8 Pi,xmu,Ek,var
      complex*16 Xi
      real*8 om(0:Iwmax)
      complex*16 g0wand(0:Iwmax),Gw(0:Iwmax),Gloc(0:Iwmax)
      complex*16 W(0:Iwmax)
      real*8 Lx,Ly,Lz,Vol,kx,ky,kz,k0x,k0y,k0z,tsc,tfcc,d
      parameter (Lx=100)
      parameter (var=0.5d0)
      parameter(d=4.d0*var)
      

      do i=0,Iwmax
         g0wand(i)=1.d0/g0wand(i)
         Gw(i)=1.d0/Gw(i)
         W(i)= Xi*om(i)+xmu- g0wand(i)+Gw(i)
      enddo 

      
 
      Ly=Lx
      Lz=Lx
      Vol = Lx*Ly*Lz
      k0x = Pi*2.d0/Lx
      k0y = Pi*2.d0/Ly
      k0z = Pi*2.d0/Lz

! fcc/cubic Gitter
! Hopping fcc
c      tfcc = 4.d0/sqrt(6.d0*d*d+12.d0)
! Hopping cubic
c      tsc  = d*tfcc/2.d0
      tsc=var*dsqrt(2.d0)/dsqrt(3.d0)
      write(6,*) 'Variance = ', var, 'D = ',tsc*3.d0
      tfcc=0.d0
! loop over momenta separated into four blocks
! uses symmetry kx -- ky -- kz
! ----------------------------------------------------

      kx = 0.d0
      do ix=1,Lx/2-1
c         if (ix/10*10.eq.ix) write(FOUT,*) 'ix=',ix
         kx = kx+k0x
         ky = 0.d0
         do iy=1,Ly/2-1
            ky = ky+k0y
            kz = 0.d0
            do iz=1,Lz/2-1
               kz = kz+k0z
c     here comes the relevant dispersion relation
               Ek = tfcc*((cos(kx)+cos(ky))*cos(kz)+cos(kx)*cos(ky))
c               if (D.le.1.d0) then
c                  Ek = Ek+tsc*(cos(2.d0*kx)+cos(2.d0*ky)+cos(2.d0*kz))
c               else  
                  Ek = Ek+tsc*(cos(kx)+cos(ky)+cos(kz))
c               endif
c               rzero=rzero+8.d0
c               rone=rone+8.d0*eq0
c               rtwo=rtwo+8.d0*eq0*eq0
                 
               
c     here comes the loop over Matsubara frequencies
               do i=0,Iwmax
                   Gloc(i)=Gloc(i)+(8.d0,0.d0)/(W(i)-Ek)
                enddo
            enddo
            
         enddo
      enddo
      
! ----------------------------------------------------

      kx = -pi
      do ix=1,2
         kx = kx+pi
         ky = 0.d0
         do iy=1,Ly/2-1
            ky = ky+k0y
            kz = 0.d0
            do iz=1,Lz/2-1
               kz = kz+k0z
                                ! here comes the relevant dispersion relation
               Ek = tfcc*((cos(kx)+cos(ky))*cos(kz)+cos(kx)*cos(ky))
               if (D.le.1.d0) then
                  Ek = Ek+tsc*(cos(2.d0*kx)+cos(2.d0*ky)+cos(2.d0*kz))
               else  
                  Ek = Ek+tsc*(cos(kx)+cos(ky)+cos(kz))
               endif
c               rzero=rzero+12.d0
c               rone=rone+12.d0*eq0
c               rtwo=rtwo+12.d0*eq0*eq0
               
! here comes the loop over Matsubara frequencies
 
               do i=0,Iwmax
                  Gloc(i) =Gloc(i)+(12.d0,0.d0)/(W(i)-Ek)
               enddo
               
            enddo
         enddo
      enddo
! ----------------------------------------------------

      kx = -pi
      do ix=1,2
         kx = kx+pi
         ky = -pi
         do iy=1,2
            ky = ky+pi
            kz = 0.d0
            do iz=1,Lz/2-1
               kz = kz+k0z
!     here comes the relevant dispersion relation
               Ek = tfcc*((cos(kx)+cos(ky))*cos(kz)+cos(kx)*cos(ky))
               if (D.le.1.d0) then
                  Ek = Ek+tsc*(cos(2.d0*kx)+cos(2.d0*ky)+cos(2.d0*kz))
               else  
                  Ek = Ek+tsc*(cos(kx)+cos(ky)+cos(kz))
               endif
c               rzero=rzero+6.d0
c               rone=rone+6.d0*eq0
c               rtwo=rtwo+6.d0*eq0*eq0
               
! here comes the loop over Matsubara frequencies

               do i=0,Iwmax
                  Gloc(i) = Gloc(i)+ (6.d0,0.d0)/(W(i)-Ek)
               enddo
               
            enddo
         enddo
      enddo
c     ----------------------------------------------------

      kx = -pi
      do ix=1,2
         kx = kx+pi
         ky = -pi
         do iy=1,2
            ky = ky+pi
            kz = -pi
            do iz=1,2
               kz = kz+pi
                                ! here comes the relevant dispersion relation
               Ek = tfcc*((cos(kx)+cos(ky))*cos(kz)+cos(kx)*cos(ky))
               if (D.le.1.d0) then
                  Ek=Ek+tsc*(cos(2.d0*kx)+cos(2.d0*ky)+cos(2.d0*kz))
               else  
                  Ek = Ek+tsc*(cos(kx)+cos(ky)+cos(kz))
               endif
c     rzero=rzero+1.d0
c     rone=rone+eq0
c     rtwo=rtwo+eq0*eq0
                 
c     here comes the loop over Matsubara frequencies
      
               do i=0,Iwmax
                  Gloc(i)=Gloc(i)+(1.d0,0.d0)/(W(i)-Ek)
               enddo
               
            enddo
         enddo
      enddo
! ----------------------------------------------------

         


      do i=0,Iwmax
         Gloc(i) = Gloc(i)/Vol
         Gloc(i)=1.d0/Gloc(i)+g0wand(i)-Gw(i)
         Gloc(i)=1.d0/Gloc(i)
      enddo

c         write(FOUT,*)
c         write(FOUT,'(a20,3f16.8)') 'first three moments:'&
c             ,real(rzero/Vol),real(rone/Vol),real(rtwo/Vol)
c         write(FOUT,*)



      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      function EllipK(x1)
c      real*8 x1,a0,a1,a2,a3,a4,b0,b1,b2,b3,b4,EllipK
c      if (x1.le.0.or.x1.ge.1) then
c            print *,x1
c            stop 99
c      endif
c      data a0,a1,a2,a3,a4 /1.38629436112,0.09666344259,0.03590092383,
c     .      0.03742563713,0.01451196212/
c      data b0,b1,b2,b3,b4 /0.5,0.12498593597,0.06880248576,
c     .      0.03328355346,0.00441787012/
c      EllipK=(((a4*x1+a3)*x1+a2)*x1+a1)*x1+a0-
c     .      ((((b4*x1+b3)*x1+b2)*x1+b1)*x1+b0)*dlog(x1)
c      return
c      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



