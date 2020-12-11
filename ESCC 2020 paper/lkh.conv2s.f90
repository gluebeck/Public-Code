      subroutine func(u,npar,fout)

      implicit real*8(a-h,o-z)
      save
      parameter(stemcells=1.d6)
      character sex
      real*8,dimension(42,86) :: pop,yinc,haz  
      real*8 u(npar),age(86),age1(86),haz_conv_2s(86),tum_conv_2s(86)
      real*8,dimension(64) :: xg1,wg1,xg2,wg2
      real*8,dimension(64) :: f1,f2,f,dfdt,tmp,dum1,dum2,s,uu,sk,haz_2s,sur_2s,hyear,nye,nye0,sfunc,SurvNu
      real*8 alphaI,alphaC,gI,gC,nu,mu0,mu1,mu2,p2,q2,rh0,rootarg

      COMMON /IFLAGS/IFLAG,IMCMC,IBOOT,ISEED,IERR  
      common /params/alphaI,alphaC,gI,gC,mu1,mu2,p2,q2

        if(iflag == 1) then

           nbio=5
           noyrs = 42
           ngs1 = 32
           ngs2 = 16
           
           do i=1,86
              age(i)=dfloat(i)-.5d0
           enddo

           age(86)=87.5d0
 
           open(unit=20,file='./SEER9_2016/SEER9-2016-WM-POP.txt',status='unknown') 
           open(unit=21,file='./SEER9_2016/SEER9-2016-WM-INC.txt',status='unknown') 

           do i=1,noyrs
                 read(20,*) (pop(i,j),j=1,86)
                 read(21,*) (yinc(i,j),j=1,86)
           enddo
           close(20); close(21)

           call legaus(0.d0,1.d0,ngs1,xg1,wg1)
           call legaus(0.d0,1.d0,ngs2,xg2,wg2)
           iflag=2

        endif

!       model parameters

        alpha=17.4d0
        tlag = u(5)

!        mu0 = dsqrt(u(1)*alphaI/(u(3)*stemcells)); mu1=mu0=r0; mu2=1.d-7*u(2)
        nu = u(1); mu1=u(2); mu2=1.d-7*u(3)
        tc = -dlog(1.d-15/nu)/nu ! for scaling highest Gauss pt

        g0 = u(4)
              
        ! year/period param

        g1 = u(nbio+1)
        g2 = u(nbio+2)

        b1 = u(nbio+3)
        b2 = u(nbio+4)
        w1 = u(nbio+5)
        w2 = u(nbio+6)

        ! range/reference year param
        hy0 = u(nbio+7)
        by0 = hy0
        gy0 = u(nbio+8)
        scale = u(nbio+9)
    
! ----  do likelihood
        fout=0.d0

        do i=20,86
           
           t=max(age(i)-tlag,.1d0)
! ----     limit for s-integration
           tlim = min(tc,t)
           s(1:ngs1)=tlim*xg1(1:ngs1) ! integrating [0,t]  !set up convolution
           uu(1:ngs1)=t-s(1:ngs1)

!     ----  compute impact of historic time
           do j=1,noyrs
              
              cyear = j+1974.d0   ! SEER 1975-2011
              byear = cyear-age(i)
              
              if(byear < 1961.d0) then
                 
              ! hy0: reference year on historical timeline

              ! forcing log nu/nu0 to zero at time s0 and year h0:  w2 = 1.d0/(hy0-byear-s0)

              ! cohort dependent cell proliferation
              dum = byear-gy0
              gfac = dexp(g1*dum*(1.d0+g2*dum))
              alphaI = alpha*gfac
              gI = g0*gfac

              ! cohort dependent mu1
              ! mu1 = r0*gfac
              
              ! for 2-stage haz
              rootarg=gI*gI+4.d0*alphaI*mu2*gfac
              p1 = 0.5d0*(-gI-dsqrt(rootarg))
              q1 = 0.5d0*(-gI+dsqrt(rootarg))
           
              f1(1:ngs1) = q1 * dexp(-p1 * uu(1:ngs1))
              f2(1:ngs1) = p1 * dexp(-q1 * uu(1:ngs1))
              f(1:ngs1)  =  f1(1:ngs1) - f2(1:ngs1)
              dfdt(1:ngs1) = p1 * q1 * (dexp(-q1 * (uu(1:ngs1))) - dexp(-p1 * (uu(1:ngs1))))

              ! birth cohort effect on nu. byear only ref by0 
              dby = byear-by0
              if(dby > 0.d0) dby = 0.d0
              bye = b1*dby*dby*(1.d0+b2*dby)
              
              ! secular/historical trend
              hyear(1:ngs1) = byear+s(1:ngs1)
              dum1(1:ngs1) = hyear(1:ngs1)-hy0
              do k=1,ngs1
                 if(dum1(k) > 0.d0) dum1(k)=0.d0
              enddo
              
              sfunc(1:ngs1) = dlog(s(1:ngs1)/scale) 
              nye(1:ngs1) = dexp(w1*sfunc(1:ngs1)*dum1(1:ngs1)*dum1(1:ngs1)*(1.d0+w2*dum1(1:ngs1))+bye+dlog(nu)) 
                             
!!!!!!  ----normalization factor for the conversion density with a (historical) time dependent nu

              do k=1,ngs1
                 sk(1:ngs2)=s(k)*xg2(1:ngs2)
                 hyear(1:ngs2) = byear+sk(1:ngs2)
                 dum1(1:ngs2) = hyear(1:ngs2)-hy0
                                  
                 do l=1,ngs2
                    if(dum1(l) > 0.d0) dum1(l)=0.d0
                 enddo
              
                 sfunc(1:ngs2) = dlog(sk(1:ngs2)/scale) 
                 nye0(1:ngs2) =  dexp(w1*sfunc(1:ngs2)*dum1(1:ngs2)*dum1(1:ngs2)*(1.d0+w2*dum1(1:ngs2))+bye+dlog(nu)) 
                 SurvNu(k) = dexp(-s(k)*sum(wg2(1:ngs2)*nye0(1:ngs2)))
              enddo
!!!!!!

              haz_2s(1:ngs1) = mu1*dfdt(1:ngs1)/f(1:ngs1)/alphaI
              sur_2s(1:ngs1) = ((q1-p1)/f(1:ngs1))**(mu1/alphaI)  

              ! haz_2s(1:ngs) = hye(1:ngs)*mu1*dfdt(1:ngs)/f(1:ngs)/alphaI
              ! sur_2s(1:ngs) = dexp(hye(1:ngs)*mu1*dlog((q1-p1)/f(1:ngs))/alphaI)

!     ----  convolving 2-stage haz

              tum_conv_2s(i) = tlim*sum(nye(1:ngs1)*wg1(1:ngs1)*SurvNu(1:ngs1)*(1.d0-sur_2s(1:ngs1)))
              haz_conv_2s(i) = tlim*sum(nye(1:ngs1)*wg1(1:ngs1)*SurvNu(1:ngs1)*haz_2s(1:ngs1)*sur_2s(1:ngs1))/ &
              (1.d0-tum_conv_2s(i))      
           
              haz(j,i)   =  haz_conv_2s(i)

              xlam =  haz(j,i) * pop(j,i) ! Poisson mean  
              xobs =  yinc(j,i) ! observation
              fout = fout + xlam -xobs*dlog(xlam)

              endif
        
           enddo
        enddo

        if(iflag == 3) then
        print*,'post-processing ...'

        do i=30,86
           do j = 1,noyrs
              write(10,*) age(i),1974+j,yinc(j,i),pop(j,i),haz(j,i),yinc(j,i)/pop(j,i)
           enddo
        enddo
        endif
        return
        end
        





