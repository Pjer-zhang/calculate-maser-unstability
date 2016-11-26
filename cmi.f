c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+      Program CMI_Ring_Beam_IAW (Intrinsic Alfven Waves)                     +
c+      by CB Wang  02/22/2016           +
c+                                                                            +
c+      Code for calculating the growth rate of maser instability in           +
c+      cylindrical coordinate using the general expression of the growth rate.+
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c      implicit double precision (a-h,o-z)
c      parameter(pi=3.1415926535898d+00,eps=1.d-8)
c      common /para1/f,en2,isign
c      common /para2/alfa,beta,up0,uz0,fc,anorm
c      fce=90.d0
c      fpe=1.0d0*fce
c      fpe2=fpe*fpe
c      fce2=fce*fce
c      pitch=20.d0
c      u0=0.3
c      up0=u0*dsin(pitch*pi/180.d0)
c      uz0=u0*dcos(pitch*pi/180.d0)
c      alfa=0.1d0*up0
c      beta=0.2d0*uz0
c      isign=0
c      open(10,file='gr.dat')
c      do theta=1.d0,180.d0,1.d0
c      do fratio=0.6d0,2.5d0,1.d-2
c      if (isign .eq. 0) then
c      f_cut=fpe
c      elseif (isign .eq. 1) then
c      f_cut=(dsqrt(4.d0*fpe2+fce2)+fce)/2.d0
c      endif
c      f=fce*fratio
c      if (f .gt. f_cut) then
c      call growth(theta,fpe,fce,gr)
c      else
c      gr=0.d0
c      endif
c      write(10,50) fpe/fce,theta,fratio,gr
c      enddo
c      enddo
c50      format(1x,4(g12.3))
c      close(10)
c      end


c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+    Programe For calculating the maxium growth rate at different values of   +
c+    fpe/fce. The maxium growth rate, its corresponding wave frequency and    +
c+    propagation angle are calculated.                                        +
c+                                                                             +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        implicit double precision (a-h,o-z)
        parameter(pi=3.1415926535898d+00,eps=1.d-28)
        parameter (u0=0.3d0,imode=1)
        character*20 filename
        common /para1/f,en2,isign
        common /para2/alfa,beta,up0,uz0,fc,anorm
        common /para_pitch/pitch
        isign=0
        fce=1.d0
        fce2=fce*fce
        do pitch=90.D0, 90.D0, 1.D0
          up0=u0*dsin(pitch*pi/180.d0)
          uz0=u0*dcos(pitch*pi/180.d0)
            amiu0=dcos(pitch*pi/180.d0)
          alfa=0.1d0*u0
          beta=0.1D0
            write(*,*) 'u0, amiu0, alfa, beta=',u0, amiu0, alfa, beta
            anorm=Anormal2(up0,uz0,alfa,beta) !Constant for the normalization of distribution function.
            write(*,*) 'anorm=',anorm
            if (isign.eq.0) then
               WRITE(FILENAME,'(''gr_max_o'',I1.1,''_q'',I2.2,
     &             ''_1.dat'')') imode,Nint(pitch)
           elseif (isign.eq.1) then
               WRITE(FILENAME,'(''gr_max_x'',I1.1,''_q'',I2.2,
     &             ''_1.dat'')') imode,Nint(pitch)
           else
               WRITE(*,*) 'Warning: Uncorrect mode type (0 or 1)'
               stop
           endif
         open(10,file=filename)
         write(*,*) 'Isign=',isign,'  Nmode=',imode
         do fpfc=0.5D0,0.5d0,0.05d0
           write(*,*) 'fpfc=',fpfc
           fpe=fpfc*fce
           call maxgr(fpe,fce,fm,thetam,grm,imode)
           write(10,50) fpfc,grm,thetam,fm/fce,isign,imode
         enddo
50         format(1x,4(g12.4),2(i4))
         close(10)
      enddo
      end



      subroutine maxgr(fpe,fce,fm,thetam,grm,nh)
      implicit double precision (a-h,o-z)
          character*40 filename_contour
      parameter(nx=1001,ny=181)
      real*8 vgf(nx),vff(nx),vqf(nx)
      real*8 vgq(ny),vfq(ny),vqq(ny)
      common /para1/f,en2,isign
      common /para2/alfa,beta,up0,uz0,fc,anorm
      common /para_pitch/pitch
      fpe2=fpe*fpe
      fce2=fce*fce
      if (isign .eq. 0) then           !For O mode

          WRITE(filename_contour,'(''gr_O'',I1.1,''_pitch'',F4.1
     &        ,''_FpFc'',F4.2,''.dat'')') nh,pitch,fpe
          open(11,file=filename_contour)
          write(11,10) alfa,beta,up0,uz0,fpe
10        format(1x,5(f10.4))
          write(11,20) isign,nh,nx,ny
20        format(1x,2I4,2I8)

          f_cut=fpe
         if (nh .eq. 1) then      !For fundament wave
         f1=1.0001d0*f_cut
         f2=1.5d0*fce
        else if(nh .eq. 2 ) then      !For Harmomic wave
         f1=1.5d0*fce
         f2=2.5d0*fce
        else
         write(*,*) 'Warning: Uncorrect mode number (1 or 2)!'
         stop
        endif
      elseif (isign .eq. 1) then    !For X mode

          WRITE(filename_contour,'(''gr_X'',I1.1,''_pitch'',F4.1
     &        ,''_FpFc'',F4.2,''.dat'')') nh,pitch,fpe
          open(11,file=filename_contour)
          write(11,10) alfa,beta,up0,uz0,fpe
          write(11,20) isign,nh,nx,ny

          f_cut=(dsqrt(4.d0*fpe2+fce2)+fce)/2.d0
        if (nh .eq. 1) then      !For fundament wave
         f1=1.0001d0*f_cut
         f2=1.5d0*fce
        else if(nh .eq. 2 ) then      !For Harmomic wave
         f1=dmax1(1.5d0*fce,1.0001d0*f_cut)
         f2=2.5d0*fce
        else
         write(*,*) 'Warning: Uncorrect mode number (1 or 2)!'
         stop
        endif
      endif
      df=(f2-f1)/dble(nx-1)
      dq=180.d0/dble(ny-1)
      do iff=1,nx
        f=f1+dble(iff-1)*df
        do iq=1,ny
         theta=0.D0+dble(iq-1)*dq
               if (f .le. f_cut) then
             gr=0.d0
               else
                  call growth(theta,fpe,fce,gr)
               endif
         vgq(iq)=gr
         vfq(iq)=f
         vqq(iq)=theta
               write(11,30) f,theta,gr
30             format(1x,3(g15.7))
        enddo
        call piksrt(ny,vgq,vfq,vqq)
             vgf(iff)=vgq(ny)
             vff(iff)=vfq(ny)
             vqf(iff)=vqq(ny)
      enddo
      call piksrt(nx,vgf,vff,vqf)
      fm=vff(nx)
      grm=vgf(nx)
      thetam=vqf(nx)
      if (grm .le. 0.d0)then
         write(*,*) 'No emission'
      else
c         write(*,*) 'Isign=',isign,'  Nmode=',nh
         write(*,*) 's=f/fce=',fm/fce
         write(*,*) 'grm=',grm
         write(*,*) 'thetam=', thetam
         write(*,*)
      endif
      close(11)
      return
      end




c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+    Subroutine for calculating the growth rate of X or O mode.                                                         +
c+--------------------------------------------------------------------------------------------------------------------------------------+
c+      theta --> direction angle of the wave propagation                                                                         +
c+    fpe   --> plasma frequency of the back ground plasma                                                               +
c+    fce   --> electron gyrofrequency                                                                                                  +
c+    f     --> frequency of the eletromagnetic wave                                                                             +
c+    gr    --> growth rate of the mode                                                                                                +
c+                                                                                                                                                     +
c+    en    --> refraction index                                                                                                            +
c+    alfa  --> momentum dispersion in perpendicular direction                                                          +
c+    beta  --> momentum dispersion in parallel direction                                                                  +
c+    up0   --> a ring momentum                                                                                                       +
c+    uz0   --> the parrallel momentum                                                                                              +
c+    s     --> sin(theta)                                                                                                                      +
c+    c     --> cos(theta)                                                                                                                     +
c+                                                                                                                                                     +
c+    isign =  1 for X mode                                                                                                                  +
c+             0 for O mode                                                                                                                  +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine growth(theta,fpe,fce,gr)
      implicit double precision (a-h,o-z)
      parameter(eps=1.d-18)
      complex*16 func1,zint1,func2
      external func1,func2
      common /para1/f,en2,isign
      common /para2/alfa,beta,up0,uz0,fc,anorm
      common /para3/uc,aa,bb,s,c,nmode
      common /para4/ak,tt
      fpe2=fpe*fpe
      fce2=fce*fce
      fc=fce
      err=1.d-6
      lpmax=10      ! Used in integration calculation
      gr=0.d+0
      pi=4.d+0*datan(1.d0)
      pi2=pi*pi
      f2=f*f
      s=dsin(theta*pi/180.d0)
      c=dcos(theta*pi/180.d0)
      s2=s*s
      c2=c*c
      sq=f*fce*s2/dabs(2.d0*(f2-fpe2)+eps)
      sq2=sq*sq
      if (dabs(f) .le. fpe) tau=sq+dsqrt(c2+sq2+eps)
      if (dabs(f) .gt. fpe) tau=-sq-dsqrt(c2+sq2+eps)
      tau2=tau*tau
      uu=(tau2-c2)*(f2+fpe2)/((tau2+c2)*(f2-fpe2)+eps)
      if (isign .eq. 1) then      ! For X mode
      en2=1.d+0-fpe2/(f*(f+tau*fce)+eps)      !Nx^2, Nx is the refraction index
      if (en2 .le. 0.d+0) then
      gr=0.d0
      return
      endif
      if (en2 .ge. 1.d0) then
      gr=0.d0
      write(*,*)'Warning: the refraction index is great than 1!'
      write(*,*) en2
      stop
      endif
      en=dsqrt(en2+eps)      !Nx, the refraction index
      rr=1-tau*fpe2*fce*(1.d0+uu)/(2.d0*f*(f+tau*fce)**2+eps)      !Rx
      tt=-c/(tau+eps)      !Tx
      ak=fpe2*fce*s/((f2-fpe2)*(f+tau*fce)+eps)      ! Kx
      do nmode=1,3
      as=en2*c2-1.d0+dble(nmode)*dble(nmode)*fce2/f2
      if (as .le. 0.d0) then
      gr=gr+0.d0
      else
      aa2=as/(1.d0-en2*c2+eps)
      aa=dsqrt(aa2)
      bb2=aa2/(1.d0-en2*c2+eps)
      bb=dsqrt(bb2)
      uc=dble(nmode)*(fce/f)*en*c/(1.d0-en2*c2+eps)
      uz1=uc-bb
      uz2=uc+bb
*      write(*,*) 'Begin of the Integration'
      call cgcint(uz1,uz2,func2,zint1,err,lpmax) !integration of Q (define by func2)
*      write(*,*) 'End of the Integration'
      Gq=-2.d0*pi2*fpe2/(f*fce*rr*(1.d0+tt*tt)+eps)
      gr=gr+dreal(zint1)*Gq
      endif
      end do
      elseif (isign .eq. 0) then      ! For O mode
      en2=1.d0-tau*fpe2/(f*(tau*f-fce*c2)+eps)
      if (en2 .le. 0.d0) then
      gr=0.d0
      return
      endif
      if (en2 .ge. 1.d0) then
      gr=0.d0
      write(*,*)'Warning: the refraction index is great than 1!'
      stop
      endif
      en=dsqrt(en2+eps)      !No, the refraction index
      rr=1.d0+tau*fpe2*fce*c2*(1.d0-uu)      !Ro
     &          /(2.d0*f*(tau*f-fce*c2)**2+eps)
      tt=tau/(c+eps)      !To
      ak=tau*fpe2*fce*s/((f2-fpe2)*(tau*f-fce*c2)+eps)      !Ko
      do nmode=1,3
      as=en2*c2-1.d0+dble(nmode)*dble(nmode)*fce2/f2
      if (as .le. 0.d0) then
      gr=gr+0.d0
      else
      aa2=as/(1.d0-en2*c2+eps)
      aa=dsqrt(aa2)
      bb2=aa2/(1.d0-en2*c2+eps)
      bb=dsqrt(bb2)
      uc=dble(nmode)*(fce/f)*en*c/(1.d0-en2*c2+eps)
      uz1=uc-bb
      uz2=uc+bb
      call cgcint(uz1,uz2,func2,zint1,err,lpmax) !integration of Q (define by func2)
      Gq=-2.d0*pi2*fpe2/(f*fce*rr*(1.d0+tt*tt)+eps)
      gr=gr+dreal(zint1)*Gq
      endif
      end do
      endif
      return
      end


      subroutine piksrt(n,ara,arb,arc)
      implicit double precision (a-h,o-z)
      dimension ara(n),arb(n),arc(n)
      do j=2,n
      a=ara(j)
      do i=j-1,1,-1
      if(ara(i) .le. a)goto 10
      ara(i+1)=ara(i)
      arb(i+1)=arb(i)
      arc(i+1)=arc(i)
      enddo
      i=0
10      ara(i+1)=a
      enddo
      return
      end


      function ntimes(n)
      implicit double precision (a-h,o-z)
      integer n
      if (n .lt. 0) then
      write(*,*) 'ERROR: uncorrect input integer in function atimes!'
      stop
      endif
      if (n .eq. 0) then
      ntimes=1
      return
      endif
      ntimes=1
      do i=1,n
      ntimes=ntimes*i
      enddo
      return
      end

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+ For distribution function of a ring-beam without intrinsic Alfven waves,      +
c+ which can be expressed in the following form:                                 +
c+                   [                      2                    2    ]          +
c+                   [   ( u_perp - u_perp0)       ( u_z - u_z0 )     ]          +
c+     f_b(u)= A exp [- ---------------------- - -------------------  ]          +
c+                   [           2                         2          ]          +
c+                   [       Alfa                      Beta           ]          +
c+                                                                               +
c+ where, u denotes the momentum per unit mass,                                  +
c+ u_perp and u_z are the perpendicular and parallel momentum,                   +
c+ A, Alfa, Beta, u_perp0, u_z0 are constants.                                   +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      complex*16 function func1(x)      !x=uz
      implicit double precision (a-h,o-z)
      parameter(pi=3.1415926535898d+00,eps=1.d-18)
      common /para1/f,en2,isign
      common /para2/alfa,beta,up0,uz0,fc,anorm
      common /para3/uc,aa,bb,s,c,nmode
      common /para4/ak,tt
*      write(*,*) 'Calling func1'
      if ((x .gt. uc+bb) .or. (x .lt. uc-bb)) then
      func1=0.d0
      return
      endif
      en=dsqrt(en2)
      alfa2=alfa*alfa
      beta2=beta*beta
      upr=dabs(aa)*dsqrt(dabs(bb*bb-(x-uc)**2))/(dabs(bb)+eps)
      const1=anorm
      const2=(dble(nmode)*fc/f+en*c*x)
      const3=upr*dexp(-(upr-up0)**2/alfa2-(x-uz0)**2/beta2)
      const4=(upr-up0)*(1.d0-en*c*x)/alfa2+(x-uz0)*en*c*upr/beta2
c
      barg=en*upr*s*(f/fc)
      hn=(f/fc)*(ak*s+tt*(c-en*x))*bjox(nmode,barg)+bjp(nmode,barg)
c
      qfunc=hn*hn*const2*const3*const4/const1
      qfunc=hn*hn*const2*const3*const4/const1
      func1=dcmplx(qfunc,0.d0)
      return
      end


c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+ For distribution function of a ring-beam scattered by intrinsic Alfven waves, +
c+ which can be expressed in the following form:                                 +
c+                   [              2                       2         ]          +
c+                   [   ( u - u_0 )         ( miu - miu_0 )          ]          +
c+     f_b(u)= A exp [- -------------  -   -------------------------- ]          +
c+                   [           2                    2               ]          +
c+                   [       Alfa                 Beta                ]          +
c+                                                                               +
c+ where, u denotes the momentum per unit mass, miu = u_z / u = cos (theta_p),   +
c+ and u_z is parallel momentum along the backgorund magnetic field line,        +
c+ theta_p is the pitch-angle. A,  Alfa and Beta are constants.                  +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      complex*16 function func2(x)      !x=uz
      implicit double precision (a-h,o-z)
      parameter(pi=3.1415926535898d+00,eps=1.d-18)
      common /para1/f,en2,isign
      common /para2/alfa,beta,up0,uz0,fc,anorm
      common /para3/uc,aa,bb,s,c,nmode
      common /para4/ak,tt
*      write(*,*) 'Calling func2'
      if ((x .gt. uc+bb) .or. (x .lt. uc-bb)) then
      func1=0.d0
      return
      endif
      u0=dsqrt(up0*up0+uz0*uz0)
      amiu0=uz0/u0
      en=dsqrt(en2)
      alfa2=alfa*alfa
      beta2=beta*beta
      upr=dabs(aa)*dsqrt(dabs(bb*bb-(x-uc)**2))/(dabs(bb)+eps)
      u=dsqrt(upr*upr+x*x)
      amiu=x/u
      const1=anorm
      const2=(dble(nmode)*fc/f+en*c*x)
      const3=(upr**2/u**2)*dexp(-(u-u0)**2/alfa2-(amiu-amiu0)**2/beta2)
      const4=u*(u-u0)/alfa2+(en*c*u-amiu)*(amiu-amiu0)/beta2
c
      barg=en*upr*s*(f/fc)
      hn=(f/fc)*(ak*s+tt*(c-en*x))*bjox(nmode,barg)+bjp(nmode,barg)
c
      qfunc=hn*hn*const2*const3*const4/const1
      func2=dcmplx(qfunc,0.d0)
      return
      end


      REAl*8 FUNCTION Anormal1(up0,uz0,alfa,beta)
c       Return the normalization constant for a ring-beam distribution function
c       without IAW, namely, func1
      implicit double precision (a-h,o-z)
      parameter(pi=3.1415926535898d+00)
      print*,up0,uz0,alfa,beta
      alfa2=alfa*alfa
      Anormal1=pi*dsqrt(pi)*alfa2*beta
     &  *(dexp(-up0*up0/alfa2)+dsqrt(pi)*up0*(1.d0+errf(up0/alfa))/alfa)
      return
      END


      REAl*8 FUNCTION Anormal2(up0,uz0,alfa,beta)
c       Return the normalization constant for a ring-beam distribution function
c       with IAW, namely, func2
      implicit double precision (a-h,o-z)
      parameter(pi=3.1415926535898d+00)
      u0=dsqrt(up0*up0+uz0*uz0)
      amiu0=uz0/u0
      z0=-dsqrt(2.d0)*u0/alfa
      d1=dexp(z0*z0/4.0d0)*dsqrt(pi/2.0d0)*errfc(z0/dsqrt(2.d0))
      d2=dexp(z0*z0/4.0d0)*dsqrt(pi/2.0d0)*
     &   (dsqrt(pi/2.0d0)*dexp(-z0*z0/2.d0)-z0*errfc(z0/dsqrt(2.d0)))
      d3=0.5d0*(d1-z0*d2)
      Anormal2=dsqrt(pi/2.d0)*pi*beta*alfa**3*dexp(-0.5d0*u0*u0/alfa**2)
     &     *d3*(errf((1.d0-amiu0)/beta)+errf((1.d0+amiu0)/beta))
      return
      END


      subroutine cgcint(a,b,funct,zint,err,lpmax)
c
c     CGCINT: Complex integral along real axis.
c     "Borrowed" from Joe Huba.
c
      implicit double precision (a-h,o-z)
      complex*16 one3rd,s,funct,ols,zint
      dimension np1(6)
      real*8 r(364),w(364)
      logical first
      save
      data first/.true./
c
      if(first)then
        z=0.0d0
        pi=1.0d0/3.0d0
        one3rd=dcmplx(pi,z)
        pi=4.0d0*datan(1.d0)
        np1(1)=6
        k=1
        w(1)=dsin(pi/12.0d0)
        r(1)=dcos(pi/12.0d0)
        do lps=2,6
          np1(lps)=3*np1(lps-1)
          pd2np1=0.5d0*pi/float(np1(lps))
          k=k+1
          w(k)=dsin(pd2np1)
          r(k)=dcos(pd2np1)
          iup=(np1(lps)/2-3)/2
          do i=3,iup,3
            do j=1,2
              a1=float(2*(i+j)-3)*pd2np1
              k=k+1
              w(k)=dsin(a1)
              r(k)=dcos(a1)
            end do
          end do
        end do
        first=.false.
      endif
c
      ctr=0.5d0*(a+b)
      del=0.5d0*(b-a)
      lpmax=max0(lpmax,1)
      lpmax=min0(lpmax,6)
c
      lps=0
      k=lps
      x=dsqrt(0.5d0)
      s=dcmplx(0.5d0*x,z)*(funct(ctr+del*x)+funct(ctr-del*x))
c
3      lps=lps+1
      ols=s
      s=dcmplx(0.0d0,0.0d0)
c      iup=np1(lps-1)/2
      if(lps .eq. 1)then
        iup=1
      else
        iup=np1(lps-1)/2
      endif
      do 4 i=1,iup
        k=k+1
        x1=del*r(k)
        x2=del*w(k)
        s=s+dcmplx(w(k),z)*(funct(ctr+x1)+funct(ctr-x1))+
     *          dcmplx(r(k),z)*(funct(ctr+x2)+funct(ctr-x2))
4      continue
      oneo=1.0d0/float(np1(lps))
      s=one3rd*ols+dcmplx(oneo,z)*s
      t1=cdabs(s-ols)
      t2=err*cdabs(s)
      if((t1 .gt. t2) .and. (lps .lt. lpmax))goto 3
      zint=dcmplx(pi*del,z)*s
      return
      end


      REAL*8 FUNCTION errf(x)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*      write(*,*) 'Calling Errf'
*       Return the error function erf(x)
      IF(X.LT.0.D0)THEN
        ERRF=-GAMMP(0.5D0,X**2)
      ELSE
        ERRF=GAMMP(0.5D0,X**2)
      ENDIF
      RETURN
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.


      REAL*8 FUNCTION ERRFC(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*      write(*,*) 'Calling Errfc'
*       Return the complementary error function erfc(x)=1.d0-erf(x)
      IF (X .LT. 0.D0) THEN
      ERRFC=1.D0+GAMMP(0.5D0,X**2)
      ELSE
      ERRFC=GAMMQ(0.5D0,X**2)
      ENDIF
      RETURN
      END


      REAL*8 FUNCTION GAMMP(A,X)
      implicit double precision (a-h,o-z)
       IF((X .LT. 0.D0) .OR. (A .LE. 0.D0))PAUSE
      IF(X .LT. A+1.D0)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMP=GAMSER
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMP=1.D0-GAMMCF
      ENDIF
      RETURN
      END

      REAL*8 FUNCTION GAMMQ(A,X)
      implicit double precision (a-h,o-z)
       IF((X .LT. 0.D0) .OR. (A .LE. 0.D0))PAUSE
      IF(X .LT. A+1.D0)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMQ=1.D0-GAMSER
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMQ=GAMMCF
      ENDIF
      RETURN
      END

      SUBROUTINE GSER(GAMSER,A,X,GLN)
      implicit double precision (a-h,o-z)
      PARAMETER(ITMAX=100,EPS=3.0D-7)
      GLN=GAMMLN(A)
      IF(X .LE. 0.D0)THEN
        IF(X .LT. 0.D0)PAUSE
        GAMSER=0.D0
        RETURN
      ENDIF
      AP=A
      SUM=1.D0/A
      DEL=SUM
      DO 1 N=1,ITMAX
        AP=AP+1.D0
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(DABS(DEL) .LT. DABS(SUM)*EPS)GOTO 190
1      CONTINUE
      PAUSE 'A too large, ITMAX too small'
190      GAMSER=SUM*DEXP(-X+A*DLOG(X)-GLN)
      RETURN
      END

      SUBROUTINE GCF(GAMMCF,A,X,GLN)
      implicit double precision (a-h,o-z)
      PARAMETER(ITMAX=100,EPS=3.0D-7)
      GLN=GAMMLN(A)
      GOLD=0.D0
      A0=1.D0
      A1=X
      B0=0.D0
      B1=1.D0
      FAC=1.D0
      DO 1 N=1,ITMAX
        AN=N
        ANA=AN-A
        A0=(A1+A0*ANA)*FAC
        B0=(B1+B0*ANA)*FAC
        ANF=AN*FAC
        A1=X*A0+ANF*A1
        B1=X*B0+ANF*B1
        IF(A1 .NE. 0.D0) THEN
          FAC=1.D0/A1
          G=B1*FAC
          IF(DABS((G-GOLD)/G) .LT. EPS)GOTO 180
          GOLD=G
        ENDIF
1      CONTINUE
      PAUSE 'A too large, ITMAX too small'
180      GAMMCF=DEXP(-X+A*DLOG(X)-GLN)*G
      RETURN
      END

      REAL*8 FUNCTION GAMMLN(XX)
      implicit double precision (a-h,o-z)
      DIMENSION COF(6)
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *      -1.231739516D0,.120858003D-2,-5.36382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*DLOG(TMP)-TMP
      SER=ONE
      DO 1 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
1      CONTINUE
      GAMMLN=TMP+DLOG(STP*SER)
      RETURN
      END

      REAL*8 FUNCTION BJP(N,X)
      implicit double precision (a-h,o-z)
      BJP=0.5*(BJ(N-1,X)-BJ(N+1,X))
      RETURN
      END

      REAL*8 FUNCTION BJOX(N,X)
      implicit double precision (a-h,o-z)
      IF(N .EQ. 0)THEN
      BJOX=BJ(N,X)/X
      ELSE
      BJOX=0.5*(BJ(N-1,X) + BJ(N+1,X))/N
      ENDIF
      RETURN
      END

      REAL*8 FUNCTION BJ(N,X)
      implicit double precision (a-h,o-z)
      PARAMETER(IACC=40,BIGNO=1.D10,BIGNI=1.D-10)
      IF(N .EQ. 0)BJ=BESSJ0(X)
      IF(N .EQ. 1)BJ=BESSJ1(X)
      IF(N .LT. 2)RETURN
      AX=DABS(X)
      IF(AX .EQ. 0.D0)THEN
        BJ=0.D0
      ELSEIF(AX .GT. FLOAT(N))THEN
        TOX=2.D0/AX
        BJM=BESSJ0(AX)
        BJT=BESSJ1(AX)
        DO 1 J=1,N-1
        BJP=J*TOX*BJT-BJM
        BJM=BJT
1        BJT=BJP
        BJ=BJT
      ELSE
        TOX=2.D0/AX
        M=2*((N+INT(SQRT(FLOAT(IACC*N))))/2)
        BJ=0.D0
        JSUM=0.D0
        SUM=0.D0
        BJP=0.D0
        BJT=1.D0
        DO 2 J=M,1,-1
        BJM=J*TOX*BJT-BJP
        BJP=BJT
        BJT=BJM
        IF(DABS(BJT) .GT. BIGNO)THEN
            BJT=BJT*BIGNI
            BJP=BJP*BIGNI
            BJ=BJ*BIGNI
            SUM=SUM*BIGNI
          ENDIF
          IF(JSUM .NE. 0)SUM=SUM+BJT
          JSUM=1-JSUM
          IF(J .EQ. N)BJ=BJP
2        CONTINUE
        SUM=2.*SUM-BJT
        BJ=BJ/SUM
      ENDIF
      IF((X .LT. 0.D0) .AND. (MOD(N,2) .EQ. 1))BJ=-BJ
      RETURN
      END

      REAL*8 FUNCTION BESSJ0(X)
      implicit double precision (a-h,o-z)
      DATA P1,P2,P3,P4,P5/1.D0,-.1098628627D-2,.2734510407D-4,
     *	    -.2073370639D-5,.2093887211D-6/
      DATA Q1,Q2,Q3,Q4,Q5/-.1562499995D-1,.1430488765D-3,
     *	    -.6911147651D-5,.7621095161D-6,-.934945152D-7/
      DATA R1,R2,R3,R4,R5,R6/57568490574.D0,-13362590354.D0,
     *	    651619640.7D0,-11214424.18D0,77392.33017D0,-184.9052456D0/
      DATA S1,S2,S3,S4,S5,S6/57568490411.D0,1029532985.D0,9494680.718D0,
     *	    59272.64853D0,267.8532712D0,1.D0/
      IF(DABS(X) .LT. 8.D0)THEN
        Y=X**2
        BESSJ0=(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))
     *          /(S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6)))))
      ELSE
        AX=DABS(X)
        Z=8.D0/AX
        Y=Z**2
        XX=AX-.785398164D0
        BESSJ0=DSQRT(.636619772D0/AX)*(DCOS(XX)
     *          *(P1+Y*(P2+Y*(P3+Y*(P4+Y*P5))))
     *          -Z*DSIN(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))
      ENDIF
      RETURN
      END

      REAL*8 FUNCTION BESSJ1(X)
      implicit double precision (a-h,o-z)
      DATA R1,R2,R3,R4,R5,R6/72362614232.D0,-7895059235.D0,
     *	    242396853.1D0,-2972611.439D0,15704.48260D0,-30.16036606D0/
      DATA S1,S2,S3,S4,S5,S6/144725228442.D0,2300535178.D0,
     *	    18583304.74D0,99447.43394D0,376.9991397D0,1.D0/
      DATA P1,P2,P3,P4,P5/1.D0,.183105D-2,-.3516396496D-4,
     *	    .2457520174D-5,-.240337019D-6/
      DATA Q1,Q2,Q3,Q4,Q5/.04687499995D0,-.2002690873D-3,
     *	    .8449199096D-5,-.88228987D-6,.105787412D-6/
      IF(DABS(X) .LT. 8.D0)THEN
        Y=X**2
        BESSJ1=X*(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))
     *          /(S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6)))))
      ELSE
        AX=DABS(X)
        Z=8.D0/AX
        Y=Z**2
        XX=AX-2.356194491D0
        BESSJ1=DSQRT(.636619772D0/AX)*(DCOS(XX)*
     *	   (P1+Y*(P2+Y*(P3+Y*(P4+Y*P5))))
     *	   -Z*DSIN(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))*DSIGN(1.D0,X)
      ENDIF
      RETURN
      END
