***
      SUBROUTINE deltat(kw,age,tm,tn,tscls,dt,dtr)
      implicit none
*
      INTEGER kw
      REAL*8 age,tm,tn,tscls(20)
      REAL*8 dt,dtr
      REAL*8 pts1,pts2,pts3
      COMMON /POINTS/ pts1,pts2,pts3
*
*     Base new time scale for changes in radius & mass on stellar type.
*
      if(kw.le.1)then
         dt = pts1*tm
         dtr = tm - age
      elseif(kw.eq.2)then
         dt = pts1*(tscls(1) - tm)
         dtr = tscls(1) - age
      elseif(kw.eq.3)then
         if(age.lt.tscls(6))then
            dt = pts2*(tscls(4) - age)
         else
            dt = pts2*(tscls(5) - age)
         endif
         dtr = MIN(tscls(2),tn) - age
      elseif(kw.eq.4)then
         dt = pts2*tscls(3)
         dtr = MIN(tn,tscls(2) + tscls(3)) - age
      elseif(kw.eq.5)then
         if(age.lt.tscls(9))then
            dt = pts3*(tscls(7) - age)
         else
            dt = pts3*(tscls(8) - age)
         endif
         dtr = MIN(tn,tscls(13)) - age
      elseif(kw.eq.6)then
         if(age.lt.tscls(12))then
            dt = pts3*(tscls(10) - age)
         else
            dt = pts3*(tscls(11) - age)
         endif
         dt = MIN(dt,0.005d0)
         dtr = tn - age
      elseif(kw.eq.7)then
         dt = pts1*tm
         dtr = tm - age
      elseif(kw.eq.8.or.kw.eq.9)then
         if(age.lt.tscls(6))then
            dt = pts2*(tscls(4) - age)
         else
            dt = pts2*(tscls(5) - age)
         endif
         dtr = tn - age
      else
*        dt = MAX(0.1d0,age*10.d0)
         dt = MAX(0.1d0,dt*10.d0)
         dt = MIN(dt,5.0d+02)
         dtr = dt
      endif
*
      RETURN
      END
***

***
      SUBROUTINE evolv1(kw,mass,mt,r,lum,mc,rc,menv,renv,ospin,
     &                  epoch,tm,tphys,tphysf,dtp,z,zpars)
c-------------------------------------------------------------c
c
c     Evolves a single star.
c     Mass loss is an option.
c     The timestep is not constant but determined by certain criteria.
c     Plots the HRD and variables as a function of time.
c
c     Written by Jarrod Hurley 26/08/97 at the Institute of
c     Astronomy, Cambridge.
c
c-------------------------------------------------------------c
c
c     STELLAR TYPES - KW
c
c        0 - deeply or fully convective low mass MS star
c        1 - Main Sequence star
c        2 - Hertzsprung Gap
c        3 - First Giant Branch
c        4 - Core Helium Burning
c        5 - First Asymptotic Giant Branch
c        6 - Second Asymptotic Giant Branch
c        7 - Main Sequence Naked Helium star
c        8 - Hertzsprung Gap Naked Helium star
c        9 - Giant Branch Naked Helium star
c       10 - Helium White Dwarf
c       11 - Carbon/Oxygen White Dwarf
c       12 - Oxygen/Neon White Dwarf
c       13 - Neutron Star
c       14 - Black Hole
c       15 - Massless Supernova
c
c-------------------------------------------------------------c
      implicit none
*
      integer kw,it,ip,jp,j,kwold,rflag
      integer nv
      parameter(nv=50000)
*
      real*8 mass,z,aj
      real*8 epoch,tphys,tphys2,tmold,tbgold
      real*8 mt,tm,tn,tphysf,dtp,tsave
      real*8 tscls(20),lums(10),GB(10),zpars(20)
      real*8 r,lum,mc,teff,rc,menv,renv,vs(3)
      real*8 ospin,jspin,djt,djmb,k2,k3
      parameter(k3=0.21d0)
      real*8 m0,r1,lum1,mc1,rc1,menv1,renv1,k21
      real*8 dt,dtm,dtr,dr,dtdr,dms,dml,mt2,rl
      real*8 tol,tiny
      parameter(tol=1.0d-10,tiny=1.0d-14)
      real*8 ajhold,rm0,eps,alpha2
      parameter(eps=1.0d-06,alpha2=0.09d0)
      real*8 mlwind,vrotf
      external mlwind,vrotf
      logical iplot,isave
      REAL*8 neta,bwind,hewind,mxns
      COMMON /VALUE1/ neta,bwind,hewind,mxns
      REAL*8 pts1,pts2,pts3
      COMMON /POINTS/ pts1,pts2,pts3
      REAL scm(50000,14),spp(20,3)
      COMMON /SINGLE/ scm,spp
*
      dtm = 0.d0
      r = 0.d0
      lum = 0.d0
      mc = 0.d0
      mc1 = 0.d0
      rc = 0.d0
      rl = 0.d0
      if(ospin.le.0.d0)then
         ospin = 1.0d-10
         jspin = 1.0d-10
      endif
      k2 = 0.15d0
      rflag = 0
*
* Setup variables which control the output (if it is required).
*
      ip = 0
      jp = 0
      tsave = tphys
      isave = .true.
      iplot = .false.
      if(dtp.le.0.d0)then
         iplot = .true.
         isave = .false.
         tsave = tphysf
      elseif(dtp.gt.tphysf)then
         isave = .false.
         tsave = tphysf
      endif
* 
      do 10 , j = 1,nv
*
         if(neta.gt.tiny.and.j.gt.1)then
*
* Calculate mass loss from the previous timestep.
*
            dt = 1.0d+06*dtm
            dms = mlwind(kw,lum,r,mt,mc,rl,z)*dt
            if(kw.lt.10)then
               dml = mt - mc
               if(dml.lt.dms)then
                  dtm = (dml/dms)*dtm
                  dms = dml
               endif
            endif
         else
            dms = 0.d0
         endif
*
* Limit to 1% mass loss.
*
         if(dms.gt.0.01d0*mt)then
            dtm = 0.01d0*mt*dtm/dms
            dms = 0.01d0*mt
         endif
*
* Calculate the rate of angular momentum loss due to magnetic braking 
* and/or mass loss.
*
         if(j.gt.1)then
            djt = (2.d0/3.d0)*(dms/(1.0d+06*dtm))*r*r*ospin
            if(mt.gt.0.35d0.and.kw.lt.10)then
               djmb = 5.83d-16*menv*(r*ospin)**3/mt
               djt = djt + djmb
            endif
         endif
*
* Update mass and time and reset epoch for a MS (and possibly a HG) star.
*
         if(dms.gt.0.d0)then
            mt = mt - dms
            if(kw.le.2.or.kw.eq.7)then
               m0 = mass
               mc1 = mc
               mass = mt
               tmold = tm
               tbgold = tscls(1)
               CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
               if(kw.eq.2)then
                  if(GB(9).lt.mc1.or.m0.gt.zpars(3))then
                     mass = m0
                  else
                     epoch = tm + (tscls(1) - tm)*(ajhold-tmold)/
     &                            (tbgold - tmold)
                     epoch = tphys - epoch
                  endif
               else
                  epoch = tphys - ajhold*tm/tmold
               endif
            endif
         endif
         tphys2 = tphys
         tphys = tphys + dtm
*
* Find the landmark luminosities and timescales as well as setting
* the GB parameters.
*
         aj = tphys - epoch
         CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
*
* Find the current radius, luminosity, core mass and stellar type
* given the initial mass, current mass, metallicity and age
*
         kwold = kw
         CALL hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &               r,lum,kw,mc,rc,menv,renv,k2)
*
* If mass loss has occurred and no type change then check that we
* have indeed limited the radius change to 10%.
*
         if(kw.eq.kwold.and.dms.gt.0.d0.and.rflag.ne.0)then
            mt2 = mt + dms
            dml = dms/dtm
            it = 0
 20         dr = r - rm0
            if(ABS(dr).gt.0.1d0*rm0)then
               it = it + 1
               if(it.eq.20.and.kw.eq.4) goto 30
               if(it.gt.30)then
                  CALL exit(0)
                  STOP 
               endif
               dtdr = dtm/ABS(dr)
               dtm = alpha2*MAX(rm0,r)*dtdr
               if(it.ge.20) dtm = 0.5d0*dtm
               if(dtm.lt.1.0d-07*aj) goto 30
               dms = dtm*dml
               mt = mt2 - dms
               if(kw.le.2.or.kw.eq.7)then
                  mass = mt
                  CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
                  if(kw.eq.2)then
                     if(GB(9).lt.mc1.or.m0.gt.zpars(3))then
                        mass = m0
                     else
                        epoch = tm + (tscls(1) - tm)*(ajhold-tmold)/
     &                               (tbgold - tmold)
                        epoch = tphys2 - epoch
                     endif
                  else
                     epoch = tphys2 - ajhold*tm/tmold
                  endif
               endif
               tphys = tphys2 + dtm
               aj = tphys - epoch
               mc = mc1
               CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
               CALL hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &                     r,lum,kw,mc,rc,menv,renv,k2)
               goto 20
            endif
 30         continue
         endif
*
* Initialize or adjust the spin of the star.
*
         if(j.eq.1)then
            if(tphys.lt.tiny.and.ospin.lt.0.001d0)then
               ospin = 45.35d0*vrotf(mt)/r
            endif
            jspin = ospin*(k2*r*r*(mt-mc)+k3*rc*rc*mc)
         else
            jspin = MAX(1.0d-10,jspin - djt*1.0d+06*dtm)
            ospin = jspin/(k2*r*r*(mt-mc)+k3*rc*rc*mc)
         endif
*
* Test for changes in evolution type.
*
         if(j.eq.1.or.kw.ne.kwold)then
*
* Force new NS or BH to have a one second period. 
* 
            if(kw.eq.13.or.kw.eq.14)then
               ospin = 2.0d+08
               jspin = k3*rc*rc*mc*ospin
               CALL kick(kw,mass,mt,0.d0,0.d0,-1.d0,0.d0,vs)
            endif
            jp = jp + 1
            spp(jp,1) = tphys
            spp(jp,2) = float(kw)
            if(kw.eq.15)then
               spp(jp,3) = mass 
               goto 90
            else
               spp(jp,3) = mt
            endif
         endif
*
* Record values for plotting and reset epoch.
*
         epoch = tphys - aj
         if((isave.and.tphys.ge.tsave).or.iplot)then
            ip = ip + 1
            scm(ip,1) = tphys
            scm(ip,2) = float(kw)
            scm(ip,3) = mass
            scm(ip,4) = mt
            scm(ip,5) = log10(lum)
            scm(ip,6) = log10(r)
            teff = 1000.d0*((1130.d0*lum/(r**2.d0))**(1.d0/4.d0))
            scm(ip,7) = log10(teff)
            scm(ip,8) = mc
            scm(ip,9) = rc
            scm(ip,10) = menv
            scm(ip,11) = renv
            scm(ip,12) = epoch
            scm(ip,13) = ospin
            if(isave) tsave = tsave + dtp
            if(tphysf.lt.tiny)then
               ip = ip + 1
               do 35 , it = 1,13
                  scm(ip,it) = scm(ip-1,it)
 35            continue
            endif
         endif
*
         if(tphys.ge.tphysf)then
            jp = jp + 1
            spp(jp,1) = tphys
            spp(jp,2) = float(kw)
            spp(jp,3) = mt
            goto 90
         endif
*
* Record radius and current age.
*
         rm0 = r
         ajhold = aj
         if(kw.ne.kwold) kwold = kw
         CALL deltat(kw,aj,tm,tn,tscls,dtm,dtr)
*
* Check for type change.
*
         it = 0
         m0 = mass
         if((dtr-dtm).le.tol.and.kw.le.9)then
*
* Check final radius for too large a jump.
*
            aj = MAX(aj,aj*(1.d0-eps)+dtr)
            mc1 = mc 
            CALL hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &                  r1,lum1,kw,mc1,rc1,menv1,renv1,k21)
            dr = r1 - rm0
            if(ABS(dr).gt.0.1d0*rm0)then
               dtm = dtr - ajhold*eps
               dtdr = dtm/ABS(dr)
               dtm = alpha2*MAX(r1,rm0)*dtdr
               goto 40
            else
               dtm = dtr
               goto 50
            endif
         endif
*
* Limit to a 10% increase in radius assuming no further mass loss
* and thus that the pertubation functions due to small envelope mass
* will not change the radius.
*
 40      aj = ajhold + dtm
         mc1 = mc 
         CALL hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &               r1,lum1,kw,mc1,rc1,menv1,renv1,k21)
         dr = r1 - rm0
         it = it + 1
         if(it.eq.20.and.kw.eq.4) goto 50
         if(it.gt.30)then
            CALL exit(0)
            STOP 
         endif
         if(ABS(dr).gt.0.1d0*rm0)then
            dtdr = dtm/ABS(dr)
            dtm = alpha2*MAX(rm0,r1)*dtdr
            if(it.ge.20) dtm = 0.5d0*dtm
            goto 40
         endif
*
 50      continue
*
* Ensure that change of type has not occurred during radius check. 
* This is rare but may occur for HG stars of ZAMS mass > 50 Msun. 
*
         if(kw.ne.kwold)then
            kw = kwold
            mass = m0
            CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
         endif
*
* Choose minimum of time-scale and remaining interval (> 100 yrs).
*
         dtm = MAX(dtm,1.0d-07*aj)
         dtm = MIN(dtm,tsave-tphys)
*
 10   continue
*
 90   continue
*
      tphysf = tphys
      scm(ip+1,1) = -1.0
      spp(jp+1,1) = -1.0
      if(ip.ge.nv)then
         CALL exit(0)
         STOP
      endif
*
      RETURN
      END
***

***
      SUBROUTINE hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &                  r,lum,kw,mc,rc,menv,renv,k2)
*
*
*       H-R diagram for population I stars.
*       -----------------------------------
*
*       Computes the new mass, luminosity, radius & stellar type.
*       Input (MASS, AJ, TM, TN, LUMS & TSCLS) supplied by routine STAR.
*       Ref: P.P. Eggleton, M.J. Fitchett & C.A. Tout (1989) Ap.J. 347, 998.
*
*       Revised 27th March 1995 by C. A. Tout;
*       24th October 1995 to include metallicity;
*       14th November 1996 to include naked helium stars;
*       28th February 1997 to allow accretion induced supernovae.
*
*       Revised 5th April 1997 by J. R. Hurley
*       to include Z=0.001 as well as Z=0.02, convective overshooting,
*       MS hook and more elaborate CHeB
*
      implicit none
*
      integer kw,kwp
      INTEGER ceflag,tflag,ifflag,nsflag,wdflag
      COMMON /FLAGS/ ceflag,tflag,ifflag,nsflag,wdflag
*
      real*8 mass,aj,mt,tm,tn,tscls(20),lums(10),GB(10),zpars(20)
      real*8 r,lum,mc,rc,menv,renv,k2
      real*8 mch,mlp,tiny
      parameter(mch=1.44d0,mlp=12.d0,tiny=1.0d-14)
      real*8 mass0,mt0,mtc
      REAL*8 neta,bwind,hewind,mxns
      COMMON /VALUE1/ neta,bwind,hewind,mxns
* 
      real*8 thook,thg,tbagb,tau,tloop,taul,tauh,tau1,tau2,dtau,texp
      real*8 lx,ly,dell,alpha,beta,eta
      real*8 rx,ry,delr,rzams,rtms,gamma,rmin,taumin,rg
      parameter(taumin=5.0d-08)
      real*8 mcmax,mcx,mcy,mcbagb,lambda
      real*8 am,xx,fac,rdgen,mew,lum0,kap,zeta,ahe,aco
      parameter(lum0=7.0d+04,kap=-0.5d0,ahe=4.d0,aco=16.d0)
*
      real*8 thookf,tblf
      real*8 lalphf,lbetaf,lnetaf,lhookf,lgbtf,lmcgbf,lzhef,lpertf
      real*8 rzamsf,rtmsf,ralphf,rbetaf,rgammf,rhookf
      real*8 rgbf,rminf,ragbf,rzahbf,rzhef,rhehgf,rhegbf,rpertf
      real*8 mctmsf,mcgbtf,mcgbf,mcheif,mcagbf,lzahbf
      external thookf,tblf
      external lalphf,lbetaf,lnetaf,lhookf,lgbtf,lmcgbf,lzhef,lpertf
      external rzamsf,rtmsf,ralphf,rbetaf,rgammf,rhookf
      external rgbf,rminf,ragbf,rzahbf,rzhef,rhehgf,rhegbf,rpertf
      external mctmsf,mcgbtf,mcgbf,mcheif,mcagbf,lzahbf
*
*
*       ---------------------------------------------------------------------
*       MASS    Stellar mass in solar units (input: old; output: new value).
*       AJ      Current age in Myr.
*       MT      Current mass in solar units (used for R).
*       TM      Main sequence time.
*       TN      Nuclear burning time.
*       TSCLS   Time scale for different stages.
*       LUMS    Characteristic luminosity.
*       GB      Giant Branch parameters
*       ZPARS   Parameters for distinguishing various mass intervals.
*       R       Stellar radius in solar units.
*       TE      Effective temperature (suppressed).
*       KW      Classification type (0 - 15).
*       MC      Core mass.
*       ---------------------------------------------------------------------
*
*
* Make evolutionary changes to stars that have not reached KW > 5.
*
      mass0 = mass
      if(mass0.gt.100.d0) mass = 100.d0
      mt0 = mt
      if(mt0.gt.100.d0) mt = 100.d0
*
      if(kw.gt.6) goto 90
*
      tbagb = tscls(2) + tscls(3)
      thg = tscls(1) - tm
*
      rzams = rzamsf(mass)
      rtms = rtmsf(mass)
*
      if(aj.lt.tscls(1))then
*
*        Either on MS or HG
*
         rg = rgbf(mt,lums(3))
*
         if(aj.lt.tm)then
*
*           Main sequence star.
*
            mc = 0.d0
            tau = aj/tm
            thook = thookf(mass)*tscls(1)
            zeta = 0.01d0
            tau1 = MIN(1.d0,aj/thook)
            tau2 = MAX(0.d0,
     &             MIN(1.d0,(aj-(1.d0-zeta)*thook)/(zeta*thook)))
*
            dell = lhookf(mass,zpars(1))
            dtau = tau1**2 - tau2**2
            alpha = lalphf(mass)
            beta = lbetaf(mass)
            eta = lnetaf(mass)
            lx = LOG10(lums(2)/lums(1))
            if(tau.gt.taumin)then
               xx = alpha*tau + beta*tau**eta +
     &              (lx - alpha - beta)*tau**2 - dell*dtau
            else
               xx = alpha*tau + (lx - alpha)*tau**2 - dell*dtau
            endif
            lum = lums(1)*10.d0**xx
*
            delr = rhookf(mass,zpars(1))
            dtau = tau1**3 - tau2**3
            alpha = ralphf(mass)
            beta = rbetaf(mass)
            gamma = rgammf(mass)
            rx = LOG10(rtms/rzams)
* Note that the use of taumin is a slightly pedantic attempt to
* avoid floating point underflow. It IS overkill!
            if(tau.gt.taumin)then
               xx = alpha*tau + beta*tau**10 + gamma*tau**40 +
     &              (rx - alpha - beta - gamma)*tau**3 - delr*dtau
            else
               xx = alpha*tau + (rx - alpha)*tau**3 - delr*dtau
            endif
            r = rzams*10.d0**xx
*
            if(mass.lt.(zpars(1)-0.3d0))then
               kw = 0
* This following is given by Chris for low mass MS stars which will be 
* substantially degenerate. We need the Hydrogen abundance, X, which we 
* calculate from Z assuming that the helium abundance, Y, is calculated
* according to Y = 0.24 + 2*Z
               rdgen = 0.0258d0*((1.d0+zpars(11))**(5.d0/3.d0))*
     &                          (mass**(-1.d0/3.d0))
               r = MAX(rdgen,r)
            else
               kw = 1
            endif
*
         else 
*
*           Star is on the HG
*
            mcx = mc
            if(mass.le.zpars(2))then
               mc = mcgbf(lums(3),GB,lums(6))
            elseif(mass.le.zpars(3))then
               mc = mcheif(mass,zpars(2),zpars(9))
            else
               mc = mcheif(mass,zpars(2),zpars(10))
            endif
            eta = mctmsf(mass)
            tau = (aj - tm)/thg
            mc = ((1.d0 - tau)*eta + tau)*mc
            mc = MAX(mc,mcx)
*
* Test whether core mass has reached total mass.
*
            if(mc.ge.mt)then
               aj = 0.d0
               if(mass.gt.zpars(2))then
*
* Zero-age helium star
*
                  mc = 0.d0
                  mass = mt
                  kw = 7
                  CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
               else
*
* Zero-age helium white dwarf.
*
                  mc = mt
                  mass = mt
                  kw = 10
               endif
            else
               lum = lums(2)*(lums(3)/lums(2))**tau
               if(mass.le.zpars(3))then
                  rx = rg
               else
* He-ignition and end of HG occur at Rmin
                  rmin = rminf(mass)
                  ry = ragbf(mt,lums(4),zpars(2))
                  rx = MIN(rmin,ry)
                  if(mass.le.mlp)then
                     texp = log(mass/mlp)/log(zpars(3)/mlp)
                     rx = rg
                     rx = rmin*(rx/rmin)**texp
                  endif
                  tau2 = tblf(mass,zpars(2),zpars(3))
                  if(tau2.lt.tiny) rx = ry
               endif
               r = rtms*(rx/rtms)**tau
               kw = 2
            endif
*
         endif
*
* Now the GB, CHeB and AGB evolution.
*
      elseif(aj.lt.tscls(2))then
*
*        Red Giant.
*
         kw = 3
         lum = lgbtf(aj,GB(1),GB,tscls(4),tscls(5),tscls(6))
         if(mass.le.zpars(2))then
* Star has a degenerate He core which grows on the GB
            mc = mcgbf(lum,GB,lums(6))
         else
* Star has a non-degenerate He core which may grow, but
* only slightly, on the GB
            tau = (aj - tscls(1))/(tscls(2) - tscls(1))
            mcx = mcheif(mass,zpars(2),zpars(9))
            mcy = mcheif(mass,zpars(2),zpars(10))
            mc = mcx + (mcy - mcx)*tau
         endif
         r = rgbf(mt,lum)
         rg = r
         if(mc.ge.mt)then
            aj = 0.d0
            if(mass.gt.zpars(2))then
*
* Zero-age helium star
*
               mc = 0.d0
               mass = mt
               kw = 7
               CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
            else
*
* Zero-age helium white dwarf.
*
               mc = mt
               mass = mt
               kw = 10
            endif
         endif
*
      elseif(aj.lt.tbagb)then
*
*       Core helium burning star.
*
         if(kw.eq.3.and.mass.le.zpars(2))then
            mass = mt
            CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
            aj = tscls(2)
         endif
         if(mass.le.zpars(2))then
            mcx = mcgbf(lums(4),GB,lums(6))
         else
            mcx = mcheif(mass,zpars(2),zpars(10))
         endif
         tau = (aj - tscls(2))/tscls(3)
         mc = mcx + (mcagbf(mass) - mcx)*tau
*
         if(mass.le.zpars(2))then
            lx = lums(5)
            ly = lums(7)
            rx = rzahbf(mt,mc,zpars(2))
            rg = rgbf(mt,lx)
            rmin = rg*zpars(13)**(mass/zpars(2))
            texp = MIN(MAX(0.4d0,rmin/rx),2.5d0)
            ry = ragbf(mt,ly,zpars(2))
            if(rmin.lt.rx)then
               taul = (log(rx/rmin))**(1.d0/3.d0)
            else
               rmin = rx
               taul = 0.d0
            endif
            tauh = (log(ry/rmin))**(1.d0/3.d0)
            tau2 = taul*(tau - 1.d0) + tauh*tau
            r = rmin*exp(abs(tau2)**3)
            rg = rg + tau*(ry - rg)
            lum = lx*(ly/lx)**(tau**texp)
         elseif(mass.gt.zpars(3))then
*
* For HM stars He-ignition takes place at Rmin in the HG, and CHeB
* consists of a blue phase (before tloop) and a RG phase (after tloop).
*
            tau2 = tblf(mass,zpars(2),zpars(3))
            tloop = tscls(2) + tau2*tscls(3)
            rmin = rminf(mass)
            rg = rgbf(mt,lums(4))
            rx = ragbf(mt,lums(4),zpars(2))
            rmin = MIN(rmin, rx)
            if(mass.le.mlp) then
               texp = log(mass/mlp)/log(zpars(3)/mlp)
               rx = rg
               rx = rmin*(rx/rmin)**texp
            else
               rx = rmin
            end if
            texp = MIN(MAX(0.4d0,rmin/rx),2.5d0)
            lum = lums(4)*(lums(7)/lums(4))**(tau**texp)
            if(aj.lt.tloop)then
               ly = lums(4)*(lums(7)/lums(4))**(tau2**texp)
               ry = ragbf(mt,ly,zpars(2))
               taul = 0.d0
               if(ABS(rmin-rx).gt.tiny)then
                  taul = (log(rx/rmin))**(1.d0/3.d0)
               endif
               tauh = 0.d0
               if(ry.gt.rmin) tauh = (log(ry/rmin))**(1.d0/3.d0)
               tau = (aj - tscls(2))/(tau2*tscls(3))
               tau2 = taul*(tau - 1.d0) + tauh*tau
               r = rmin*exp(abs(tau2)**3)
               rg = rg + tau*(ry - rg)
            else
               r = ragbf(mt,lum,zpars(2))
               rg = r
            end if
         else
*
* For IM stars CHeB consists of a RG phase (before tloop) and a blue
* loop (after tloop).
*
            tau2 = 1.d0 - tblf(mass,zpars(2),zpars(3))
            tloop = tscls(2) + tau2*tscls(3)
            if(aj.lt.tloop)then
               tau = (tloop - aj)/(tau2*tscls(3))
               lum = lums(5)*(lums(4)/lums(5))**(tau**3)
               r = rgbf(mt,lum)
               rg = r
            else
               lx = lums(5)
               ly = lums(7)
               rx = rgbf(mt,lx)
               rmin = rminf(mt)
               texp = MIN(MAX(0.4d0,rmin/rx),2.5d0)
               ry = ragbf(mt,ly,zpars(2))
               if(rmin.lt.rx)then
                  taul = (log(rx/rmin))**(1.d0/3.d0)
               else
                  rmin = rx
                  taul = 0.d0
               endif
               tauh = (log(ry/rmin))**(1.d0/3.d0)
               tau = (aj - tloop)/(tscls(3) - (tloop - tscls(2)))
               tau2 = taul*(tau - 1.d0) + tauh*tau
               r = rmin*exp(abs(tau2)**3)
               rg = rx + tau*(ry - rx)
               lum = lx*(ly/lx)**(tau**texp)
            endif
         endif
* 
* Test whether core mass exceeds total mass.
*
         if(mc.ge.mt)then
*
* Evolved MS naked helium star.
*
            kw = 7
            xx = (aj - tscls(2))/tscls(3)
            mass = mt
            CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
            aj = xx*tm
         else
            kw = 4
         endif
*
      else
*
*        Asymptotic Red Giant.
*
* On the AGB the He core mass remains constant until at Ltp it
* is caught by the C core mass and they grow together.
*
         mcbagb = mcagbf(mass)
         mcx = mcgbtf(tbagb,GB(8),GB,tscls(7),tscls(8),tscls(9))
         mcmax = MAX(MAX(mch,0.773d0*mcbagb-0.35d0),1.05d0*mcx)
*
         if(aj.lt.tscls(13))then
            mcx = mcgbtf(aj,GB(8),GB,tscls(7),tscls(8),tscls(9))
            mc = mcbagb
            lum = lmcgbf(mcx,GB)
            if(mt.le.mc)then
*
* Evolved naked helium star as the envelope is lost but the
* star has not completed its interior burning. The star becomes
* a post-HeMS star.
*
               kw = 9
               mt = mc
               mass = mt
               mc = mcx
               CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
               if(mc.le.GB(7))then
                  aj = tscls(4) - (1.d0/((GB(5)-1.d0)*GB(8)*GB(4)))*
     &                            (mc**(1.d0-GB(5)))
               else
                  aj = tscls(5) - (1.d0/((GB(6)-1.d0)*GB(8)*GB(3)))*
     &                            (mc**(1.d0-GB(6)))
               endif
               aj = MAX(aj,tm)
               goto 90
            else
               kw = 5
            endif
         else
            kw = 6
            mc = mcgbtf(aj,GB(2),GB,tscls(10),tscls(11),tscls(12))
            lum = lmcgbf(mc,GB)
*
* Approximate 3rd Dredge-up on AGB by limiting Mc.
*
            lambda = MIN(0.9d0,0.3d0+0.001d0*mass**5)
            tau = tscls(13)
            mcx = mcgbtf(tau,GB(2),GB,tscls(10),tscls(11),tscls(12))
            mcy = mc
            mc = mc - lambda*(mcy-mcx)
            mcx = mc
            mcmax = MIN(mt,mcmax)   
         endif
         r = ragbf(mt,lum,zpars(2))
         rg = r
*
* Mc,x represents the C core mass and we now test whether it
* exceeds either the total mass or the maximum allowed core mass.
*
         if(mcmax-mcx.lt.tiny)then
            aj = 0.d0
            mc = mcmax
            if(mc.lt.mch)then
               if(ifflag.ge.1)then
*
* Invoke WD IFMR from HPE, 1995, MNRAS, 272, 800. 
*
                  if(zpars(14).ge.1.0d-08)then
                     mc = MIN(0.36d0+0.104d0*mass,0.58d0+0.061d0*mass)
                     mc = MAX(0.54d0+0.042d0*mass,mc)
                     if(mass.lt.1.d0) mc = 0.46d0
                  else
                     mc = MIN(0.29d0+0.178d0*mass,0.65d0+0.062d0*mass)
                     mc = MAX(0.54d0+0.073d0*mass,mc)
                  endif
                  mc = MIN(mch,mc)
               endif
*
               mt = mc
               if(mcbagb.lt.1.6d0)then
*     
* Zero-age Carbon/Oxygen White Dwarf
*
                  kw = 11
               else
*     
* Zero-age Oxygen/Neon White Dwarf
*
                  kw = 12
               endif
               mass = mt
*
            else
               if(mcbagb.lt.1.6d0)then
*
* Star is not massive enough to ignite C burning.
* so no remnant is left after the SN
*
                  kw = 15
                  aj = 0.d0
                  mt = 0.d0
                  lum = 1.0d-10
                  r = 1.0d-10
               else
                  if(nsflag.eq.0)then
                     mt = 1.17d0 + 0.09d0*mc
                  elseif(nsflag.ge.1)then
*
* Use NS/BH mass given by Belczynski et al. 2002, ApJ, 572, 407. 
*
                     if(mc.lt.2.5d0)then
                        mcx = 0.161767d0*mc + 1.067055d0
                     else
                        mcx = 0.314154d0*mc + 0.686088d0
                     endif
                     if(mc.le.5.d0)then
                        mt = mcx
                     elseif(mc.lt.7.6d0)then
                        mt = mcx + (mc - 5.d0)*(mt - mcx)/2.6d0
                     endif
                  endif
                  mc = mt
                  if(mt.le.mxns)then
*
* Zero-age Neutron star
*
                     kw = 13
                  else
*
* Zero-age Black hole
*
                     kw = 14
                  endif  
               endif
            endif
         endif
*
      endif
*
 90   continue
*
      if(kw.ge.7.and.kw.le.9)then
*
* Naked Helium Star
*
         rzams = rzhef(mt)
         rx = rzams
         if(aj.lt.tm)then
*
* Main Sequence
*
            kw = 7
            tau = aj/tm
            am = MAX(0.d0,0.85d0-0.08d0*mass)
            lum = lums(1)*(1.d0+0.45d0*tau+am*tau**2)
            am = MAX(0.d0,0.4d0-0.22d0*LOG10(mt))
            r = rx*(1.d0+am*(tau-tau**6))
            rg = rx
* Star has no core mass and hence no memory of its past
* which is why we subject mass and mt to mass loss for
* this phase.
            mc = 0.d0
            if(mt.lt.zpars(10)) kw = 10
         else
*
* Helium Shell Burning
*
            kw = 8
            lum = lgbtf(aj,GB(8),GB,tscls(4),tscls(5),tscls(6))
            r = rhehgf(mt,lum,rx,lums(2))
            rg = rhegbf(lum)
            if(r.ge.rg)then
               kw = 9
               r = rg
            endif
            mc = mcgbf(lum,GB,lums(6))
            mtc = MIN(mt,1.45d0*mt-0.31d0)
            mcmax = MIN(mtc,MAX(mch,0.773d0*mass-0.35d0))
            if(mcmax-mc.lt.tiny)then
               aj = 0.d0
               mc = mcmax
               if(mc.lt.mch)then
                  if(mass.lt.1.6d0)then
*     
* Zero-age Carbon/Oxygen White Dwarf
*
                     mt = MAX(mc,(mc+0.31d0)/1.45d0)
                     kw = 11
                  else
*     
* Zero-age Oxygen/Neon White Dwarf
*
                     mt = mc
                     kw = 12
                  endif
                  mass = mt
               else
                  if(mass.lt.1.6d0)then
*
* Star is not massive enough to ignite C burning.
* so no remnant is left after the SN
*
                     kw = 15
                     aj = 0.d0
                     mt = 0.d0
                     lum = 1.0d-10
                     r = 1.0d-10
                  else
                     if(nsflag.eq.0)then
                        mt = 1.17d0 + 0.09d0*mc
                     elseif(nsflag.ge.1)then
                        if(mc.lt.2.5d0)then
                           mcx = 0.161767d0*mc + 1.067055d0
                        else
                           mcx = 0.314154d0*mc + 0.686088d0
                        endif
                        if(mc.le.5.d0)then
                           mt = mcx
                        elseif(mc.lt.7.6d0)then
                           mt = mcx + (mc - 5.d0)*(mt - mcx)/2.6d0
                        endif
                     endif
                     mc = mt
                     if(mt.le.mxns)then
*
* Zero-age Neutron star
*
                        kw = 13
                     else
*
* Zero-age Black hole
*
                        kw = 14
                     endif
                  endif  
               endif
            endif
         endif
      endif
*
      if(kw.ge.10.and.kw.le.12)then
*
*        White dwarf.
*
         mc = mt
         if(mc.ge.mch)then
*
* Accretion induced supernova with no remnant
* unless WD is ONe in which case we assume a NS 
* of minimum mass is the remnant. 
*
            if(kw.eq.12)then
               kw = 13
               aj = 0.d0
               mt = 1.3d0
            else
               kw = 15
               aj = 0.d0
               mt = 0.d0
               lum = 1.0d-10
               r = 1.0d-10
            endif
         else
*
            if(kw.eq.10)then
               xx = ahe
            else
               xx = aco
            endif
*
            if(wdflag.eq.0)then
*
* Mestel cooling
*
               lum = 635.d0*mt*zpars(14)/(xx*(aj+0.1d0))**1.4d0
*
            elseif(wdflag.ge.1)then
*
* modified-Mestel cooling
*
               if(aj.lt.9000.0)then
                  lum = 300.d0*mt*zpars(14)/(xx*(aj+0.1d0))**1.18d0
               else
                  fac = (9000.1d0*xx)**5.3d0
                  lum = 300.d0*fac*mt*zpars(14)/(xx*(aj+0.1d0))**6.48d0
               endif
*
            endif
*
            r = 0.0115d0*SQRT(MAX(1.48204d-06,(mch/mt)**(2.d0/3.d0)
     &                                      - (mt/mch)**(2.d0/3.d0)))
            r = MIN(0.1d0,r)
            if(mt.lt.0.0005d0) r = 0.09d0
            if(mt.lt.0.000005d0) r = 0.009d0
*
         endif
      endif
*
      if(kw.eq.13)then
*
*        Neutron Star.
*
         mc = mt
         if(mc.gt.mxns)then
*
* Accretion induced Black Hole?
*
            kw = 14
            aj = 0.d0
         else
            lum = 0.02d0*(mt**0.67d0)/(MAX(aj,0.1d0))**2
            r = 1.4d-05
         endif
      endif
*
      if(kw.eq.14)then
*
*        Black hole
*
         mc = mt
         lum = 1.0d-10
         r = 4.24d-06*mt
      endif
*
* Calculate the core radius and the luminosity and radius of the
* remnant that the star will become.
*
      tau = 0.d0
      if(kw.le.1.or.kw.eq.7)then
         rc = 0.d0
      elseif(kw.le.3)then
         if(mass.gt.zpars(2))then
            lx = lzhef(mc)
            rx = rzhef(mc)
            rc = rx
         else
            if(wdflag.eq.0)then
               lx = 635.d0*mc*zpars(14)/((ahe*0.1d0)**1.4d0)
            elseif(wdflag.ge.1)then
               lx = 300.d0*mc*zpars(14)/((ahe*0.1d0)**1.18d0)
            endif
            rx = 0.0115d0*SQRT(MAX(1.48204d-06,
     &           (mch/mc)**(2.d0/3.d0)-(mc/mch)**(2.d0/3.d0)))
            rc = 5.d0*rx
         endif
      elseif(kw.eq.4)then
         tau = (aj - tscls(2))/tscls(3)
         kwp = 7
         CALL star(kwp,mc,mc,tm,tn,tscls,lums,GB,zpars)
         am = MAX(0.d0,0.85d0-0.08d0*mc)
         lx = lums(1)*(1.d0+0.45d0*tau+am*tau**2)
         rx = rzhef(mc)
         am = MAX(0.d0,0.4d0-0.22d0*LOG10(mc))
         rx = rx*(1.d0+am*(tau-tau**6))
         CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
         rc = rx
      elseif(kw.eq.5)then
         kwp = 9
         if(tn.gt.tbagb) tau = 3.d0*(aj-tbagb)/(tn-tbagb)
         CALL star(kwp,mc,mc,tm,tn,tscls,lums,GB,zpars)
         lx = lmcgbf(mcx,GB)
         if(tau.lt.1.d0) lx = lums(2)*(lx/lums(2))**tau
         rx = rzhef(mc)
         rx = MIN(rhehgf(mc,lx,rx,lums(2)),rhegbf(lx))
         CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
         rc = rx
      elseif(kw.le.9)then
         if(wdflag.eq.0)then
            lx = 635.d0*mc*zpars(14)/((aco*0.1d0)**1.4d0)
         elseif(wdflag.ge.1)then
            lx = 300.d0*mc*zpars(14)/((aco*0.1d0)**1.18d0)
         endif
         rx = 0.0115d0*SQRT(MAX(1.48204d-06,
     &        (mch/mc)**(2.d0/3.d0) - (mc/mch)**(2.d0/3.d0)))
         rc = 5.d0*rx
      else
         rc = r
         menv = 1.0d-10
         renv = 1.0d-10
         k2 = 0.21d0
      endif
*
* Perturb the luminosity and radius due to small envelope mass.
*
      if(kw.ge.2.and.kw.le.9.and.kw.ne.7)then
         mew = ((mt-mc)/mt)*MIN(5.d0,MAX(1.2d0,(lum/lum0)**kap))
         if(kw.ge.8) mew = ((mtc-mc)/mtc)*5.d0
         if(mew.lt.1.d0)then
            xx = lpertf(mt,mew)
            lum = lx*(lum/lx)**xx
            if(r.le.rx)then
               xx = 0.d0
            else
               xx = rpertf(mt,mew,r,rx)
            endif
            r = rx*(r/rx)**xx
         endif
         rc = MIN(rc,r)
      endif
*
* Calculate mass and radius of convective envelope, and envelope
* gyration radius.
*
      if(kw.lt.10)then
         CALL mrenv(kw,mass,mt,mc,lum,r,rc,aj,tm,lums(2),lums(3),
     &              lums(4),rzams,rtms,rg,menv,renv,k2)
      endif
*
      if(mass.gt.99.99d0)then
         mass = mass0
      endif
      if(mt.gt.99.99d0)then
         mt = mt0
      endif
*
      return
      end
***

***
      SUBROUTINE kick(kw,m1,m1n,m2,ecc,sep,jorb,vs)
      implicit none
*
      integer kw,k
      INTEGER idum
      COMMON /VALUE3/ idum
      INTEGER idum2,iy,ir(32)
      COMMON /RAND3/ idum2,iy,ir
      integer bhflag
      real*8 m1,m2,m1n,ecc,sep,jorb,ecc2
      real*8 pi,twopi,gmrkm,yearsc,rsunkm
      parameter(yearsc=3.1557d+07,rsunkm=6.96d+05)
      real*8 mm,em,dif,der,del,r
      real*8 u1,u2,vk,v(4),s,theta,phi
      real*8 sphi,cphi,stheta,ctheta,salpha,calpha
      real*8 vr,vr2,vk2,vn2,hn2
      real*8 mu,cmu,vs(3),v1,v2,mx1,mx2
      real*8 sigma
      COMMON /VALUE4/ sigma,bhflag
      real ran3,xx
      external ran3
*
      do k = 1,3
         vs(k) = 0.d0
      enddo
*     if(kw.eq.14.and.bhflag.eq.0) goto 95
*
      pi = ACOS(-1.d0)
      twopi = 2.d0*pi
* Conversion factor to ensure velocities are in km/s using mass and
* radius in solar units.
      gmrkm = 1.906125d+05
*
* Find the initial separation by randomly choosing a mean anomaly.
      if(sep.gt.0.d0.and.ecc.ge.0.d0)then
         xx = RAN3(idum)
         mm = xx*twopi
         em = mm
 2       dif = em - ecc*SIN(em) - mm
         if(ABS(dif/mm).le.1.0d-04) goto 3
         der = 1.d0 - ecc*COS(em)
         del = dif/der
         em = em - del
         goto 2
 3       continue
         r = sep*(1.d0 - ecc*COS(em))
*
* Find the initial relative velocity vector.
         salpha = SQRT((sep*sep*(1.d0-ecc*ecc))/(r*(2.d0*sep-r)))
         calpha = (-1.d0*ecc*SIN(em))/SQRT(1.d0-ecc*ecc*COS(em)*COS(em))
         vr2 = gmrkm*(m1+m2)*(2.d0/r - 1.d0/sep)
         vr = SQRT(vr2)
      else
         vr = 0.d0
         vr2 = 0.d0
         salpha = 0.d0
         calpha = 0.d0
      endif
*
* Generate Kick Velocity using Maxwellian Distribution (Phinney 1992).
* Use Henon's method for pairwise components (Douglas Heggie 22/5/97).
      do 20 k = 1,2
         u1 = RAN3(idum)
         u2 = RAN3(idum)
* Generate two velocities from polar coordinates S & THETA.
         s = sigma*SQRT(-2.d0*LOG(1.d0 - u1))
         theta = twopi*u2
         v(2*k-1) = s*COS(theta)
         v(2*k) = s*SIN(theta)
 20   continue
      vk2 = v(1)**2 + v(2)**2 + v(3)**2
      vk = SQRT(vk2)
      if((kw.eq.14.and.bhflag.eq.0).or.kw.lt.0)then
         vk2 = 0.d0
         vk = 0.d0
         if(kw.lt.0) kw = 13
      endif
      sphi = -1.d0 + 2.d0*u1
      phi = ASIN(sphi)
      cphi = COS(phi)
      stheta = SIN(theta)
      ctheta = COS(theta)
*     WRITE(66,*)' KICK VK PHI THETA ',vk,phi,theta
      if(sep.le.0.d0.or.ecc.lt.0.d0) goto 90
*
* Determine the magnitude of the new relative velocity.
      vn2 = vk2+vr2-2.d0*vk*vr*(ctheta*cphi*salpha-stheta*cphi*calpha)
* Calculate the new semi-major axis.
      sep = 2.d0/r - vn2/(gmrkm*(m1n+m2))
      sep = 1.d0/sep
*     if(sep.le.0.d0)then
*        ecc = 1.1d0
*        goto 90
*     endif
* Determine the magnitude of the cross product of the separation vector
* and the new relative velocity.
      v1 = vk2*sphi*sphi
      v2 = (vk*ctheta*cphi-vr*salpha)**2
      hn2 = r*r*(v1 + v2)
* Calculate the new eccentricity.
      ecc2 = 1.d0 - hn2/(gmrkm*sep*(m1n+m2))
      ecc2 = MAX(ecc2,0.d0)
      ecc = SQRT(ecc2)
* Calculate the new orbital angular momentum taking care to convert
* hn to units of Rsun^2/yr.
      jorb = (m1n*m2/(m1n+m2))*SQRT(hn2)*(yearsc/rsunkm)
* Determine the angle between the new and old orbital angular
* momentum vectors.
      cmu = (vr*salpha-vk*ctheta*cphi)/SQRT(v1 + v2)
      mu = ACOS(cmu)
* Calculate the components of the velocity of the new centre-of-mass.
 90   continue
      if(ecc.le.1.0)then
* Calculate the components of the velocity of the new centre-of-mass.
         mx1 = vk*m1n/(m1n+m2)
         mx2 = vr*(m1-m1n)*m2/((m1n+m2)*(m1+m2))
         vs(1) = mx1*ctheta*cphi + mx2*salpha
         vs(2) = mx1*stheta*cphi + mx2*calpha
         vs(3) = mx1*sphi
      else
* Calculate the relative hyperbolic velocity at infinity (simple method).
         sep = r/(ecc-1.d0)
*        cmu = SQRT(ecc-1.d0)
*        mu = ATAN(cmu)
         mu = ACOS(1.d0/ecc)
         vr2 = gmrkm*(m1n+m2)/sep
         vr = SQRT(vr2)
         vs(1) = vr*SIN(mu)
         vs(2) = vr*COS(mu)
         vs(3) = 0.d0
         ecc = MIN(ecc,99.99d0)
      endif
*
 95   continue
*
      RETURN
      END
***

***
      real*8 FUNCTION mlwind(kw,lum,r,mt,mc,rl,z)
      implicit none
      integer kw
      real*8 lum,r,mt,mc,rl,z
      real*8 dml,dms,dmt,p0,x,mew,lum0,kap,neta,bwind,hewind,mxns
      parameter(lum0=7.0d+04,kap=-0.5d0)
      common /value1/ neta,bwind,hewind,mxns
*
* Calculate stellar wind mass loss.
*
* Apply mass loss of Nieuwenhuijzen & de Jager, A&A, 1990, 231, 134,
* for massive stars over the entire HRD.
      dms = 0.d0
      if(lum.gt.4000.d0)then
         x = MIN(1.d0,(lum-4000.d0)/500.d0)
         dms = 9.6d-15*x*(r**0.81d0)*(lum**1.24d0)*(mt**0.16d0)
         dms = dms*(z/0.02d0)**(1.d0/2.d0)
      endif
      if(kw.ge.2.and.kw.le.9)then
* 'Reimers' mass loss
         dml = neta*4.0d-13*r*lum/mt
         if(rl.gt.0.d0) dml = dml*(1.d0 + bwind*(MIN(0.5d0,(r/rl)))**6)
* Apply mass loss of Vassiliadis & Wood, ApJ, 1993, 413, 641, 
* for high pulsation periods on AGB.
         if(kw.eq.5.or.kw.eq.6)then
            p0 = -2.07d0 - 0.9d0*log10(mt) + 1.94d0*log10(r)
            p0 = 10.d0**p0
            p0 = MIN(p0,2000.d0)
            dmt = -11.4d0+0.0125d0*(p0-100.d0*MAX(mt-2.5d0,0.d0))
            dmt = 10.d0**dmt
            dmt = 1.d0*MIN(dmt,1.36d-09*lum)
            dml = MAX(dml,dmt)
         endif
         if(kw.gt.6)then
            dms = MAX(dml,1.0d-13*hewind*lum**(3.d0/2.d0))
         else
            dms = MAX(dml,dms)
            mew = ((mt-mc)/mt)*MIN(5.d0,MAX(1.2d0,(lum/lum0)**kap))
* reduced WR-like mass loss for small H-envelope mass
            if(mew.lt.1.d0)then
               dml = 1.0d-13*lum**(3.d0/2.d0)*(1.d0 - mew)
               dms = MAX(dml,dms)
            end if
* LBV-like mass loss beyond the Humphreys-Davidson limit.
            x = 1.0d-5*r*sqrt(lum)
            if(lum.gt.6.0d+05.and.x.gt.1.d0)then
               dml = 0.1d0*(x-1.d0)**3*(lum/6.0d+05-1.d0)
               dms = dms + dml
            endif
         endif
      endif
*
      mlwind = dms
*
      return
      end
***

***
      SUBROUTINE mrenv(kw,mass,mt,mc,lum,rad,rc,aj,tm,ltms,lbgb,lhei,
     &                 rzams,rtms,rg,menv,renv,k2e)
      implicit none
      integer kw
      real*8 mass,mt,mc,lum,rad,rc,aj,tm
      real*8 k2e,menv,menvg,menvt,menvz,renv,renvg,renvt,renvz
      real*8 A,B,C,D,E,F,x,y
      real*8 k2bgb,k2g,k2z,logm,logmt,lbgb,ltms,lhei,rg,rtms,rzams
      real*8 teff,tebgb,tetms,tau,tauenv,tautms
*
* A function to estimate the mass and radius of the convective envelope,
* as well as the gyration radius of the envelope.
* N.B. Valid only for Z=0.02!
*
* The following input is needed from HRDIAG:
*   kw = stellar type
*   mass = zero-age stellar mass
*   mt = actual mass
*   mc = core mass (not really needed, can also be done outside subroutine)
*   lum = luminosity
*   rad = radius
*   rc = core radius (not really needed...)
*   aj = age
*   tm = main-sequence lifetime
*   ltms = luminosity at TMS, lums(2)
*   lbgb = luminosity at BGB, lums(3)
*   lhei = luminosity at He ignition, lums(4)
*   rzams = radius at ZAMS
*   rtms = radius at TMS
*   rg = giant branch or Hayashi track radius, approporaite for the type. 
*        For kw=1 or 2 this is radius at BGB, and for kw=4 either GB or 
*        AGB radius at present luminosity.
*
      logm = log10(mass)
      A = MIN(0.81d0,MAX(0.68d0,0.68d0+0.4d0*logm))
      C = MAX(-2.5d0,MIN(-1.5d0,-2.5d0+5.d0*logm))
      D = -0.1d0
      E = 0.025d0
*
* Zero-age and BGB values of k^2.
*
      k2z = MIN(0.21d0,MAX(0.09d0-0.27d0*logm,0.037d0+0.033d0*logm))
      if(logm.gt.1.3d0) k2z = k2z - 0.055d0*(logm-1.3d0)**2
      k2bgb = MIN(0.15d0,MIN(0.147d0+0.03d0*logm,0.162d0-0.04d0*logm))
*
      if(kw.ge.3.and.kw.le.6)then
*
* Envelope k^2 for giant-like stars; this will be modified for non-giant
* CHeB stars or small envelope mass below.
* Formula is fairly accurate for both FGB and AGB stars if M <= 10, and
* gives reasonable values for higher masses. Mass dependence is on actual
* rather than ZA mass, expected to work for mass-losing stars (but not
* tested!). The slightly complex appearance is to insure continuity at 
* the BGB, which depends on the ZA mass.
*
         logmt = log10(mt)
         F = 0.208d0 + 0.125d0*logmt - 0.035d0*logmt**2
         B = 1.0d+04*mt**(3.d0/2.d0)/(1.d0+0.1d0*mt**(3.d0/2.d0))
         x = ((lum-lbgb)/B)**2
         y = (F - 0.033d0*log10(lbgb))/k2bgb - 1.d0
         k2g = (F - 0.033d0*log10(lum) + 0.4d0*x)/(1.d0+y*(lbgb/lum)+x)
      elseif(kw.eq.9)then
*
* Rough fit for for HeGB stars...
*
         B = 3.0d+04*mt**(3.d0/2.d0)
         x = (MAX(0.d0,lum/B-0.5d0))**2
         k2g = (k2bgb + 0.4d0*x)/(1.d0 + 0.4d0*x)
      else
         k2g = k2bgb
      endif
*
      if(kw.le.2)then
         menvg = 0.5d0
         renvg = 0.65d0
      elseif(kw.eq.3.and.lum.lt.3.d0*lbgb)then
*
* FGB stars still close to the BGB do not yet have a fully developed CE.
*
         x = MIN(3.d0,lhei/lbgb)
         tau = MAX(0.d0,MIN(1.d0,(x-lum/lbgb)/(x-1.d0)))
         menvg = 1.d0 - 0.5d0*tau**2
         renvg = 1.d0 - 0.35d0*tau**2
      else
         menvg = 1.d0
         renvg = 1.d0
      endif
*
      if(rad.lt.rg)then
*
* Stars not on the Hayashi track: MS and HG stars, non-giant CHeB stars,
* HeMS and HeHG stars, as well as giants with very small envelope mass.
*
         
         if(kw.le.6)then
*
* Envelope k^2 fitted for MS and HG stars.
* Again, pretty accurate for M <= 10 but less so for larger masses.
* [Note that this represents the whole star on the MS, so there is a 
* discontinuity in stellar k^2 between MS and HG - okay for stars with a 
* MS hook but low-mass stars should preferably be continous...]
*
* For other types of star not on the Hayashi track we use the same fit as 
* for HG stars, this is not very accurate but has the correct qualitative 
* behaviour. For CheB stars this is an overestimate because they appear
* to have a more centrally concentrated envelope than HG stars.
*
            k2e = (k2z-E)*(rad/rzams)**C + E*(rad/rzams)**D
         elseif(kw.eq.7)then
* Rough fit for naked He MS stars.
            tau = aj/tm
            k2e = 0.08d0 - 0.03d0*tau
         elseif(kw.le.9)then
* Rough fit for HeHG stars.
            k2e = 0.08d0*rzams/rad
         endif
*
* tauenv measures proximity to the Hayashi track in terms of Teff.
* If tauenv>0 then an appreciable convective envelope is present, and
* k^2 needs to be modified.
*
         if(kw.le.2)then
            teff = sqrt(sqrt(lum)/rad)
            tebgb = sqrt(sqrt(lbgb)/rg)
            tauenv = MAX(0.d0,MIN(1.d0,(tebgb/teff-A)/(1.d0-A)))
         else
            tauenv = MAX(0.d0,MIN(1.d0,(sqrt(rad/rg)-A)/(1.d0-A)))
         endif
*
         if(tauenv.gt.0.d0)then
            menv = menvg*tauenv**5
            renv = renvg*tauenv**(5.d0/4.d0)
            if(kw.le.1)then
* Zero-age values for CE mass and radius.
               x = MAX(0.d0,MIN(1.d0,(0.1d0-logm)/0.55d0))
               menvz = 0.18d0*x + 0.82d0*x**5
               renvz = 0.4d0*x**(1.d0/4.d0) + 0.6d0*x**10
               y = 2.d0 + 8.d0*x
* Values for CE mass and radius at start of the HG.
               tetms = sqrt(sqrt(ltms)/rtms)
               tautms = MAX(0.d0,MIN(1.d0,(tebgb/tetms-A)/(1.d0-A)))
               menvt = menvg*tautms**5
               renvt = renvg*tautms**(5.d0/4.d0)
* Modified expressions during MS evolution.
               tau = aj/tm
               if(tautms.gt.0.d0)then
                  menv = menvz + tau**y*menv*(menvt - menvz)/menvt
                  renv = renvz + tau**y*renv*(renvt - renvz)/renvt
               else
                  menv = 0.d0
                  renv = 0.d0
               endif
               k2e = k2e + tau**y*tauenv**3*(k2g - k2e)
            else
               k2e = k2e + tauenv**3*(k2g - k2e)
            endif
         else
            menv = 0.d0
            renv = 0.d0
         endif
      else
*
* All other stars should be true giants.
*
         menv = menvg
         renv = renvg
         k2e = k2g
      endif
*
      menv = menv*(mt - mc)
      renv = renv*(rad - rc)
      menv = MAX(menv,1.0d-10)
      renv = MAX(renv,1.0d-10)
*
      return
      end
***

***
      REAL FUNCTION ran3(IDUM)
*
* Random number generator from Numerical Recipes, Press et al. pg 272.
*
      IMPLICIT NONE
      INTEGER j,k,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      PARAMETER(im1=2147483563,im2=2147483399,ia1=40014,ia2=40692)
      PARAMETER(iq1=53668,iq2=52774,ir1=12211,ir2=3791,ntab=32)
      INTEGER idum
      INTEGER idum2,iy,ir(ntab)
      REAL am
      COMMON /RAND3/ idum2,iy,ir
      DATA idum2/123456789/, iy/0/, ir/ntab*0/
*
      am = 1.0/float(im1)
      imm1 = im1 - 1
      ndiv = 1 + imm1/ntab
*
      if(idum.le.0)then
         idum = MAX(-idum,1)
         idum2 = idum
         do 11 , j = ntab+8,1,-1
            k = idum/iq1
            idum = ia1*(idum-k*iq1)-k*ir1
            if(idum.lt.0) idum = idum + im1
            if(j.le.ntab) ir(j) = idum
 11      continue
         iy = ir(1)
      endif
      k = idum/iq1
      idum = ia1*(idum-k*iq1)-k*ir1
      if(idum.lt.0) idum = idum + im1
      k = idum2/iq2
      idum2 = ia2*(idum2-k*iq2)-k*ir2
      if(idum2.lt.0) idum2 = idum2 + im2
      j = 1 + iy/ndiv
      iy = ir(j) - idum2
      ir(j) = idum
      if(iy.lt.1) iy = iy + imm1
      ran3 = am*iy
*
      RETURN
      END
***

***
      SUBROUTINE star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
*
*
*       Stellar luminosity & evolution time. 
*       ------------------------------------
*
      implicit none
*
      integer kw
*
      real*8 mass,mt,tm,tn,tscls(20),lums(10),GB(10),zpars(20)
      real*8 tgb,tbagb,mch,mcmax,mc1,mc2,mcbagb,dx,am
      real*8 lambda,tau,mtc,mass0
      parameter(mch=1.44d0)
*
      real*8 lzamsf,lzahbf,lzhef
      real*8 tbgbf,thookf,tHef,themsf,mcgbf,mcagbf,mcheif,mcgbtf
      real*8 ltmsf,lbgbf,lHeIf,lHef,lbagbf,lmcgbf
      external lzamsf,lzahbf,lzhef
      external tbgbf,thookf,tHef,themsf,mcgbf,mcagbf,mcheif,mcgbtf
      external ltmsf,lbgbf,lHeIf,lHef,lbagbf,lmcgbf
*
*       Computes the characteristic luminosities at different stages (LUMS),
*       and various timescales (TSCLS).
*       Ref: P.P. Eggleton, M.J. Fitchett & C.A. Tout (1989) Ap.J. 347, 998.
*
*       Revised 27th March 1995 by C. A. Tout
*       and 24th October 1995 to include metallicity
*       and 13th December 1996 to include naked helium stars
*
*       Revised 5th April 1997 by J. R. Hurley
*       to include Z=0.001 as well as Z=0.02, convective overshooting,
*       MS hook and more elaborate CHeB. It now also sets the Giant
*       Branch parameters relevant to the mass of the star.
*
*       ------------------------------------------------------------
*       Times: 1; BGB              2; He ignition   3; He burning
*              4; Giant t(inf1)    5; Giant t(inf2) 6; Giant t(Mx)
*              7; FAGB t(inf1)     8; FAGB t(inf2)  9; FAGB  t(Mx)
*             10; SAGB t(inf1)    11; SAGB t(inf2) 12; SAGB  t(Mx)
*             13; TP              14; t(Mcmax)     
*
*       LUMS:  1; ZAMS             2; End MS        3; BGB
*              4; He ignition      5; He burning    6; L(Mx)
*              7; BAGB             8; TP
*
*       GB:    1; effective A(H)   2; A(H,He)       3; B
*              4; D                5; p             6; q
*              7; Mx               8; A(He)         9; Mc,BGB
*
*       ------------------------------------------------------------
*
*
      mass0 = mass
      if(mass0.gt.100.d0) mass = 100.d0
*
      if(kw.ge.7.and.kw.le.9) goto 90
      if(kw.ge.10) goto 95
*
* MS and BGB times
*
      tscls(1) = tbgbf(mass)
      tm = MAX(zpars(8),thookf(mass))*tscls(1)
*
* Zero- and terminal age main sequence luminosity
*
      lums(1) = lzamsf(mass)
      lums(2) = ltmsf(mass)
*
* Set the GB parameters
*
      GB(1) = MAX(-4.8d0,MIN(-5.7d0+0.8d0*mass,-4.1d0+0.14d0*mass))
      GB(1) = 10.d0**GB(1)
      GB(2) = 1.27d-05
      GB(8) = 8.0d-05
      GB(3) = MAX(3.0d+04,500.d0 + 1.75d+04*mass**0.6d0)
      if(mass.le.2.0)then
         GB(4) = zpars(6)
         GB(5) = 6.d0
         GB(6) = 3.d0
      elseif(mass.lt.2.5)then
         dx = zpars(6) - (0.975d0*zpars(6) - 0.18d0*2.5d0)
         GB(4) = zpars(6) - dx*(mass - 2.d0)/(0.5d0)
         GB(5) = 6.d0 - (mass - 2.d0)/(0.5d0)
         GB(6) = 3.d0 - (mass - 2.d0)/(0.5d0)
      else
         GB(4) = MAX(-1.d0,0.5d0*zpars(6) - 0.06d0*mass)
         GB(4) = MAX(GB(4),0.975d0*zpars(6) - 0.18d0*mass)
         GB(5) = 5.d0
         GB(6) = 2.d0
      endif
      GB(4) = 10.d0**GB(4)
      GB(7) = (GB(3)/GB(4))**(1.d0/(GB(5)-GB(6)))
*
* Change in slope of giant L-Mc relation.
      lums(6) = GB(4)*GB(7)**GB(5)
*
* HeI ignition luminosity
      lums(4) = lHeIf(mass,zpars(2)) 
      lums(7) = lbagbf(mass,zpars(2))
*
      if(mass.lt.0.1d0.and.kw.le.1)then
         tscls(2) = 1.1d0*tscls(1)
         tscls(3) = 0.1d0*tscls(1)
         lums(3) = lbgbf(mass) 
         goto 96
      endif
*
      if(mass.le.zpars(3))then
* Base of the giant branch luminosity
         lums(3) = lbgbf(mass) 
* Set GB timescales 
         tscls(4) = tscls(1) + (1.d0/((GB(5)-1.d0)*GB(1)*GB(4)))*
     &              ((GB(4)/lums(3))**((GB(5)-1.d0)/GB(5)))
         tscls(6) = tscls(4) - (tscls(4) - tscls(1))*((lums(3)/lums(6))
     &              **((GB(5)-1.d0)/GB(5)))
         tscls(5) = tscls(6) + (1.d0/((GB(6)-1.d0)*GB(1)*GB(3)))*
     &              ((GB(3)/lums(6))**((GB(6)-1.d0)/GB(6)))
* Set Helium ignition time
         if(lums(4).le.lums(6))then
            tscls(2) = tscls(4) - (1.d0/((GB(5)-1.d0)*GB(1)*GB(4)))*
     &                      ((GB(4)/lums(4))**((GB(5)-1.d0)/GB(5)))
         else
            tscls(2) = tscls(5) - (1.d0/((GB(6)-1.d0)*GB(1)*GB(3)))*
     &                      ((GB(3)/lums(4))**((GB(6)-1.d0)/GB(6)))
         endif
         tgb = tscls(2) - tscls(1)
         if(mass.le.zpars(2))then
            mc1 = mcgbf(lums(4),GB,lums(6))
            mc2 = mcagbf(mass)
            lums(5) = lzahbf(mass,mc1,zpars(2))
            tscls(3) = tHef(mass,mc1,zpars(2))
         else
            lums(5) = lHef(mass)*lums(4)
            tscls(3) = tHef(mass,1.d0,zpars(2))*tscls(1)
         endif
      else
* Note that for M>zpars(3) there is no GB as the star goes from
* HG -> CHeB -> AGB. So in effect tscls(1) refers to the time of
* Helium ignition and not the BGB.
         tscls(2) = tscls(1)
         tscls(3) = tHef(mass,1.d0,zpars(2))*tscls(1)
* This now represents the luminosity at the end of CHeB, ie. BAGB
         lums(5) = lums(7)
* We set lums(3) to be the luminosity at the end of the HG
         lums(3) = lums(4)
      endif
*
* Set the core mass at the BGB.
*
      if(mass.le.zpars(2))then
         GB(9) = mcgbf(lums(3),GB,lums(6))
      elseif(mass.le.zpars(3))then
         GB(9) = mcheif(mass,zpars(2),zpars(9))
      else
         GB(9) = mcheif(mass,zpars(2),zpars(10))
      endif
*
* FAGB time parameters
*
      tbagb = tscls(2) + tscls(3)
      tscls(7) = tbagb + (1.d0/((GB(5)-1.d0)*GB(8)*GB(4)))*
     &               ((GB(4)/lums(7))**((GB(5)-1.d0)/GB(5)))
      tscls(9) = tscls(7) - (tscls(7) - tbagb)*((lums(7)/lums(6))
     &                                    **((GB(5)-1.d0)/GB(5)))
      tscls(8) = tscls(9) + (1.d0/((GB(6)-1.d0)*GB(8)*GB(3)))*
     &               ((GB(3)/lums(6))**((GB(6)-1.d0)/GB(6)))
*
* Now to find Ltp and ttp using Mc,He,tp
*
      mcbagb = mcagbf(mass)
      mc1 = mcbagb
      if(mc1.ge.0.8d0.and.mc1.lt.2.25d0)then
* The star undergoes dredge-up at Ltp causing a decrease in Mc,He
         mc1 = 0.44d0*mc1 + 0.448d0
      endif
      lums(8) = lmcgbf(mc1,GB)
      if(mc1.le.GB(7))then
         tscls(13) = tscls(7) - (1.d0/((GB(5)-1.d0)*GB(8)*GB(4)))*
     &                   (mc1**(1.d0-GB(5)))
      else
         tscls(13) = tscls(8) - (1.d0/((GB(6)-1.d0)*GB(8)*GB(3)))*
     &                   (mc1**(1.d0-GB(6)))
      endif
*
* SAGB time parameters
*
      if(mc1.le.GB(7))then
         tscls(10) = tscls(13) + (1.d0/((GB(5)-1.d0)*GB(2)*GB(4)))*
     &               ((GB(4)/lums(8))**((GB(5)-1.d0)/GB(5)))
         tscls(12) = tscls(10) - (tscls(10) - tscls(13))*
     &               ((lums(8)/lums(6))**((GB(5)-1.d0)/GB(5)))
         tscls(11) = tscls(12) + (1.d0/((GB(6)-1.d0)*GB(2)*GB(3)))*
     &               ((GB(3)/lums(6))**((GB(6)-1.d0)/GB(6)))
      else
         tscls(10) = tscls(7)
         tscls(12) = tscls(9)
         tscls(11) = tscls(13) + (1.d0/((GB(6)-1.d0)*GB(2)*GB(3)))*
     &               ((GB(3)/lums(8))**((GB(6)-1.d0)/GB(6)))
      endif
*
* Get an idea of when Mc,C = Mc,C,max on the AGB
      tau = tscls(2) + tscls(3)
      mc2 = mcgbtf(tau,GB(8),GB,tscls(7),tscls(8),tscls(9))
      mcmax = MAX(MAX(mch,0.773d0*mcbagb - 0.35d0),1.05d0*mc2)
*
      if(mcmax.le.mc1)then
         if(mcmax.le.GB(7))then
            tscls(14) = tscls(7) - (1.d0/((GB(5)-1.d0)*GB(8)*GB(4)))*
     &                      (mcmax**(1.d0-GB(5)))
         else
            tscls(14) = tscls(8) - (1.d0/((GB(6)-1.d0)*GB(8)*GB(3)))*
     &                      (mcmax**(1.d0-GB(6)))
         endif
      else
* Star is on SAGB and we need to increase mcmax if any 3rd
* dredge-up has occurred.
         lambda = MIN(0.9d0,0.3d0+0.001d0*mass**5)
         mcmax = (mcmax - lambda*mc1)/(1.d0 - lambda)
         if(mcmax.le.GB(7))then
            tscls(14) = tscls(10) - (1.d0/((GB(5)-1.d0)*GB(2)*GB(4)))*
     &                      (mcmax**(1.d0-GB(5)))
         else
            tscls(14) = tscls(11) - (1.d0/((GB(6)-1.d0)*GB(2)*GB(3)))*
     &                      (mcmax**(1.d0-GB(6)))
         endif
      endif
      tscls(14) = MAX(tbagb,tscls(14))
      if(mass.ge.100.d0)then
         tn = tscls(2)
         goto 100
      endif
*
* Calculate the nuclear timescale - the time of exhausting
* nuclear fuel without further mass loss.
* This means we want to find when Mc = Mt which defines Tn and will
* be used in determining the timestep required. Note that after some 
* stars reach Mc = Mt there will be a Naked Helium Star lifetime
* which is also a nuclear burning period but is not included in Tn.
*
      if(ABS(mt-mcbagb).lt.1.0d-14.and.kw.lt.5)then
         tn = tbagb
      else
* Note that the only occurence of Mc being double-valued is for stars
* that have a dredge-up. If Mt = Mc where Mc could be the value taken
* from CHeB or from the AGB we need to check the current stellar type.
         if(mt.gt.mcbagb.or.(mt.ge.mc1.and.kw.gt.4))then
            if(kw.eq.6)then
               lambda = MIN(0.9d0,0.3d0+0.001d0*mass**5)
               mc1 = (mt - lambda*mc1)/(1.d0 - lambda)
            else
               mc1 = mt
            endif
            if(mc1.le.GB(7))then
               tn = tscls(10) - (1.d0/((GB(5)-1.d0)*GB(2)*GB(4)))*
     &                         (mc1**(1.d0-GB(5)))
            else
               tn = tscls(11) - (1.d0/((GB(6)-1.d0)*GB(2)*GB(3)))*
     &                         (mc1**(1.d0-GB(6)))
            endif
         else
            if(mass.gt.zpars(3))then
               mc1 = mcheif(mass,zpars(2),zpars(10))
               if(mt.le.mc1)then
                  tn = tscls(2)
               else
                  tn = tscls(2) + tscls(3)*((mt - mc1)/(mcbagb - mc1))
               endif
            elseif(mass.le.zpars(2))then
               mc1 = mcgbf(lums(3),GB,lums(6))
               mc2 = mcgbf(lums(4),GB,lums(6))
               if(mt.le.mc1)then
                  tn = tscls(1)
               elseif(mt.le.mc2)then
                  if(mt.le.GB(7))then
                     tn = tscls(4) - (1.d0/((GB(5)-1.d0)*GB(1)*GB(4)))*
     &                               (mt**(1.d0-GB(5)))
                  else
                     tn = tscls(5) - (1.d0/((GB(6)-1.d0)*GB(1)*GB(3)))*
     &                               (mt**(1.d0-GB(6)))
                  endif
               else
                  tn = tscls(2) + tscls(3)*((mt - mc2)/(mcbagb - mc2))
               endif
            else
               mc1 = mcheif(mass,zpars(2),zpars(9))
               mc2 = mcheif(mass,zpars(2),zpars(10))
               if(mt.le.mc1)then
                  tn = tscls(1)
               elseif(mt.le.mc2)then
                  tn = tscls(1) + tgb*((mt - mc1)/(mc2 - mc1))
               else
                  tn = tscls(2) + tscls(3)*((mt - mc2)/(mcbagb - mc2))
               endif
            endif
         endif
      endif
      tn = MIN(tn,tscls(14))
*
      goto 100
*
 90   continue
*
* Calculate Helium star Main Sequence lifetime.
*
      tm = themsf(mass)
      tscls(1) = tm
*
* Zero- and terminal age Helium star main sequence luminosity
*
      lums(1) = lzhef(mass)
      am = MAX(0.d0,0.85d0-0.08d0*mass)
      lums(2) = lums(1)*(1.d0+0.45d0+am)
*
* Set the Helium star GB parameters
*
      GB(8) = 8.0d-05
      GB(3) = 4.1d+04
      GB(4) = 5.5d+04/(1.d0+0.4d0*mass**4)
      GB(5) = 5.d0
      GB(6) = 3.d0
      GB(7) = (GB(3)/GB(4))**(1.d0/(GB(5)-GB(6)))
* Change in slope of giant L-Mc relation.
      lums(6) = GB(4)*GB(7)**GB(5)
*
*** Set Helium star GB timescales 
*
      mc1 = mcgbf(lums(2),GB,lums(6))
      tscls(4) = tm + (1.d0/((GB(5)-1.d0)*GB(8)*GB(4)))*
     &                                          mc1**(1.d0-GB(5))
      tscls(6) = tscls(4) - (tscls(4) - tm)*((GB(7)/mc1)
     &                                     **(1.d0-GB(5)))
      tscls(5) = tscls(6) + (1.d0/((GB(6)-1.d0)*GB(8)*GB(3)))*
     &                                       GB(7)**(1.d0-GB(6))
*
* Get an idea of when Mc = MIN(Mt,Mc,C,max) on the GB
      mtc = MIN(mt,1.45d0*mt-0.31d0)
      if(mtc.le.0.d0) mtc = mt
      mcmax = MIN(mtc,MAX(mch,0.773d0*mass-0.35d0))
      if(mcmax.le.GB(7))then
         tscls(14) = tscls(4) - (1.d0/((GB(5)-1.d0)*GB(8)*GB(4)))*
     &                   (mcmax**(1.d0-GB(5)))
      else
         tscls(14) = tscls(5) - (1.d0/((GB(6)-1.d0)*GB(8)*GB(3)))*
     &                   (mcmax**(1.d0-GB(6)))
      endif
      tscls(14) = MAX(tscls(14),tm)
      tn = tscls(14)
*
      goto 100
*
 95   continue
      tm = 1.0d+10
      tscls(1) = tm
 96   continue
      tn = 1.0d+10
*
 100  continue
      mass = mass0
*
      return
      end
***

***
      SUBROUTINE zcnsts(z,zpars)
* 
      implicit none
      integer kw
*
      real*8 z,zpars(20)
      real*8 tm,tn,tscls(20),lums(10),GB(10)
      real*8 lzs,dlzs,lz,lzd,dum1,m1,m2,rr,rb,mhefl,lhefl,thefl,lx
      real*8 tbgbf,thef,lbagbf,lheif,lhef,lzahbf
      real*8 rgbf,ragbf,rminf,mcgbf
      real*8 msp(200),gbp(200),c(5)
      external tbgbf,thef,lbagbf,lheif,lhef,lzahbf
      external rgbf,ragbf,rminf,mcgbf
*
      common /MSCFF/ msp
      common /GBCFF/ gbp
*
*
*       Initialization data set for fitting formulae.
*       -----------------------------------------------------------
*
      real*8 xz(76),xt(31),xl(72),xr(119),xg(112),xh(99)
*
* data for Lzams(1->35) and Rzams(36->76)
*
      data xz /
     &  3.970417d-01, -3.2913574d-01, 3.4776688d-01, 3.7470851d-01, 
     &  9.011915d-02,
     &  8.527626d+00,-2.441225973d+01, 5.643597107d+01, 3.706152575d+01, 
     &  5.4562406d+00,
     &  2.5546d-04, -1.23461d-03, -2.3246d-04,  4.5519d-04, 
     &  1.6176d-04,
     &  5.432889d+00, -8.62157806d+00, 1.344202049d+01, 
     &  1.451584135d+01, 3.39793084d+00,
     &  5.563579d+00,-1.032345224d+01, 1.944322980d+01, 
     &  1.897361347d+01, 4.16903097d+00,
     &  7.8866060d-01, -2.90870942d+00,  6.54713531d+00, 
     &  4.05606657d+00, 5.3287322d-01,
     &  5.86685d-03, -1.704237d-02, 3.872348d-02, 2.570041d-02, 
     &  3.83376d-03,
     &  1.715359d+00, 6.2246212d-01, -9.2557761d-01, -1.16996966d+00,
     & -3.0631491d-01,
     &  6.597788d+00, -4.2450044d-01,-1.213339427d+01,-1.073509484d+01,
     & -2.51487077d+00,
     &  1.008855000d+01, -7.11727086d+00,-3.167119479d+01,
     & -2.424848322d+01,-5.33608972d+00,
     &  1.012495d+00, 3.2699690d-01, -9.23418d-03, -3.876858d-02,
     & -4.12750d-03,
     &  7.490166d-02, 2.410413d-02, 7.233664d-02, 3.040467d-02, 
     &  1.97741d-03, 1.077422d-02,
     &  3.082234d+00, 9.447205d-01, -2.15200882d+00, -2.49219496d+00,
     & -6.3848738d-01,
     &  1.784778d+01, -7.4534569d+00,-4.896066856d+01,-4.005386135d+01,
     & -9.09331816d+00,
     &  2.2582d-04, -1.86899d-03, 3.88783d-03, 1.42402d-03,-7.671d-05/
*
* data for Tbgb(1->17) and Thook(18->31)
*
      data xt /1.593890d+03, 2.053038d+03, 1.231226d+03,
     &         2.327785d+02, 2.706708d+03, 1.483131d+03,
     &         5.772723d+02, 7.411230d+01, 1.466143d+02,
     &        -1.048442d+02,-6.795374d+01,-1.391127d+01,
     &         4.141960d-02, 4.564888d-02, 2.958542d-02,
     &         5.571483d-03, 3.426349d-01,
     &         1.949814d+01, 1.758178d+00,-6.008212d+00,
     &        -4.470533d+00, 4.903830d+00, 5.212154d-02,
     &         3.166411d-02,-2.750074d-03,-2.271549d-03,
     &         1.312179d+00,-3.294936d-01, 9.231860d-02,
     &         2.610989d-02, 8.073972d-01/
*
* data for Ltms(1->27), Lalpha(28->43), Lbeta(44->56) and Lhook(57->72)
*
      data xl /1.031538d+00,-2.434480d-01, 7.732821d+00,
     &         6.460705d+00, 1.374484d+00, 1.043715d+00,
     &        -1.577474d+00,-5.168234d+00,-5.596506d+00,
     &        -1.299394d+00, 7.859573d+02,-8.542048d+00,
     &        -2.642511d+01,-9.585707d+00, 3.858911d+03,
     &         2.459681d+03,-7.630093d+01,-3.486057d+02,
     &        -4.861703d+01, 2.888720d+02, 2.952979d+02,
     &         1.850341d+02, 3.797254d+01, 7.196580d+00,
     &         5.613746d-01, 3.805871d-01, 8.398728d-02,
     &         2.321400d-01, 1.828075d-03,-2.232007d-02,
     &        -3.378734d-03, 1.163659d-02, 3.427682d-03,
     &         1.421393d-03,-3.710666d-03, 1.048020d-02,
     &        -1.231921d-02,-1.686860d-02,-4.234354d-03,
     &         1.555590d+00,-3.223927d-01,-5.197429d-01,
     &        -1.066441d-01,
     &         3.855707d-01,-6.104166d-01, 5.676742d+00,
     &         1.060894d+01, 5.284014d+00, 3.579064d-01,
     &        -6.442936d-01, 5.494644d+00, 1.054952d+01,
     &         5.280991d+00, 9.587587d-01, 8.777464d-01,
     &         2.017321d-01,
     &         1.910302d-01, 1.158624d-01, 3.348990d-02,
     &         2.599706d-03, 3.931056d-01, 7.277637d-02,
     &        -1.366593d-01,-4.508946d-02, 3.267776d-01,
     &         1.204424d-01, 9.988332d-02, 2.455361d-02,
     &         5.990212d-01, 5.570264d-02, 6.207626d-02,
     &         1.777283d-02/
*
* data for Rtms(1->40), Ralpha(41->64), Rbeta(65->83), Rgamma(84->103)
* and Rhook(104->119)
*
      data xr /2.187715d-01,-2.154437d+00,-3.768678d+00,
     &        -1.975518d+00,-3.021475d-01, 1.466440d+00,
     &         1.839725d+00, 6.442199d+00, 4.023635d+00,
     &         6.957529d-01, 2.652091d+01, 8.178458d+01,
     &         1.156058d+02, 7.633811d+01, 1.950698d+01,
     &         1.472103d+00,-2.947609d+00,-3.312828d+00,
     &        -9.945065d-01, 3.071048d+00,-5.679941d+00,
     &        -9.745523d+00,-3.594543d+00,-8.672073d-02,
     &         2.617890d+00, 1.019135d+00,-3.292551d-02,
     &        -7.445123d-02, 1.075567d-02, 1.773287d-02,
     &         9.610479d-03, 1.732469d-03, 1.476246d+00,
     &         1.899331d+00, 1.195010d+00, 3.035051d-01,
     &         5.502535d+00,-6.601663d-02, 9.968707d-02,
     &         3.599801d-02,
     &         4.907546d-01,-1.683928d-01,-3.108742d-01,
     &        -7.202918d-02, 4.537070d+00,-4.465455d+00,
     &        -1.612690d+00,-1.623246d+00, 1.796220d+00,
     &         2.814020d-01, 1.423325d+00, 3.421036d-01,
     &         2.256216d+00, 3.773400d-01, 1.537867d+00,
     &         4.396373d-01, 1.564231d-03, 1.653042d-03,
     &        -4.439786d-03,-4.951011d-03,-1.216530d-03,
     &         5.210157d+00,-4.143695d+00,-2.120870d+00,
     &         1.071489d+00,-1.164852d-01,-8.623831d-02,
     &        -1.582349d-02, 7.108492d-01, 7.935927d-01,
     &         3.926983d-01, 3.622146d-02, 3.478514d+00,
     &        -2.585474d-02,-1.512955d-02,-2.833691d-03,
     &         3.969331d-03, 4.539076d-03, 1.720906d-03,
     &         1.897857d-04, 9.132108d-01,-1.653695d-01,
     &         3.636784d-02,
     &         1.192334d-02, 1.083057d-02, 1.230969d+00,
     &         1.551656d+00,-1.668868d-01, 5.818123d-01,
     &        -1.105027d+01,-1.668070d+01, 7.615495d-01,
     &         1.068243d-01,-2.011333d-01,-9.371415d-02,
     &        -1.015564d-01,-2.161264d-01,-5.182516d-02,
     &        -3.868776d-01,-5.457078d-01,-1.463472d-01,
     &         9.409838d+00, 1.522928d+00,
     &         7.330122d-01, 5.192827d-01, 2.316416d-01,
     &         8.346941d-03, 1.172768d+00,-1.209262d-01,
     &        -1.193023d-01,-2.859837d-02, 3.982622d-01,
     &        -2.296279d-01,-2.262539d-01,-5.219837d-02,
     &         3.571038d+00,-2.223625d-02,-2.611794d-02,
     &        -6.359648d-03/
*
* data for Lbgb(1->24), Lbagb(25->44), Rgb(45->66), Ragb(67->100),
* Mchei(101->102) and Mcbagb(103->112)
*
      data xg /9.511033d+01, 6.819618d+01,-1.045625d+01,
     &        -1.474939d+01, 3.113458d+01, 1.012033d+01,
     &        -4.650511d+00,-2.463185d+00, 1.413057d+00,
     &         4.578814d-01,-6.850581d-02,-5.588658d-02,
     &         3.910862d+01, 5.196646d+01, 2.264970d+01,
     &         2.873680d+00, 4.597479d+00,-2.855179d-01,
     &         2.709724d-01, 6.682518d+00, 2.827718d-01,
     &        -7.294429d-02, 4.637345d+00, 9.301992d+00,
     &         1.626062d+02,-1.168838d+01,-5.498343d+00,
     &         3.336833d-01,-1.458043d-01,-2.011751d-02,
     &         7.425137d+01, 1.790236d+01, 3.033910d+01,
     &         1.018259d+01, 9.268325d+02,-9.739859d+01,
     &        -7.702152d+01,-3.158268d+01, 1.127018d+01,
     &         1.622158d+00,-1.443664d+00,-9.474699d-01,
     &         2.474401d+00, 3.892972d-01,
     &         9.960283d-01, 8.164393d-01, 2.383830d+00,
     &         2.223436d+00, 8.638115d-01, 1.231572d-01,
     &         2.561062d-01, 7.072646d-02,-5.444596d-02,
     &        -5.798167d-02,-1.349129d-02, 1.157338d+00,
     &         1.467883d+00, 4.299661d+00, 3.130500d+00,
     &         6.992080d-01, 1.640687d-02, 4.022765d-01,
     &         3.050010d-01, 9.962137d-01, 7.914079d-01,
     &         1.728098d-01,
     &         1.125124d+00, 1.306486d+00, 3.622359d+00,
     &         2.601976d+00, 3.031270d-01,-1.343798d-01,
     &         3.349489d-01, 4.531269d-03, 1.131793d-01,
     &         2.300156d-01, 7.632745d-02, 1.467794d+00,
     &         2.798142d+00, 9.455580d+00, 8.963904d+00,
     &         3.339719d+00, 4.426929d-01, 4.658512d-01,
     &         2.597451d-01, 9.048179d-01, 7.394505d-01,
     &         1.607092d-01, 1.110866d+00, 9.623856d-01,
     &         2.735487d+00, 2.445602d+00, 8.826352d-01,
     &         1.140142d-01,-1.584333d-01,-1.728865d-01,
     &        -4.461431d-01,-3.925259d-01,-1.276203d-01,
     &        -1.308728d-02,
     &         9.796164d-02, 1.350554d+00,
     &         1.445216d-01,-6.180219d-02, 3.093878d-02,
     &         1.567090d-02, 1.304129d+00, 1.395919d-01,
     &         4.142455d-03,-9.732503d-03, 5.114149d-01,
     &        -1.160850d-02/
*
* data for Lhei(1->14), Lhe(15->25) Rmin(26->43), The(44->65), 
* Tbl(66->79), Lzahb(80->87) and Rzahb(88->99)
*
      data xh /2.751631d+03, 3.557098d+02,
     &        -3.820831d-02, 5.872664d-02, 1.5d+01,
     &         1.071738d+02,-8.970339d+01,-3.949739d+01,
     &         7.348793d+02,-1.531020d+02,-3.793700d+01,
     &         9.219293d+00,-2.005865d+00,-5.561309d-01,
     &         2.917412d+00, 1.575290d+00, 5.751814d-01,
     &         6.371760d-01, 3.880523d-01, 4.916389d+00,
     &         2.862149d+00, 7.844850d-01, 3.629118d+00,
     &        -9.112722d-01, 1.042291d+00,
     &         1.609901d+01, 7.391573d+00, 2.277010d+01,
     &         8.334227d+00, 1.747500d-01, 6.271202d-02,
     &        -2.324229d-02,-1.844559d-02, 2.752869d+00,
     &         2.729201d-02, 4.996927d-01, 2.496551d-01,
     &        -9.138012d-02,-3.671407d-01, 3.518506d+00,
     &         1.112440d+00,-4.556216d-01,-2.179426d-01,
     &         1.314955d+02, 2.009258d+01,-5.143082d-01,
     &        -1.379140d+00, 1.823973d+01,-3.074559d+00,
     &        -4.307878d+00, 1.5d1,
     &         2.327037d+00, 2.403445d+00, 1.208407d+00,
     &         2.087263d-01, 1.079113d-01, 1.762409d-02,
     &         1.096601d-02, 3.058818d-03, 2.327409d+00,
     &         6.901582d-01,-2.158431d-01,-1.084117d-01,
     &         1.997378d+00,-8.126205d-01,
     &         2.214315d+00,-1.975747d+00, 4.805428d-01,
     &         2.471620d+00,-5.401682d+00, 3.247361d+00,
     &         5.072525d+00, 1.146189d+01, 6.961724d+00,
     &         1.316965d+00, 5.139740d+00, 1.127733d+00,
     &         2.344416d-01,-3.793726d-01,
     &         5.496045d+01,-1.289968d+01, 6.385758d+00,
     &         1.832694d+00,-5.766608d-02, 5.696128d-02,
     &         1.211104d+02, 1.647903d+00,
     &         2.063983d+00, 7.363827d-01, 2.654323d-01,
     &        -6.140719d-02, 2.214088d+02, 2.187113d+02,
     &         1.170177d+01,-2.635340d+01, 2.003160d+00,
     &         9.388871d-01, 9.656450d-01, 2.362266d-01/
*
***

*
      data c /3.040581d-01, 8.049509d-02, 8.967485d-02,
     &        8.780198d-02, 2.219170d-02/
*
*       ------------------------------------------------------------
*
*      zpars:  1; M below which hook doesn't appear on MS, Mhook.
*              2; M above which He ignition occurs non-degenerately, Mhef.
*              3; M above which He ignition occurs on the HG, Mfgb. 
*              4; M below which C/O ignition doesn't occur, Mup.
*              5; M above which C ignites in the centre, Mec.
*              6; value of log D for M<= zpars(3)
*              7; value of x for Rgb propto M^(-x)
*              8; value of x for tMS = MAX(tHOOK,x*tBGB)
*              9; constant for McHeIf when computing Mc,BGB, mchefl.
*             10; constant for McHeIf when computing Mc,HeI, mchefl.
*             11; hydrogen abundance.
*             12; helium abundance.
*             13; constant x in rmin = rgb*x**y used by LM CHeB.
*             14; z**0.4 to be used for WD L formula.
*
*       ------------------------------------------------------------
*
      lzs = log10(z/0.02d0)
      dlzs = 1.d0/(z*log(10.d0))
      lz = log10(z)
      lzd = lzs + 1.d0
*
      zpars(1) = 1.0185d0 + lzs*(0.16015d0 + lzs*0.0892d0)
      zpars(2) = 1.995d0 + lzs*(0.25d0 + lzs*0.087d0)
      zpars(3) = 16.5d0*z**0.06d0/(1.d0 + (1.0d-04/z)**1.27d0)
      zpars(4) = MAX(6.11044d0 + 1.02167d0*lzs, 5.d0)
      zpars(5) = zpars(4) + 1.8d0
      zpars(6) = 5.37d0 + lzs*0.135d0
      zpars(7) = c(1) + lzs*(c(2) + lzs*(c(3) + lzs*(c(4) + lzs*c(5))))
      zpars(8) = MAX(0.95d0,MAX(0.95d0-(10.d0/3.d0)*(z-0.01d0),
     &           MIN(0.99d0,0.98d0-(100.d0/7.d0)*(z-0.001d0))))
***
* Lzams
      msp(1) = xz(1)+lzs*(xz(2)+lzs*(xz(3)+lzs*(xz(4)+lzs*xz(5))))
      msp(2) = xz(6)+lzs*(xz(7)+lzs*(xz(8)+lzs*(xz(9)+lzs*xz(10))))
      msp(3) = xz(11)+lzs*(xz(12)+lzs*(xz(13)+lzs*(xz(14)+lzs*xz(15))))
      msp(4) = xz(16)+lzs*(xz(17)+lzs*(xz(18)+lzs*(xz(19)+lzs*xz(20))))
      msp(5) = xz(21)+lzs*(xz(22)+lzs*(xz(23)+lzs*(xz(24)+lzs*xz(25))))
      msp(6) = xz(26)+lzs*(xz(27)+lzs*(xz(28)+lzs*(xz(29)+lzs*xz(30))))
      msp(7) = xz(31)+lzs*(xz(32)+lzs*(xz(33)+lzs*(xz(34)+lzs*xz(35))))
* Rzams
      msp(8) = xz(36)+lzs*(xz(37)+lzs*(xz(38)+lzs*(xz(39)+lzs*xz(40))))
      msp(9) = xz(41)+lzs*(xz(42)+lzs*(xz(43)+lzs*(xz(44)+lzs*xz(45))))
      msp(10) = xz(46)+lzs*(xz(47)+lzs*(xz(48)+lzs*(xz(49)+lzs*xz(50))))
      msp(11) = xz(51)+lzs*(xz(52)+lzs*(xz(53)+lzs*(xz(54)+lzs*xz(55))))
      msp(12) = xz(56)+lzs*(xz(57)+lzs*(xz(58)+lzs*(xz(59)+lzs*xz(60))))
      msp(13) = xz(61)
      msp(14) = xz(62)+lzs*(xz(63)+lzs*(xz(64)+lzs*(xz(65)+lzs*xz(66))))
      msp(15) = xz(67)+lzs*(xz(68)+lzs*(xz(69)+lzs*(xz(70)+lzs*xz(71))))
      msp(16) = xz(72)+lzs*(xz(73)+lzs*(xz(74)+lzs*(xz(75)+lzs*xz(76))))
* Tbgb 
      msp(17) = xt(1)+lzs*(xt(2)+lzs*(xt(3)+lzs*xt(4)))
      msp(18) = xt(5)+lzs*(xt(6)+lzs*(xt(7)+lzs*xt(8)))
      msp(19) = xt(9)+lzs*(xt(10)+lzs*(xt(11)+lzs*xt(12)))
      msp(20) = xt(13)+lzs*(xt(14)+lzs*(xt(15)+lzs*xt(16)))
      msp(21) = xt(17)
* dTbgb/dz
      msp(117) = dlzs*(xt(2)+lzs*(2.d0*xt(3)+3.d0*lzs*xt(4)))
      msp(118) = dlzs*(xt(6)+lzs*(2.d0*xt(7)+3.d0*lzs*xt(8)))
      msp(119) = dlzs*(xt(10)+lzs*(2.d0*xt(11)+3.d0*lzs*xt(12)))
      msp(120) = dlzs*(xt(14)+lzs*(2.d0*xt(15)+3.d0*lzs*xt(16)))
* Thook
      msp(22) = xt(18)+lzs*(xt(19)+lzs*(xt(20)+lzs*xt(21)))
      msp(23) = xt(22)
      msp(24) = xt(23)+lzs*(xt(24)+lzs*(xt(25)+lzs*xt(26)))
      msp(25) = xt(27)+lzs*(xt(28)+lzs*(xt(29)+lzs*xt(30)))
      msp(26) = xt(31)
* Ltms 
      msp(27) = xl(1)+lzs*(xl(2)+lzs*(xl(3)+lzs*(xl(4)+lzs*xl(5))))
      msp(28) = xl(6)+lzs*(xl(7)+lzs*(xl(8)+lzs*(xl(9)+lzs*xl(10))))
      msp(29) = xl(11)+lzs*(xl(12)+lzs*(xl(13)+lzs*xl(14)))
      msp(30) = xl(15)+lzs*(xl(16)+lzs*(xl(17)+lzs*(xl(18)+lzs*xl(19))))
      msp(27) = msp(27)*msp(30)
      msp(28) = msp(28)*msp(30)
      msp(31) = xl(20)+lzs*(xl(21)+lzs*(xl(22)+lzs*xl(23)))
      msp(32) = xl(24)+lzs*(xl(25)+lzs*(xl(26)+lzs*xl(27)))
* Lalpha
      m2 = 2.d0
      msp(33) = xl(28)+lzs*(xl(29)+lzs*(xl(30)+lzs*xl(31)))
      msp(34) = xl(32)+lzs*(xl(33)+lzs*(xl(34)+lzs*xl(35)))
      msp(35) = xl(36)+lzs*(xl(37)+lzs*(xl(38)+lzs*xl(39)))
      msp(36) = xl(40)+lzs*(xl(41)+lzs*(xl(42)+lzs*xl(43)))
      msp(37) = MAX(0.9d0,1.1064d0+lzs*(0.415d0+0.18d0*lzs))
      msp(38) = MAX(1.d0,1.19d0+lzs*(0.377d0+0.176d0*lzs))
      if(z.gt.0.01d0)then
         msp(37) = MIN(msp(37),1.d0)
         msp(38) = MIN(msp(38),1.1d0)
      endif
      msp(39) = MAX(0.145d0,0.0977d0-lzs*(0.231d0+0.0753d0*lzs))
      msp(40) = MIN(0.24d0+lzs*(0.18d0+0.595d0*lzs),0.306d0+0.053d0*lzs)
      msp(41) = MIN(0.33d0+lzs*(0.132d0+0.218d0*lzs),
     &              0.3625d0+0.062d0*lzs)
      msp(42) = (msp(33)+msp(34)*m2**msp(36))/
     &          (m2**0.4d0+msp(35)*m2**1.9d0)
* Lbeta
      msp(43) = xl(44)+lzs*(xl(45)+lzs*(xl(46)+lzs*(xl(47)+lzs*xl(48))))
      msp(44) = xl(49)+lzs*(xl(50)+lzs*(xl(51)+lzs*(xl(52)+lzs*xl(53))))
      msp(45) = xl(54)+lzs*(xl(55)+lzs*xl(56))
      msp(46) = MIN(1.4d0,1.5135d0+0.3769d0*lzs)
      msp(46) = MAX(0.6355d0-0.4192d0*lzs,MAX(1.25d0,msp(46)))
* Lhook
      msp(47) = xl(57)+lzs*(xl(58)+lzs*(xl(59)+lzs*xl(60)))
      msp(48) = xl(61)+lzs*(xl(62)+lzs*(xl(63)+lzs*xl(64)))
      msp(49) = xl(65)+lzs*(xl(66)+lzs*(xl(67)+lzs*xl(68)))
      msp(50) = xl(69)+lzs*(xl(70)+lzs*(xl(71)+lzs*xl(72)))
      msp(51) = MIN(1.4d0,1.5135d0+0.3769d0*lzs)
      msp(51) = MAX(0.6355d0-0.4192d0*lzs,MAX(1.25d0,msp(51)))
* Rtms
      msp(52) = xr(1)+lzs*(xr(2)+lzs*(xr(3)+lzs*(xr(4)+lzs*xr(5))))
      msp(53) = xr(6)+lzs*(xr(7)+lzs*(xr(8)+lzs*(xr(9)+lzs*xr(10))))
      msp(54) = xr(11)+lzs*(xr(12)+lzs*(xr(13)+lzs*(xr(14)+lzs*xr(15))))
      msp(55) = xr(16)+lzs*(xr(17)+lzs*(xr(18)+lzs*xr(19)))
      msp(56) = xr(20)+lzs*(xr(21)+lzs*(xr(22)+lzs*xr(23)))
      msp(52) = msp(52)*msp(54)
      msp(53) = msp(53)*msp(54)
      msp(57) = xr(24)
      msp(58) = xr(25)+lzs*(xr(26)+lzs*(xr(27)+lzs*xr(28)))
      msp(59) = xr(29)+lzs*(xr(30)+lzs*(xr(31)+lzs*xr(32)))
      msp(60) = xr(33)+lzs*(xr(34)+lzs*(xr(35)+lzs*xr(36)))
      msp(61) = xr(37)+lzs*(xr(38)+lzs*(xr(39)+lzs*xr(40)))
*
      msp(62) = MAX(0.097d0-0.1072d0*(lz+3.d0),MAX(0.097d0,MIN(0.1461d0,
     &              0.1461d0+0.1237d0*(lz+2.d0))))
      msp(62) = 10.d0**msp(62)
      m2 = msp(62) + 0.1d0
      msp(63) = (msp(52)+msp(53)*msp(62)**msp(55))/
     &          (msp(54)+msp(62)**msp(56))
      msp(64) = (msp(57)*m2**3+msp(58)*m2**msp(61)+
     &           msp(59)*m2**(msp(61)+1.5d0))/(msp(60)+m2**5)
* Ralpha
      msp(65) = xr(41)+lzs*(xr(42)+lzs*(xr(43)+lzs*xr(44)))
      msp(66) = xr(45)+lzs*(xr(46)+lzs*(xr(47)+lzs*xr(48)))
      msp(67) = xr(49)+lzs*(xr(50)+lzs*(xr(51)+lzs*xr(52)))
      msp(68) = xr(53)+lzs*(xr(54)+lzs*(xr(55)+lzs*xr(56)))
      msp(69) = xr(57)+lzs*(xr(58)+lzs*(xr(59)+lzs*(xr(60)+lzs*xr(61))))
      msp(70) = MAX(0.9d0,MIN(1.d0,1.116d0+0.166d0*lzs))
      msp(71) = MAX(1.477d0+0.296d0*lzs,MIN(1.6d0,-0.308d0-1.046d0*lzs))
      msp(71) = MAX(0.8d0,MIN(0.8d0-2.d0*lzs,msp(71)))
      msp(72) = xr(62)+lzs*(xr(63)+lzs*xr(64))
      msp(73) = MAX(0.065d0,0.0843d0-lzs*(0.0475d0+0.0352d0*lzs))
      msp(74) = 0.0736d0+lzs*(0.0749d0+0.04426d0*lzs)
      if(z.lt.0.004d0) msp(74) = MIN(0.055d0,msp(74))
      msp(75) = MAX(0.091d0,MIN(0.121d0,0.136d0+0.0352d0*lzs))
      msp(76) = (msp(65)*msp(71)**msp(67))/(msp(66) + msp(71)**msp(68))
      if(msp(70).gt.msp(71))then
         msp(70) = msp(71)
         msp(75) = msp(76)
      endif
* Rbeta
      msp(77) = xr(65)+lzs*(xr(66)+lzs*(xr(67)+lzs*xr(68)))
      msp(78) = xr(69)+lzs*(xr(70)+lzs*(xr(71)+lzs*xr(72)))
      msp(79) = xr(73)+lzs*(xr(74)+lzs*(xr(75)+lzs*xr(76)))
      msp(80) = xr(77)+lzs*(xr(78)+lzs*(xr(79)+lzs*xr(80)))
      msp(81) = xr(81)+lzs*(xr(82)+lzs*lzs*xr(83))
      if(z.gt.0.01d0) msp(81) = MAX(msp(81),0.95d0)
      msp(82) = MAX(1.4d0,MIN(1.6d0,1.6d0+lzs*(0.764d0+0.3322d0*lzs)))
* Rgamma
      msp(83) = MAX(xr(84)+lzs*(xr(85)+lzs*(xr(86)+lzs*xr(87))),
     &              xr(96)+lzs*(xr(97)+lzs*xr(98)))
      msp(84) = MIN(0.d0,xr(88)+lzs*(xr(89)+lzs*(xr(90)+lzs*xr(91))))
      msp(84) = MAX(msp(84),xr(99)+lzs*(xr(100)+lzs*xr(101)))
      msp(85) = xr(92)+lzs*(xr(93)+lzs*(xr(94)+lzs*xr(95)))
      msp(85) = MAX(0.d0,MIN(msp(85),7.454d0+9.046d0*lzs))
      msp(86) = MIN(xr(102)+lzs*xr(103),MAX(2.d0,-13.3d0-18.6d0*lzs))
      msp(87) = MIN(1.5d0,MAX(0.4d0,2.493d0+1.1475d0*lzs))
      msp(88) = MAX(1.d0,MIN(1.27d0,0.8109d0-0.6282d0*lzs))
      msp(88) = MAX(msp(88),0.6355d0-0.4192d0*lzs)
      msp(89) = MAX(5.855420d-02,-0.2711d0-lzs*(0.5756d0+0.0838d0*lzs))
* Rhook
      msp(90) = xr(104)+lzs*(xr(105)+lzs*(xr(106)+lzs*xr(107)))
      msp(91) = xr(108)+lzs*(xr(109)+lzs*(xr(110)+lzs*xr(111)))
      msp(92) = xr(112)+lzs*(xr(113)+lzs*(xr(114)+lzs*xr(115)))
      msp(93) = xr(116)+lzs*(xr(117)+lzs*(xr(118)+lzs*xr(119)))
      msp(94) = MIN(1.25d0,
     &          MAX(1.1d0,1.9848d0+lzs*(1.1386d0+0.3564d0*lzs)))
      msp(95) = 0.063d0 + lzs*(0.0481d0 + 0.00984d0*lzs)
      msp(96) = MIN(1.3d0,MAX(0.45d0,1.2d0+2.45d0*lzs))
* Lneta
      if(z.gt.0.0009d0)then
         msp(97) = 10.d0
      else
         msp(97) = 20.d0
      endif
* Lbgb 
      gbp(1) = xg(1)+lzs*(xg(2)+lzs*(xg(3)+lzs*xg(4)))
      gbp(2) = xg(5)+lzs*(xg(6)+lzs*(xg(7)+lzs*xg(8)))
      gbp(3) = xg(9)+lzs*(xg(10)+lzs*(xg(11)+lzs*xg(12)))
      gbp(4) = xg(13)+lzs*(xg(14)+lzs*(xg(15)+lzs*xg(16)))
      gbp(5) = xg(17)+lzs*(xg(18)+lzs*xg(19))
      gbp(6) = xg(20)+lzs*(xg(21)+lzs*xg(22))
      gbp(3) = gbp(3)**gbp(6)
      gbp(7) = xg(23)
      gbp(8) = xg(24)
* Lbagb
* set gbp(16) = 1.d0 until it is reset later with an initial
* call to Lbagbf using mass = zpars(2) and mhefl = 0.0
      gbp(9) = xg(25) + lzs*(xg(26) + lzs*xg(27))
      gbp(10) = xg(28) + lzs*(xg(29) + lzs*xg(30))
      gbp(11) = 15.d0
      gbp(12) = xg(31)+lzs*(xg(32)+lzs*(xg(33)+lzs*xg(34)))
      gbp(13) = xg(35)+lzs*(xg(36)+lzs*(xg(37)+lzs*xg(38)))
      gbp(14) = xg(39)+lzs*(xg(40)+lzs*(xg(41)+lzs*xg(42)))
      gbp(15) = xg(43)+lzs*xg(44)
      gbp(12) = gbp(12)**gbp(15)
      gbp(14) = gbp(14)**gbp(15)
      gbp(16) = 1.d0
* Rgb
      gbp(17) = -4.6739d0-0.9394d0*lz
      gbp(17) = 10.d0**gbp(17)
      gbp(17) = MAX(gbp(17),-0.04167d0+55.67d0*z)
      gbp(17) = MIN(gbp(17),0.4771d0-9329.21d0*z**2.94d0)
      gbp(18) = MIN(0.54d0,0.397d0+lzs*(0.28826d0+0.5293d0*lzs))
      gbp(19) = MAX(-0.1451d0,-2.2794d0-lz*(1.5175d0+0.254d0*lz))
      gbp(19) = 10.d0**gbp(19)
      if(z.gt.0.004d0)then
         gbp(19) = MAX(gbp(19),0.7307d0+14265.1d0*z**3.395d0)
      endif
      gbp(20) = xg(45)+lzs*(xg(46)+lzs*(xg(47)+lzs*(xg(48)+
     &          lzs*(xg(49)+lzs*xg(50)))))
      gbp(21) = xg(51)+lzs*(xg(52)+lzs*(xg(53)+lzs*(xg(54)+lzs*xg(55))))
      gbp(22) = xg(56)+lzs*(xg(57)+lzs*(xg(58)+lzs*(xg(59)+
     &          lzs*(xg(60)+lzs*xg(61)))))
      gbp(23) = xg(62)+lzs*(xg(63)+lzs*(xg(64)+lzs*(xg(65)+lzs*xg(66))))
* Ragb
      gbp(24) = MIN(0.99164d0-743.123d0*z**2.83d0,
     &              1.0422d0+lzs*(0.13156d0+0.045d0*lzs))
      gbp(25) = xg(67)+lzs*(xg(68)+lzs*(xg(69)+lzs*(xg(70)+
     &          lzs*(xg(71)+lzs*xg(72)))))
      gbp(26) = xg(73)+lzs*(xg(74)+lzs*(xg(75)+lzs*(xg(76)+lzs*xg(77))))
      gbp(27) = xg(78)+lzs*(xg(79)+lzs*(xg(80)+lzs*(xg(81)+
     &          lzs*(xg(82)+lzs*xg(83)))))
      gbp(28) = xg(84)+lzs*(xg(85)+lzs*(xg(86)+lzs*(xg(87)+lzs*xg(88))))
      gbp(29) = xg(89)+lzs*(xg(90)+lzs*(xg(91)+lzs*(xg(92)+
     &          lzs*(xg(93)+lzs*xg(94)))))
      gbp(30) = xg(95)+lzs*(xg(96)+lzs*(xg(97)+lzs*(xg(98)+
     &          lzs*(xg(99)+lzs*xg(100)))))
      m1 = zpars(2) - 0.2d0
      gbp(31) = gbp(29) + gbp(30)*m1
      gbp(32) = MIN(gbp(25)/zpars(2)**gbp(26),gbp(27)/zpars(2)**gbp(28))
* Mchei
      gbp(33) = xg(101)**4
      gbp(34) = xg(102)*4.d0
* Mcagb
      gbp(35) = xg(103)+lzs*(xg(104)+lzs*(xg(105)+lzs*xg(106)))
      gbp(36) = xg(107)+lzs*(xg(108)+lzs*(xg(109)+lzs*xg(110)))
      gbp(37) = xg(111)+lzs*xg(112)
      gbp(35) = gbp(35)**4
      gbp(36) = gbp(36)*4.d0
      gbp(37) = gbp(37)**4
* Lhei
* set gbp(41) = -1.d0 until it is reset later with an initial
* call to Lheif using mass = zpars(2) and mhefl = 0.0
      gbp(38) = xh(1)+lzs*xh(2)
      gbp(39) = xh(3)+lzs*xh(4)
      gbp(40) = xh(5)
      gbp(41) = -1.d0
      gbp(42) = xh(6)+lzs*(xh(7)+lzs*xh(8))
      gbp(43) = xh(9)+lzs*(xh(10)+lzs*xh(11))
      gbp(44) = xh(12)+lzs*(xh(13)+lzs*xh(14))
      gbp(42) = gbp(42)**2
      gbp(44) = gbp(44)**2
* Lhe
      gbp(45) = xh(15)+lzs*(xh(16)+lzs*xh(17))
      if(lzs.gt.-1.d0)then
         gbp(46) = 1.d0 - xh(19)*(lzs+1.d0)**xh(18)
      else
         gbp(46) = 1.d0
      endif
      gbp(47) = xh(20)+lzs*(xh(21)+lzs*xh(22))
      gbp(48) = xh(23)+lzs*(xh(24)+lzs*xh(25))
      gbp(45) = gbp(45)**gbp(48)
      gbp(47) = gbp(47)**gbp(48)
      gbp(46) = gbp(46)/zpars(3)**0.1d0+(gbp(46)*gbp(47)-gbp(45))/
     &          zpars(3)**(gbp(48)+0.1d0)
* Rmin
      gbp(49) = xh(26)+lzs*(xh(27)+lzs*(xh(28)+lzs*xh(29)))
      gbp(50) = xh(30)+lzs*(xh(31)+lzs*(xh(32)+lzs*xh(33)))
      gbp(51) = xh(34)+lzs*(xh(35)+lzs*(xh(36)+lzs*xh(37)))
      gbp(52) = 5.d0+xh(38)*z**xh(39)
      gbp(53) = xh(40)+lzs*(xh(41)+lzs*(xh(42)+lzs*xh(43)))
      gbp(49) = gbp(49)**gbp(53)
      gbp(51) = gbp(51)**(2.d0*gbp(53))
* The
* set gbp(57) = -1.d0 until it is reset later with an initial
* call to Thef using mass = zpars(2), mc = 0.0  and mhefl = 0.0
      gbp(54) = xh(44)+lzs*(xh(45)+lzs*(xh(46)+lzs*xh(47)))
      gbp(55) = xh(48)+lzs*(xh(49)+lzs*xh(50))
      gbp(55) = MAX(gbp(55),1.d0)
      gbp(56) = xh(51)
      gbp(57) = -1.d0
      gbp(58) = xh(52)+lzs*(xh(53)+lzs*(xh(54)+lzs*xh(55)))
      gbp(59) = xh(56)+lzs*(xh(57)+lzs*(xh(58)+lzs*xh(59)))
      gbp(60) = xh(60)+lzs*(xh(61)+lzs*(xh(62)+lzs*xh(63)))
      gbp(61) = xh(64)+lzs*xh(65)
      gbp(58) = gbp(58)**gbp(61)
      gbp(60) = gbp(60)**5
* Tbl
      dum1 = zpars(2)/zpars(3)
      gbp(62) = xh(66)+lzs*xh(67)
      gbp(62) = -gbp(62)*log10(dum1)
      gbp(63) = xh(68)
      if(lzd.gt.0.d0) then
         gbp(64) = 1.d0-lzd*(xh(69)+lzd*(xh(70)+lzd*xh(71)))
      else
         gbp(64) = 1.d0
      end if
      gbp(65) = 1.d0-gbp(64)*dum1**gbp(63)
      gbp(66) = 1.d0 - lzd*(xh(77) + lzd*(xh(78) + lzd*xh(79)))
      gbp(67) = xh(72) + lzs*(xh(73) + lzs*(xh(74) + lzs*xh(75)))
      gbp(68) = xh(76)
* Lzahb
      gbp(69) = xh(80) + lzs*(xh(81) + lzs*xh(82))
      gbp(70) = xh(83) + lzs*(xh(84) + lzs*xh(85))
      gbp(71) = 15.d0
      gbp(72) = xh(86)
      gbp(73) = xh(87)
* Rzahb
      gbp(75) = xh(88) + lzs*(xh(89) + lzs*(xh(90) + lzs*xh(91)))
      gbp(76) = xh(92) + lzs*(xh(93) + lzs*(xh(94) + lzs*xh(95)))
      gbp(77) = xh(96) + lzs*(xh(97) + lzs*(xh(98) + lzs*xh(99)))
***
* finish Lbagb
      mhefl = 0.d0
      lx = lbagbf(zpars(2),mhefl)
      gbp(16) = lx
* finish LHeI
      dum1 = 0.d0
      lhefl = lheif(zpars(2),mhefl)
      gbp(41) = (gbp(38)*zpars(2)**gbp(39)-lhefl)/
     &          (EXP(zpars(2)*gbp(40))*lhefl)
* finish THe
      thefl = thef(zpars(2),dum1,mhefl)*tbgbf(zpars(2))
      gbp(57) = (thefl-gbp(54))/(gbp(54)*EXP(gbp(56)*zpars(2)))
* finish Tblf
      rb = ragbf(zpars(3),lheif(zpars(3),zpars(2)),mhefl)
      rr = 1.d0 - rminf(zpars(3))/rb
      rr = MAX(rr,1.0d-12)
      gbp(66) = gbp(66)/(zpars(3)**gbp(67)*rr**gbp(68))
* finish Lzahb
      gbp(74) = lhefl*lHef(zpars(2))
***
      kw = 0
      tm = 0.d0
      tn = 0.d0
      CALL star(kw,zpars(2),zpars(2),tm,tn,tscls,lums,GB,zpars)
      zpars(9) = mcgbf(lums(3),GB,lums(6))
      zpars(10) = mcgbf(lums(4),GB,lums(6))
* set the hydrogen and helium abundances
      zpars(11) = 0.76d0 - 3.d0*z
      zpars(12) = 0.24d0 + 2.d0*z
* set constant for low-mass CHeB stars
      zpars(13) = rminf(zpars(2))/
     &            rgbf(zpars(2),lzahbf(zpars(2),zpars(9),zpars(2)))
* 
      zpars(14) = z**0.4d0
*
      return
      end
***

***
      real*8 FUNCTION lzamsf(m)
      implicit none
      real*8 m,mx,a(200)
      common /MSCFF/ a
*
* A function to evaluate Lzams
* ( from Tout et al., 1996, MNRAS, 281, 257 ).
*
      mx = SQRT(m)
      lzamsf = (a(1)*m**5*mx + a(2)*m**11)/
     &         (a(3) + m**3 + a(4)*m**5 + a(5)*m**7 +
     &          a(6)*m**8 + a(7)*m**9*mx)
*
      return
      end
***
      real*8 FUNCTION rzamsf(m)
      implicit none
      real*8 m,mx,a(200)
      common /MSCFF/ a
*
* A function to evaluate Rzams
* ( from Tout et al., 1996, MNRAS, 281, 257 ).
*
      mx = SQRT(m)
      rzamsf = ((a(8)*m**2 + a(9)*m**6)*mx + a(10)*m**11 +
     &          (a(11) + a(12)*mx)*m**19)/
     &         (a(13) + a(14)*m**2 + 
     &          (a(15)*m**8 + m**18 + a(16)*m**19)*mx)
*
      return
      end
***
      real*8 FUNCTION tbgbf(m)
      implicit none
      real*8 m,a(200)
      common /MSCFF/ a
*
* A function to evaluate the lifetime to the BGB or to
* Helium ignition if no FGB exists.
* (JH 24/11/97)
*
      tbgbf = (a(17) + a(18)*m**4 + a(19)*m**(11.d0/2.d0) + m**7)/
     &        (a(20)*m**2 + a(21)*m**7)
*
      return
      end
***
      real*8 FUNCTION tbgbdf(m)
      implicit none
      real*8 m,mx,f,df,g,dg,a(200)
      common /MSCFF/ a
*
* A function to evaluate the derivitive of the lifetime to the BGB
* (or to Helium ignition if no FGB exists) wrt mass.
* (JH 24/11/97)
*
      mx = SQRT(m)
      f = a(17) + a(18)*m**4 + a(19)*m**5*mx + m**7
      df = 4.d0*a(18)*m**3 + 5.5d0*a(19)*m**4*mx + 7.d0*m**6
      g = a(20)*m**2 + a(21)*m**7
      dg = 2.d0*a(20)*m + 7.d0*a(21)*m**6
      tbgbdf = (df*g - f*dg)/(g*g)
*
      return
      end
***
      real*8 FUNCTION tbgdzf(m)
      implicit none
      real*8 m,mx,f,df,g,dg,a(200)
      common /MSCFF/ a
*
* A function to evaluate the derivitive of the lifetime to the BGB
* (or to Helium ignition if no FGB exists) wrt Z.
* (JH 14/12/98)
*
      mx = m**5*SQRT(m)
      f = a(17) + a(18)*m**4 + a(19)*mx + m**7
      df = a(117) + a(118)*m**4 + a(119)*mx
      g = a(20)*m**2 + a(21)*m**7
      dg = a(120)*m**2
      tbgdzf = (df*g - f*dg)/(g*g)
*
      return
      end
***
      real*8 FUNCTION thookf(m)
      implicit none
      real*8 m,a(200)
      common /MSCFF/ a
*
* A function to evaluate the lifetime to the end of the MS
* hook ( for those models that have one ) as a fraction of 
* the lifetime to the BGB
* Note that this function is only valid for M > Mhook.
* (JH 24/11/97)
*
      thookf = 1.d0 - 0.01d0*MAX(a(22)/m**a(23),a(24)+a(25)/m**a(26))
      thookf = MAX(thookf,0.5d0)
*
      return
      end
***
      real*8 FUNCTION ltmsf(m)
      implicit none
      real*8 m,a(200)
      common /MSCFF/ a
*
* A function to evaluate the luminosity at the end of the MS
* (JH 24/11/97)
*
      ltmsf = (a(27)*m**3 + a(28)*m**4 + a(29)*m**(a(32)+1.8d0))/
     &        (a(30) + a(31)*m**5 + m**a(32))
* 
      return
      end
***
      real*8 FUNCTION lalphf(m)
      implicit none
      real*8 m,mcut,a(200)
      common /MSCFF/ a
*
* A function to evaluate the Luminosity alpha coefficent.
* (JH 24/11/97)
*
      mcut = 2.d0
      if(m.ge.mcut)then
         lalphf = (a(33) + a(34)*m**a(36))/(m**0.4d0 + a(35)*m**1.9d0)
      else
         if(m.le.0.5d0)then
            lalphf = a(39)
         elseif(m.le.0.7d0)then
            lalphf = a(39) + ((0.3d0 - a(39))/0.2d0)*(m - 0.5d0)
         elseif(m.le.a(37))then
            lalphf = 0.3d0 + ((a(40)-0.3d0)/(a(37)-0.7d0))*(m - 0.7d0)
         elseif(m.le.a(38))then
            lalphf = a(40) + ((a(41)-a(40))/(a(38)-a(37)))*(m - a(37))
         else
            lalphf = a(41) + ((a(42)-a(41))/(mcut-a(38)))*(m - a(38))
         endif
      endif
*
      return
      end
***
      real*8 FUNCTION lbetaf(m)
      implicit none
      real*8 m,a1,a(200)
      common /MSCFF/ a
*
* A function to evaluate the Luminosity beta coefficent.
* (JH 24/11/97)
*
      lbetaf = a(43) - a(44)*m**a(45)
      lbetaf = MAX(lbetaf,0.d0)
      if(m.gt.a(46).and.lbetaf.gt.0.d0)then
         a1 = a(43) - a(44)*a(46)**a(45)
         lbetaf = a1 - 10.d0*a1*(m - a(46))
         lbetaf = MAX(lbetaf,0.d0)
      endif
*
      return
      end
***
      real*8 FUNCTION lnetaf(m)
      implicit none
      real*8 m,a(200)
      common /MSCFF/ a
*
* A function to evaluate the Luminosity neta exponent.
* (JH 24/11/97)
*
      if(m.le.1.d0)then
         lnetaf = 10.d0
      elseif(m.ge.1.1d0)then
         lnetaf = 20.d0
      else
         lnetaf = 10.d0 + 100.d0*(m - 1.d0)
      endif
      lnetaf = MIN(lnetaf,a(97))
*
      return
      end
***
      real*8 FUNCTION lhookf(m,mhook)
      implicit none
      real*8 m,mhook,a2,a(200)
      common /MSCFF/ a
*
* A function to evalute the luminosity at the start of
* the MS hook ( for those stars that have one ).
* Note that this function is only valid for M > Mhook.
* (JH 24/11/97)
*
      if(m.le.mhook)then
         lhookf = 0.d0
      elseif(m.ge.a(51))then
         lhookf = MIN(a(47)/m**a(48),a(49)/m**a(50))
      else
         a2 = MIN(a(47)/a(51)**a(48),a(49)/a(51)**a(50))
         lhookf = a2*((m-mhook)/(a(51)-mhook))**0.4d0
      endif
*
      return
      end
***
      real*8 FUNCTION rtmsf(m)
      implicit none
      real*8 m,m2,rchk,a(200)
      common /MSCFF/ a
      real*8 rzamsf
      external rzamsf
*
* A function to evaluate the radius at the end of the MS
* Note that a safety check is added to ensure Rtms > Rzams
* when extrapolating the function to low masses. 
* (JH 24/11/97)
*
      m2 = a(62) + 0.1d0
      if(m.le.a(62))then
         rchk = 1.5d0*rzamsf(m)
         rtmsf = MAX(rchk,(a(52) + a(53)*m**a(55))/(a(54) + m**a(56)))
      elseif(m.ge.m2)then
         rtmsf = (a(57)*m**3+a(58)*m**a(61)+a(59)*m**(a(61)+1.5d0))/
     &           (a(60) + m**5)
      else
         rtmsf = a(63) + ((a(64) - a(63))/0.1d0)*(m - a(62))
      endif
* 
      return
      end
***
      real*8 FUNCTION ralphf(m)
      implicit none
      real*8 m,a5,a(200)
      common /MSCFF/ a
*
* A function to evaluate the radius alpha coefficent.
* (JH 24/11/97)
*
      if(m.le.0.5d0)then
         ralphf = a(73)
      elseif(m.le.0.65d0)then
         ralphf = a(73) + ((a(74) - a(73))/0.15d0)*(m - 0.5d0)
      elseif(m.le.a(70))then
         ralphf = a(74) + ((a(75)-a(74))/(a(70)-0.65d0))*(m - 0.65d0)
      elseif(m.le.a(71))then
         ralphf = a(75) + ((a(76) - a(75))/(a(71) - a(70)))*(m - a(70))
      elseif(m.le.a(72))then
         ralphf = (a(65)*m**a(67))/(a(66) + m**a(68))
      else
         a5 = (a(65)*a(72)**a(67))/(a(66) + a(72)**a(68))
         ralphf = a5 + a(69)*(m - a(72))
      endif
*
      return
      end
***
      real*8 FUNCTION rbetaf(m)
      implicit none
      real*8 m,m2,m3,b2,b3,a(200)
      common /MSCFF/ a
*
* A function to evaluate the radius beta coefficent.
* (JH 24/11/97)
*
      m2 = 2.d0
      m3 = 16.d0
      if(m.le.1.d0)then
         rbetaf = 1.06d0
      elseif(m.le.a(82))then
         rbetaf = 1.06d0 + ((a(81)-1.06d0)/(a(82)-1.d0))*(m-1.d0)
      elseif(m.le.m2)then
         b2 = (a(77)*m2**(7.d0/2.d0))/(a(78) + m2**a(79))
         rbetaf = a(81) + ((b2-a(81))/(m2-a(82)))*(m-a(82))
      elseif(m.le.m3)then
         rbetaf = (a(77)*m**(7.d0/2.d0))/(a(78) + m**a(79))
      else
         b3 = (a(77)*m3**(7.d0/2.d0))/(a(78) + m3**a(79))
         rbetaf = b3 + a(80)*(m - m3)
      endif
      rbetaf = rbetaf - 1.d0
*
      return
      end
***
      real*8 FUNCTION rgammf(m)
      implicit none
      real*8 m,m1,b1,a(200)
      common /MSCFF/ a
*
* A function to evaluate the radius gamma coefficent.
* (JH 24/11/97)
*
      m1 = 1.d0
      if(m.gt.(a(88)+0.1d0))then
         rgammf = 0.d0
      else
         b1 = MAX(0.d0,a(83) + a(84)*(m1-a(85))**a(86))
         if(m.le.m1)then
            rgammf = a(83) + a(84)*ABS(m-a(85))**a(86)
         elseif(m.le.a(88))then
            rgammf = b1 + (a(89) - b1)*((m - m1)/(a(88) - m1))**a(87)
         else
            if(a(88).gt.m1) b1 = a(89)
            rgammf = b1 - 10.d0*b1*(m - a(88))
         endif
         rgammf = MAX(rgammf,0.d0)
      endif
*
      return
      end
***
      real*8 FUNCTION rhookf(m,mhook)
      implicit none
      real*8 m,mhook,m2,b2,a(200)
      common /MSCFF/ a
*
* A function to evalute the radius at the start of
* the MS hook ( for those stars that have one ).
* Note that this function is only valid for M > Mhook.
* (JH 24/11/97)
*
      if(m.le.mhook)then
         rhookf = 0.d0
      elseif(m.le.a(94))then
         rhookf = a(95)*SQRT((m-mhook)/(a(94)-mhook))
      elseif(m.le.2.d0)then
         m2 = 2.d0
         b2 = (a(90) + a(91)*m2**(7.d0/2.d0))/
     &        (a(92)*m2**3 + m2**a(93)) - 1.d0
         rhookf = a(95) + (b2-a(95))*((m-a(94))/(m2-a(94)))**a(96)
      else
         rhookf = (a(90) + a(91)*m**(7.d0/2.d0))/
     &            (a(92)*m**3 + m**a(93)) - 1.d0
      endif
*
      return
      end
***
      real*8 FUNCTION lbgbf(m)
      real*8 m,a(200)
      common /GBCFF/ a
*
* A function to evaluate the luminosity at the end of the 
* FGB ( for those models that have one )
* Note that this function is only valid for LM & IM stars
* (JH 24/11/97)
*
      lbgbf = (a(1)*m**a(5) + a(2)*m**a(8))/
     &        (a(3) + a(4)*m**a(7) + m**a(6))
* 
      return
      end
***
      real*8 FUNCTION lbgbdf(m)
      real*8 m,a(200)
      real*8 f,df,g,dg
      common /GBCFF/ a
*
* A function to evaluate the derivitive of the Lbgb function.
* Note that this function is only valid for LM & IM stars
* (JH 24/11/97)
*
      f = a(1)*m**a(5) + a(2)*m**a(8)
      df = a(5)*a(1)*m**(a(5)-1.d0) + a(8)*a(2)*m**(a(8)-1.d0)
      g = a(3) + a(4)*m**a(7) + m**a(6)
      dg = a(7)*a(4)*m**(a(7)-1.d0) + a(6)*m**(a(6)-1.d0)
*
      lbgbdf = (df*g - f*dg)/(g*g)
* 
      return
      end
***
      real*8 FUNCTION lbagbf(m,mhefl)
      implicit none
      real*8 m,mhefl,a4,a(200)
      common /GBCFF/ a
*
* A function to evaluate the BAGB luminosity. (OP 21/04/98)
* Continuity between LM and IM functions is ensured by setting
* gbp(16) = lbagbf(mhefl,0.0) with gbp(16) = 1.0.
*
      a4 = (a(9)*mhefl**a(10) - a(16))/(exp(mhefl*a(11))*a(16))
*
      if(m.lt.mhefl)then
         lbagbf = a(9)*m**a(10)/(1.d0 + a4*exp(m*a(11)))
      else
         lbagbf = (a(12) + a(13)*m**(a(15)+1.8d0))/(a(14) + m**a(15))
      endif
*
      return
      end
***
      real*8 FUNCTION rgbf(m,lum)
      implicit none
      real*8 m,lum,a1,a(200)
      common /GBCFF/ a
*
* A function to evaluate radius on the GB.
* (JH 24/11/97)
*
      a1 = MIN(a(20)/m**a(21),a(22)/m**a(23))
      rgbf = a1*(lum**a(18) + a(17)*lum**a(19))
*
      return
      end
***
      real*8 FUNCTION rgbdf(m,lum)
      implicit none
      real*8 m,lum,a1,a(200)
      common /GBCFF/ a
*
* A function to evaluate radius derivitive on the GB (as f(L)).
* (JH 24/11/97)
*
      a1 = MIN(a(20)/m**a(21),a(22)/m**a(23))
      rgbdf = a1*(a(18)*lum**(a(18)-1.d0) + 
     &            a(17)*a(19)*lum**(a(19)-1.d0))
*
      return
      end
***
      real*8 FUNCTION ragbf(m,lum,mhelf)
      implicit none
      real*8 m,lum,mhelf,m1,a1,a4,xx,a(200)
      common /GBCFF/ a
*
* A function to evaluate radius on the AGB.
* (JH 24/11/97)
*
      m1 = mhelf - 0.2d0
      if(m.ge.mhelf)then
         xx = a(24)
      elseif(m.ge.m1)then
         xx = 1.d0 + 5.d0*(a(24)-1.d0)*(m-m1)
      else
         xx = 1.d0
      endif
      a4 = xx*a(19)
      if(m.le.m1)then
         a1 = a(29) + a(30)*m
      elseif(m.ge.mhelf)then
         a1 = MIN(a(25)/m**a(26),a(27)/m**a(28))
      else
         a1 = a(31) + 5.d0*(a(32)-a(31))*(m-m1)
      endif
*
      ragbf = a1*(lum**a(18) + a(17)*lum**a4)
*
      return
      end
***
      real*8 FUNCTION ragbdf(m,lum,mhelf)
      implicit none
      real*8 m,lum,mhelf,m1,a1,a4,xx,a(200)
      common /GBCFF/ a
*
* A function to evaluate radius derivitive on the AGB (as f(L)).
* (JH 24/11/97)
*
      m1 = mhelf - 0.2d0
      if(m.ge.mhelf)then
         xx = a(24)
      elseif(m.ge.m1)then
         xx = 1.d0 + 5.d0*(a(24)-1.d0)*(m-m1)
      else
         xx = 1.d0
      endif
      a4 = xx*a(19)
      if(m.le.m1)then
         a1 = a(29) + a(30)*m
      elseif(m.ge.mhelf)then
         a1 = MIN(a(25)/m**a(26),a(27)/m**a(28))
      else
         a1 = a(31) + 5.d0*(a(32)-a(31))*(m-m1)
      endif
*
      ragbdf = a1*(a(18)*lum**(a(18)-1.d0) + 
     &             a(17)*a4*lum**(a4-1.d0))
*
      return
      end
***
      real*8 FUNCTION mctmsf(m)
      implicit none
      real*8 m,m525
*
* A function to evaluate core mass at the end of the MS as a 
* fraction of the BGB value, i.e. this must be multiplied by 
* the BGB value (see below) to give the actual core mass (JH 5/9/99)
*
      m525 = m**(21.d0/4.d0)
      mctmsf = (1.586d0 + m525)/(2.434d0 + 1.02d0*m525)
*
      return
      end
***
      real*8 FUNCTION mcheif(m,mhefl,mchefl)
      implicit none
      real*8 m,mhefl,mchefl,mcbagb,a3,a(200)
      common /GBCFF/ a
      real*8 mcagbf
      external mcagbf
*
* A function to evaluate core mass at BGB or He ignition
* (depending on mchefl) for IM & HM stars  (OP 25/11/97)
*
      mcbagb = mcagbf(m)
      a3 = mchefl**4 - a(33)*mhefl**a(34)
      mcheif = MIN(0.95d0*mcbagb,(a3 + a(33)*m**a(34))**(1.d0/4.d0))
*
      return
      end
***
      real*8 FUNCTION mheif(mc,mhefl,mchefl)
      implicit none
      real*8 mc,mhefl,mchefl,m1,m2,a3,a(200)
      common /GBCFF/ a
      real*8 mbagbf
      external mbagbf
*
* A function to evaluate mass at BGB or He ignition
* (depending on mchefl) for IM & HM stars by inverting
* mcheif
*
      m1 = mbagbf(mc/0.95d0)
      a3 = mchefl**4 - a(33)*mhefl**a(34)
      m2 = ((mc**4 - a3)/a(33))**(1.d0/a(34))
      mheif = MAX(m1,m2)
*
      return
      end
***
      real*8 FUNCTION mcagbf(m)
      implicit none
      real*8 m,a(200)
      common /GBCFF/ a
*
* A function to evaluate core mass at the BAGB  (OP 25/11/97)
*
      mcagbf = (a(37) + a(35)*m**a(36))**(1.d0/4.d0)
*
      return
      end
***
      real*8 FUNCTION mbagbf(mc)
      implicit none
      real*8 mc,mc4,a(200)
      common /GBCFF/ a
*
* A function to evaluate mass at the BAGB by inverting mcagbf.
*
      mc4 = mc**4
      if(mc4.gt.a(37))then
         mbagbf = ((mc4 - a(37))/a(35))**(1.d0/a(36))
      else
         mbagbf = 0.d0
      endif
*
      return
      end
***
      real*8 FUNCTION mcgbtf(t,A,GB,tinf1,tinf2,tx)
      implicit none
      real*8 t,A,GB(10),tinf1,tinf2,tx
*
* A function to evaluate Mc given t for GB, AGB and NHe stars
*
      if(t.le.tx)then
         mcgbtf = ((GB(5)-1.d0)*A*GB(4)*(tinf1 - t))**
     &                              (1.d0/(1.d0-GB(5)))
      else
         mcgbtf = ((GB(6)-1.d0)*A*GB(3)*(tinf2 - t))**
     &                              (1.d0/(1.d0-GB(6)))
      endif
*
      return
      end
***
      real*8 FUNCTION lgbtf(t,A,GB,tinf1,tinf2,tx)
      implicit none
      real*8 t,A,GB(10),tinf1,tinf2,tx
*
* A function to evaluate L given t for GB, AGB and NHe stars
*
      if(t.le.tx)then
         lgbtf = GB(4)*(((GB(5)-1.d0)*A*GB(4)*(tinf1 - t))**
     &                              (GB(5)/(1.d0-GB(5))))
      else
         lgbtf = GB(3)*(((GB(6)-1.d0)*A*GB(3)*(tinf2 - t))**
     &                              (GB(6)/(1.d0-GB(6))))
      endif
*
      return
      end
***
      real*8 FUNCTION mcgbf(lum,GB,lx)
      implicit none
      real*8 lum,GB(10),lx
*
* A function to evaluate Mc given L for GB, AGB and NHe stars
*
      if(lum.le.lx)then
         mcgbf = (lum/GB(4))**(1.d0/GB(5))
      else
         mcgbf = (lum/GB(3))**(1.d0/GB(6))
      endif
*
      return
      end
***
      real*8 FUNCTION lmcgbf(mc,GB)
      implicit none
      real*8 mc,GB(10)
*
* A function to evaluate L given Mc for GB, AGB and NHe stars
*
      if(mc.le.GB(7))then
         lmcgbf = GB(4)*(mc**GB(5))
      else
         lmcgbf = GB(3)*(mc**GB(6))
      endif
*
      return
      end
***
      real*8 FUNCTION lHeIf(m,mhefl)
      implicit none
      real*8 m,mhefl,a(200)
      common /GBCFF/ a
*
* A function to evaluate He-ignition luminosity  (OP 24/11/97)
* Continuity between the LM and IM functions is ensured with a first
* call setting lhefl = lHeIf(mhefl,0.0)
*
      if(m.lt.mhefl)then
         lHeIf = a(38)*m**a(39)/(1.d0 + a(41)*EXP(m*a(40)))
      else
         lHeIf = (a(42) + a(43)*m**3.8d0)/(a(44) + m**2)
      endif
*
      return
      end
***
      real*8 FUNCTION lHef(m)
      implicit none
      real*8 m,a(200)
      common /GBCFF/ a
*
* A function to evaluate the ratio LHe,min/LHeI  (OP 20/11/97)
* Note that this function is everywhere <= 1, and is only valid
* for IM stars
*
      lHef = (a(45) + a(46)*m**(a(48)+0.1d0))/(a(47) + m**a(48))
*
      return
      end
***
      real*8 FUNCTION rminf(m)
      implicit none
      real*8 m,mx,a(200)
      common /GBCFF/ a
*
* A function to evaluate the minimum radius during He-burning
* for IM & HM stars  (OP 20/11/97)
*
      mx = m**a(53)
      rminf = (a(49)*m + (a(50)*m)**a(52)*mx)/(a(51) + mx)
*
      return
      end
***
      real*8 FUNCTION tHef(m,mc,mhefl)
      implicit none
      real*8 m,mc,mhefl,mm,a(200)
      common /GBCFF/ a
      real*8 themsf
      external themsf
*
* A function to evaluate the He-burning lifetime.  (OP 26/11/97)
* For IM & HM stars, tHef is relative to tBGB.
* Continuity between LM and IM stars is ensured by setting
* thefl = tHef(mhefl,0.0,,0.0), and the call to themsf ensures
* continuity between HB and NHe stars as Menv -> 0.
*
      if(m.le.mhefl)then
         mm = MAX((mhefl - m)/(mhefl - mc),1.0d-12)
         tHef = (a(54) + (themsf(mc) - a(54))*mm**a(55))*
     &          (1.d0 + a(57)*EXP(m*a(56)))
      else
         mm = m**5
         tHef = (a(58)*m**a(61) + a(59)*mm)/(a(60) + mm)
      endif
*
      return
      end
***
      real*8 FUNCTION tblf(m,mhefl,mfgb)
      implicit none
      real*8 m,mhefl,mfgb,mr,m1,m2,r1,a(200)
      common /GBCFF/ a
      real*8 lheif,rminf,ragbf
      external lheif,rminf,ragbf
*
* A function to evaluate the blue-loop fraction of the He-burning
* lifetime for IM & HM stars  (OP 28/01/98)
*
      mr = mhefl/mfgb
      if(m.le.mfgb) then
         m1 = m/mfgb
         m2 = log10(m1)/log10(mr)
         m2 = max(m2,1.0d-12)
         tblf = a(64)*m1**a(63) + a(65)*m2**a(62)
      else
         r1 = 1.d0 - rminf(m)/ragbf(m,lheif(m,mhefl),mhefl)
         r1 = max(r1,1.0d-12)
         tblf = a(66)*m**a(67)*r1**a(68)
      end if
      tblf = MIN(1.d0,MAX(0.d0,tblf))
      if(tblf.lt.1.0d-10) tblf = 0.d0
*
      return
      end
***
      real*8 FUNCTION lzahbf(m,mc,mhefl)
      implicit none
      real*8 m,mc,mhefl,mm,a4,a5,a(200)
      common /GBCFF/ a
      real*8 lzhef
      external lzhef
*
* A function to evaluate the ZAHB luminosity for LM stars. (OP 28/01/98)
* Continuity with LHe,min for IM stars is ensured by setting
* lx = lHeif(mhefl,z,0.0,1.0)*lHef(mhefl,z,mfgb), and the call to lzhef
* ensures continuity between the ZAHB and the NHe-ZAMS as Menv -> 0.
*
      a5 = lzhef(mc)
      a4 = (a(69) + a5 - a(74))/((a(74) - a5)*exp(a(71)*mhefl))
      mm = MAX((m-mc)/(mhefl - mc),1.0d-12)
      lzahbf = a5 + (1.d0 + a(72))*a(69)*mm**a(70)/
     &         ((1.d0 + a(72)*mm**a(73))*(1.d0 + a4*EXP(m*a(71))))
*
      return
      end
***
      real*8 FUNCTION rzahbf(m,mc,mhefl)
      implicit none
      real*8 m,mc,mhefl,rx,ry,mm,f,a(200)
      common /GBCFF/ a
      real*8 rzhef,rgbf,lzahbf
*
* A function to evaluate the ZAHB radius for LM stars. (OP 28/01/98)
* Continuity with R(LHe,min) for IM stars is ensured by setting
* lx = lHeif(mhefl,z,0.0,1.0)*lHef(mhefl,z,mfgb), and the call to rzhef
* ensures continuity between the ZAHB and the NHe-ZAMS as Menv -> 0.
*
      rx = rzhef(mc)
      ry = rgbf(m,lzahbf(m,mc,mhefl))
      mm = MAX((m-mc)/(mhefl - mc),1.0d-12)
      f = (1.d0 + a(76))*mm**a(75)/(1.d0 + a(76)*mm**a(77))
      rzahbf = (1.d0 - f)*rx + f*ry
*
      return
      end
***
      real*8 FUNCTION lzhef(m)
      implicit none
      real*8 m,m15
*
* A function to evaluate Naked Helium star 'ZAMS' luminosity
*
      m15 = m*SQRT(m)
      lzhef = 1.5262d+04*m**(41.d0/4.d0)/
     &        (0.0469d0 + m**6*(31.18d0 + m15*(29.54d0 + m15)))
*
      return
      end
***
      real*8 FUNCTION rzhef(m)
      implicit none
      real*8 m
*
* A function to evaluate Helium star 'ZAMS' radius
*
      rzhef = 0.2391d0*m**4.6d0/(0.0065d0 + (0.162d0 + m)*m**3)
*
      return
      end
***
      real*8 FUNCTION themsf(m)
      implicit none
      real*8 m
*
* A function to evaluate Helium star main sequence lifetime
*
      themsf = (0.4129d0 + 18.81d0*m**4 + 1.853d0*m**6)/m**(13.d0/2.d0)
*
      return
      end
***
      real*8 FUNCTION rhehgf(m,lum,rx,lx)
      implicit none
      real*8 m,lum,rx,lx,cm
*
* A function to evaluate Helium star radius on the Hertzsprung gap 
* from its mass and luminosity. 
*
      cm = 2.0d-03*m**(5.d0/2.d0)/(2.d0 + m**5)
      rhehgf = rx*(lum/lx)**0.2d0 + 0.02d0*(EXP(cm*lum) - EXP(cm*lx))
*
      return
      end
***
      real*8 FUNCTION rhegbf(lum)
      implicit none
      real*8 lum
*
* A function to evaluate Helium star radius on the giant branch. 
*
      rhegbf = 0.08d0*lum**(3.d0/4.d0)
*
      return
      end
***
      real*8 FUNCTION lpertf(m,mew)
      implicit none
      real*8 m,mew
      real*8 b,c
*
* A function to obtain the exponent that perturbs luminosity.
*
      b = 0.002d0*MAX(1.d0,2.5d0/m)
      c = 3.d0
      lpertf = ((1.d0 + b**c)*((mew/b)**c))/(1.d0+(mew/b)**c)
*
      return
      end
***
      real*8 FUNCTION rpertf(m,mew,r,rc)
      implicit none
      real*8 m,mew,r,rc
      real*8 a,b,c,q,fac,facmax
*
* A function to obtain the exponent that perturbs radius.
*
      if(mew.le.0.d0)then
         rpertf = 0.d0
      else
         a = 0.1d0
         b = 0.006d0*MAX(1.d0,2.5d0/m)
         c = 3.d0
         q = log(r/rc)
         fac = a/q
         facmax = -14.d0/log10(mew)
         fac = MIN(fac,facmax)
         rpertf = ((1.d0 + b**c)*((mew/b)**c)*(mew**fac))/
     &            (1.d0+(mew/b)**c)
      endif
*
      return
      end
***
      real*8 FUNCTION vrotf(m)
      implicit none
      real*8 m
*
      vrotf = 330.d0*m**3.3d0/(15.d0 + m**3.45d0)
*
      return
      end
***
      real*8 FUNCTION celamf(kw,m,lum,rad,rzams,menv,fac)
      implicit none
      integer kw
      real*8 m,lum,rad,rzams,menv,fac
      real*8 lam1,lam2,m1,logm,logl
      real*8 aa,bb,cc,dd
*
* A function to estimate lambda for common-envelope.
*
      if(fac.ge.0.d0)then
*
* No fits yet for naked He stars...
*
         if(kw.gt.6)then
            celamf = 0.5d0
            goto 90
         endif
*
         if(menv.gt.0.d0)then
* Formulae for giant-like stars; also used for HG and CHeB stars close
* to the Hayashi track.
            logl = log10(lum)
            logm = log10(m)
            if(kw.le.5)then
               m1 = m
               if(kw.gt.3) m1 = 100.d0
               lam1 = 3.d0/(2.4d0 + 1.d0/m1**(3.d0/2.d0)) - 0.15d0*logl
               lam1 = MIN(lam1,0.8d0)
            else
               lam1 = -3.5d0 - 0.75d0*logm + logl
            endif
            if(kw.gt.3)then
               lam2 = MIN(0.9d0,0.58d0 + 0.75d0*logm) - 0.08d0*logl
               if(kw.lt.6)then
                  lam1 = MIN(lam2,lam1)
               else
                  lam1 = MAX(lam2,lam1)
                  lam1 = MIN(lam1,1.d0)
               endif
            endif
            lam1 = 2.d0*lam1
            if(fac.gt.0.d0)then
* Use a fraction FAC of the ionization energy in the energy balance.
               if(kw.le.3)then
                  aa = MIN(1.2d0*(logm - 0.25d0)**2 - 0.7d0,-0.5d0)
               else
                  aa = MAX(-0.2d0 - logm,-0.5d0)
               endif
               bb = MAX(3.d0 - 5.d0*logm,1.5d0)
               cc = MAX(3.7d0 + 1.6d0*logm,3.3d0 + 2.1d0*logm)
               lam2 = aa + ATAN(bb*(cc - logl))
               if(kw.le.3)then
                  dd = MAX(0.d0,MIN(0.15d0,0.15d0 - 0.25d0*logm))
                  lam2 = lam2 + dd*(logl - 2.d0)
               endif
               lam2 = MAX(lam2,1.d-2)
               lam2 = MAX(1.d0/lam2,lam1)
               if(fac.ge.1.d0)then
                  lam1 = lam2
               else
                  lam1 = lam1 + fac*(lam2 - lam1)
               endif
            endif
         endif
*
         if(menv.lt.1.d0)then
* Formula for HG stars; also reasonable for CHeB stars in blue loop.
            lam2 = 0.42d0*(rzams/rad)**0.4d0
* Alternatively:
*           lam2 = 0.3d0*(rtms/rad)**0.4d0
            lam2 = 2.d0*lam2
         endif
*
         if(menv.le.0.d0)then
            celamf = lam2
         elseif(menv.ge.1.d0)then
            celamf = lam1
         else
* Interpolate between HG and GB values depending on conv. envelope mass.
            celamf = lam2 + sqrt(menv)*(lam1 - lam2)
         endif
*
 90      continue
*
      else
         celamf = -1.d0*fac
      endif
*
      return
      end
***

***
      PROGRAM sse
c-------------------------------------------------------------c
c
c     Evolves a single star.
c     Mass loss is an option.
c     The timestep is not constant but determined by certain criteria.
c
c     Written by Jarrod Hurley 26/08/97 at the Institute of
c     Astronomy, Cambridge.
c
c-------------------------------------------------------------c
c
c     STELLAR TYPES - KW
c
c        0 - deeply or fully convective low mass MS star
c        1 - Main Sequence star
c        2 - Hertzsprung Gap
c        3 - First Giant Branch
c        4 - Core Helium Burning
c        5 - First Asymptotic Giant Branch
c        6 - Second Asymptotic Giant Branch
c        7 - Main Sequence Naked Helium star
c        8 - Hertzsprung Gap Naked Helium star
c        9 - Giant Branch Naked Helium star
c       10 - Helium White Dwarf
c       11 - Carbon/Oxygen White Dwarf
c       12 - Oxygen/Neon White Dwarf
c       13 - Neutron Star
c       14 - Black Hole
c       15 - Massless Supernova
c
c-------------------------------------------------------------c
      implicit none
*
*
* const_bse.h
*
      INTEGER idum
      COMMON /VALUE3/ idum
      INTEGER idum2,iy,ir(32)
      COMMON /RAND3/ idum2,iy,ir
      INTEGER ktype(0:14,0:14)
      COMMON /TYPES/ ktype
      INTEGER ceflag,tflag,ifflag,nsflag,wdflag
      COMMON /FLAGS/ ceflag,tflag,ifflag,nsflag,wdflag
      INTEGER bhflag
*
      REAL*8 neta,bwind,hewind,mxns,alpha1,lambda
      REAL*8 sigma,beta,xi,acc2,epsnov,eddfac,gamma
      COMMON /VALUE1/ neta,bwind,hewind,mxns
      COMMON /VALUE2/ alpha1,lambda
      COMMON /VALUE4/ sigma,bhflag
      COMMON /VALUE5/ beta,xi,acc2,epsnov,eddfac,gamma
      REAL*8 pts1,pts2,pts3
      COMMON /POINTS/ pts1,pts2,pts3
      REAL*8 dmmax,drmax
      COMMON /TSTEPC/ dmmax,drmax
      REAL scm(50000,14),spp(20,3)
      COMMON /SINGLE/ scm,spp
      REAL bcm(50000,34),bpp(80,10)
      COMMON /BINARY/ bcm,bpp
*

*
      integer kw,j,k
*
      real*8 mass,mt,z,zpars(20)
      real*8 epoch,tms,tphys,tphysf,dtp
      real*8 r,lum,ospin
      real*8 mc,rc,menv,renv
*
************************************************************************
* Input:
*
* mass is in solar units.
* z is metallicity in the range 0.0001 -> 0.03 where 0.02 is Population I.
* tphysf is the maximum evolution time in Myr.
*
* neta is the Reimers mass-loss coefficent (neta*4x10^-13; 0.5 normally).
* bwind is the binary enhanced mass loss parameter (inactive for single).
* hewind is a helium star mass loss factor (1.0 normally). 
* sigma is the dispersion in the Maxwellian for the SN kick speed (190 km/s). 
*
* ifflag > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800 (0). 
* wdflag > 0 uses modified-Mestel cooling for WDs (0). 
* bhflag > 0 allows velocity kick at BH formation (0). 
* nsflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407 (1). 
* mxns is the maximum NS mass (1.8, nsflag=0; 3.0, nsflag=1). 
* idum is the random number seed used in the kick routine. 
*
* Next come the parameters that determine the timesteps chosen in each
* evolution phase:
*                 pts1 - MS                  (0.05)
*                 pts2 - GB, CHeB, AGB, HeGB (0.01)
*                 pts3 - HG, HeMS            (0.02)
* as decimal fractions of the time taken in that phase.
*
* If you enter a negative mass then parameters for an evolved star are
* required in the order of:
* initial mass, current mass, type, current time & epoch,
* otherwise the star will start on the ZAMS.
*
      READ(*,*)mass,z,tphysf
      READ(*,*)neta,bwind,hewind,sigma
      READ(*,*)ifflag,wdflag,bhflag,nsflag,mxns
      READ(*,*)pts1,pts2,pts3
*
************************************************************************
*
* Set parameters which depend on the metallicity 
*
      CALL zcnsts(z,zpars)
      if(idum.gt.0) idum = -idum
*
      if(mass.gt.0.0)then
*
* Initialize the star to begin on the ZAMS.
*
         mt = mass
         kw = 1
         tphys = 0.d0
         epoch = 0.d0
      else
*
* Obtain parameters of evolved star (from stdin).
*
         READ(*,*)mass,mt,kw,tphys,epoch
      endif
*
* Set the initial spin of the star.
*
      ospin = 0.d0
*
* Set the data-save parameter.
*
      dtp = 0.d0
* 
      CALL evolv1(kw,mass,mt,r,lum,mc,rc,menv,renv,ospin,
     &            epoch,tms,tphys,tphysf,dtp,z,zpars)
*
************************************************************************
* Output:
*
      j = 0
      if(scm(1,1).lt.0.0) goto 50
*
* The scm array stores the stellar parameters at the specified output 
* times. The parameters are (in order of storage):
*    
*    Time, stellar type, initial mass, current mass, log10(L), log10(r),
*    log10(Teff), core mass, epoch and spin.
*
 30   j = j + 1
      if(scm(j,1).lt.0.0)then
         scm(j-1,1) = scm(j,1)
         j = j - 1
      endif
      WRITE(*,99)(scm(j,k),k=1,8),scm(j,10),scm(j,12),scm(j,13)
      if(scm(j,1).ge.0.0) goto 30
 99   FORMAT(8f10.4,1p,e12.4,0p,f12.4,1p,e12.4)
*
 50   continue
*
************************************************************************
*
      STOP
      END
***