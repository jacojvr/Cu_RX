! ! ! ! compile using F2Py:
! ! ! ! f2py --fcompiler=gfortran -c ./umat_1d.f -m pythonrx
! !
      subroutine umat1d(stress, statev, dstran,
     & dtime, temp, nstatv, props, nprops)
!f2py integer intent(hide),depend(statev) :: nstatv = len(statev)
!f2py real*8 intent(in,out),dimension(nstatv) :: statev
!f2py integer intent(hide),depend(props) :: nprops = len(props)
!f2py real*8 intent(in),dimension(nprops) :: props
!f2py real*8 intent(in) :: dstran,dtime,temp
!f2py real*8 intent(in,out) :: stress
c
      implicit real*8(a-h,o-z)
      dimension statev(nstatv),props(nprops),tempstatev(nstatv)
      parameter(zero=0.d0, one=1.d0, two=2.d0, three=3.d0, six=6.d0,
     & enumax=.4999d0, newton=10, toler=1.0d-8)
c
c elastic properties
c
      emu0=props(1) ! mu0 at temperature T0
      ed0=props(2) ! D0 at temperature T0
      et0=props(3) ! Temperature T0
      enu=min(dabs(props(4)),enumax) 
!  
      if(temp.gt.zero) then
      eg = emu0 - ed0/(dexp(et0/temp)-one)
      else
      eg = emu0 ! use reference value if no thermal effects set
      endif
      eg2=two*eg
      eg3=three*eg
      emod = eg2*(one+enu)
c
c initialize temporary statevariable array
      do i =1,nstatv
      tempstatev(i) = statev(i)
      end do
c
c calculate predictor stress and elastic strain
c
      stress=stress+emod*dstran
c
c calculate trial
c
      strial=dabs(stress)
c
c check yield condition
      call tdbrec(syiel0,dsy,statev,tempstatev,nstatv,
     & zero,dtime,temp,props,nprops)
!      write(*,*) syiel0,' ',strial
      
c
c determine if actively yielding
c
      if (strial.gt.(one+toler)*syiel0) then
!      write(*,*) "F: plasticity"
c
c actively yielding
c separate the hydrostatic from the deviatoric stress
c calculate the flow direction
c
      flow=stress/strial
c
c solve for equivalent von mises stress
c and equivalent plastic strain increment using newton iteration
c
      syield=syiel0
      cop=toler
      do kewton=1, newton
      rhs=strial-emod*cop-syield
      cop=dabs(cop+rhs/(emod+dsy))
      deqpl = cop
!      call tdbrec(syieldn,dsy,statev,tempstatev,
!      & nstatv,deqpl+toler,dtime,temp,props,nprops)
      call tdbrec(syield,dsy,statev,tempstatev,
     & nstatv,deqpl,dtime,temp,props,nprops)
!      write(*,*) ' '
!      write(*,*) 'FD: ',(syieldn-syield)/toler
!      write(*,*) 'AN: ',dsy
      if(abs(rhs).lt.toler*syiel0) goto 10
      end do
c
c write warning message to .msg file
c
!      write(7,2) newton
!  2    format(//,30x,'***warning - plasticity algorithm did not ',
!      & 'converge after ',i3,' iterations')
 10   continue
c
c update stress, elastic and plastic strains and
c equivalent plastic strain
      stress=flow*syield!*c2
      
!      else if((strial.lt.toler).and.(statev(11).gt.zero))then
!      !! shift ISV's
!      do k1=1,nstatv
!      statev(k1) = tempstatev(k1)
!      tempstatev(k1) = zero
!      end do
!      do k1=1,4
!      tempstatev(k1) = statev(k1)
!      tempstatev(k1+4) = statev(k1)
!      end do
      
      endif
      
      do k1=1,nstatv
      statev(k1) = tempstatev(k1)
      end do
c
      return
      end
c

      subroutine getR(R,dRdxj,fxnvec,Reps,xj,xi,fxinfo,
     &  depl,dtime,temp,props,nprops)
c
      implicit real*8(a-h,o-z)
      dimension props(nprops),R(2),dRdxj(2,2),fxinfo(3),Reps(2),
     &  xj(2),xi(2),fxnvec(4)
c
      parameter(zero=0.d0,half=0.5d0,one=1.d0,two=2.d0,
     & toler=1.d-10,ratelim=1.d-8)
c
c Elastic Properties:
      emu0 = props(1)
      ed0 = props(2) 
      et0 = props(3) 
      enu = props(4)
c Reference stress values
      siga = props(5)
      sig0 = props(6)
c Scaling function
      a0e = props(28)
      rate0 = props(29)
      !qe = props(30)
      !pe = props(31)
c Evolution of misorient:      
      cld = one!props(10)
c stage 4
      cg = props(19)
      rg = props(20)
c Storage
      c1 = props(12)
c Dynamic recovery:
      c20 = props(13)
      a02 = props(14)
      rate02 = props(15)
c Thermal recovery
      c30 = props(16)
      r3 = props(17)
      a03 = props(18)
c Recrystallisation
      cx0 = props(21)
      a0x = props(22)
      cxl = props(23)
      rxl = props(24)
      rxa = props(25)
      rxb = props(26)
      cxc = props(27)
c
c get info from previous values: contained in fxinfo:
      fxc = fxinfo(1)
      fxcr = dabs(fxinfo(2))
      fxnp = fxinfo(3)
c
      rate = depl/dtime
      if(rate.lt.ratelim)then
      rate = ratelim
      endif
c
      if(temp.gt.et0)then
      emu = emu0 - ed0/(dexp(et0/temp)-one)
      r2m = -temp/(emu*a02)
      r3c = c30*dexp(-a03/temp)*dtime
      r5c0 = cx0*dexp(-a0x/temp)*emu*dtime
      else
c Use constants:
      emu = emu0
      r2m = -a02
      r3c = c30*dtime
      r5c0 = cx0*emu*dtime
      endif
c ISV's at previous convergence and current guess:
      x1p = max(xi(1),one)
      x2p = max(xi(2),zero)
      x1 = max(xj(1),one)
      x2 = max(xj(2),zero)
c
c Growth of next recrystallised volume:
      cldbar = min(x2,one)
      r5c1 = (one-dexp(-cxl*(cldbar)**rxl))
      r5c = r5c0*r5c1*(x1)
c interfacial area 
      fxn = fxnp
      rxg = fxc*((fxn/fxc)**rxa)*((one-fxn/fxc)**rxb)*
     & (one+cxc*(one-fxc))
      drxg = rxa*fxc*((fxn/fxc)**(rxa-one))*
     & ((one-fxn/fxc)**rxb)*(one+cxc*(one-fxc)) -
     & rxb*fxc*((fxn/fxc)**rxa)*((one-fxn/fxc)**(rxb-one))*
     & (one+cxc*(one-fxc))
      fxnr = dabs(r5c*rxg)
      ffxn = fxn - fxnp - fxnr
c resolve residual
      icount = 0
      do while((icount.lt.15).and.(dabs(ffxn).gt.toler))
      icount = icount+1
      dffxn = one - half*r5c*drxg
      if(dabs(dffxn).lt.toler)then
      dffxn = toler
      endif
      fxn = min(dabs(fxn - ffxn/dffxn),fxc-toler)
      rxg = fxc*((fxn/fxc)**rxa)*((one-fxn/fxc)**rxb)*
     & (one+cxc*(one-fxc))
      drxg = rxa*fxc*((fxn/fxc)**(rxa-one))*
     & ((one-fxn/fxc)**rxb)*(one+cxc*(one-fxc)) -
     & rxb*fxc*((fxn/fxc)**rxa)*((one-fxn/fxc)**(rxb-one))*
     & (one+cxc*(one-fxc))
      fxnr = dabs(r5c*rxg)
      ffxn = fxn - fxnp - fxnr
      end do
c ----------------------------------------------------------------
c PARTIAL : CHANGE OF fxn WITH RESPECT TO x1 AND x2:
c partial gradients d(fxn)/d(x1)   
      ddfxdr = (one-r5c*drxg)
      if(dabs(ddfxdr).lt.toler)then
      ddfxdr = toler
      endif 
      dfxndx1 = half*(r5c0*r5c1*rxg)/ddfxdr
      dmdx2 = zero
      if(x2.lt.one)then
      dmdx2 = one
      endif
      dr5c1dm = rxl*cxl*dexp(-cxl*(cldbar)**rxl)*
     & (cldbar)**(rxl-one)
      dfxndm = half*(r5c0*dr5c1dm*x1*rxg)/ddfxdr
c partial gradients d(fxn)/d(x2)
      dfxndx2 = dfxndm*dmdx2
c ----------------------------------------------------------------
c Residual equations on the ISV values:
      c2 = c20*(rate/rate02)**r2m
      dc2de = c20*r2m*((rate/rate0)**(r2m-one))/
     & (dtime*rate0)
      sqx1p = dsqrt(x1p)
      sqx1 = dsqrt(x1)
      hx2 = cld
      hx1TEMP = -r3c*(x1**r3+x1p**r3)
      hx1 = (cg)*(x2)**rg+c1*sqx1-c2*x1
      if(fxn.ge.fxc)then
      rx0 = zero
      else
      rx0 = one/(fxc - fxn)
      endif
      if(x1p.eq.one)then
      rx0=zero
      hx2 = zero
      endif
      rfxc = fxcr*rx0
      drx = rfxc*rx0
      R1 = x1-x1p-hx1*depl-hx1TEMP+x1*rfxc
      R2 = x2-x2p-hx2*depl+x2*rfxc
      R = (/ R1 , R2 /)
c ----------------------------------------------------------------
c PARTIAL : CHANGE IN RESIDUAL WITH RESPECT TO x1 AND x2:
      dhrdrTEMP=-r3*r3c*x1**(r3-one)
      dhrdr=half*c1/sqx1-c2
c  # partial gradients d(fR1)/d(x1)
      dR1dr0=one-dhrdr*depl-dhrdrTEMP+rfxc
      dR1dx1=dR1dr0+x1*drx*dfxndx1
c  # partial gradients d(fR1)/d(x2)
      dR1dx2=-depl*rg*(cg)*(x2)**(rg-one)+
     & x1*drx*dfxndx2
c partial gradients d(F2)/d(x1)
      dR2dx1 = x2*drx*dfxndx1
c partial gradients d(F2)/d(x2)
      dR2dl0= one + rfxc
      dR2dx2 = dR2dl0 + x2*drx*dfxndx2
      dRdxj=reshape((/dR1dx1,dR2dx1,dR1dx2,dR2dx2/),(/2,2/))
c
      dR1de = hx1 - depl*dc2de*x1
      dR2de = hx2
      Reps = (/dR1de,dR2de/)  
c
c exchange new supplementary info using fxinfo:
      fxnvec(1) = fxn
      fxnvec(2) = fxnr
      fxnvec(3) = dfxndx1
      fxnvec(4) = dfxndx2
c      
      return
      end
c
c
c
c ----------------------------------------------------------------
      subroutine tdbrec(sy,dsy,statev,tempstatev,
     & nstatv,depl,dtime,temp,props,nprops)
!      
! c ----------------------------------------------------------------
!      subroutine fisotropic(sy,dsy,depl,dtime,temp,
!      & statev,tempstatev,nstatv,props,nprops)
c
      implicit real*8(a-h,o-z)
      logical checkrx
      dimension props(nprops),statev(nstatv),tempstatev(nstatv),
     & xi(2),xj(2),R(2),dRdX(2,2),fxinfo(3),Reps(2),dXdR(2,2),
     & fxnvec(4),xjupd(2)
     
      parameter(zero=0.d0,half=0.5d0,one=1.d0,two=2.d0,
     & toler=1.d-4,x10=one,x20=1.d-10,fxn0=1.d-3,ratelim=1.d-8)
c
c Elastic Properties:
      emu0 = props(1)
      ed0 = props(2) 
      et0 = props(3) 
      enu = props(4)
c Reference stress values
      siga = props(5)
      sig0 = props(6)
c Scaling function
      a0e = props(28)
      rate0 = props(29)
      qe = props(30)
      pe = props(31)
c
      rate = depl/dtime
      if(rate.lt.ratelim)then
      rate = ratelim
      endif
c
      if(temp.gt.et0)then
      emu = emu0 - ed0/(dexp(et0/temp)-one)
      sfe0 = temp/(a0e*emu)
      else
      emu = emu0
      sfe0 = one/a0e
      endif
      emusf = emu/emu0
c 
      sfel = dlog(rate0/rate)*sfe0
      sfe = dabs(one-sfel**(one/qe))**(one/pe)
c ISV shift
      nrrx = (nstatv-3)/4 
      if(statev(7).gt.(0.999d0))then
      ixvf=0
      do while(ixvf.lt.nrrx)
      lstskip = 4*(ixvf)+3
      statev(lstskip+1)=statev(lstskip+5)
      statev(lstskip+2)=statev(lstskip+6)
      statev(lstskip+3)=statev(lstskip+7)
      statev(lstskip+4)=statev(lstskip+8)
      enddo
      statev(lstskip+5)=zero
      statev(lstskip+6)=zero
      statev(lstskip+7)=zero
      statev(lstskip+8)=zero
      endif
c ----------------------------------------------------------------      
      fxc = one
      fxcp = one
      fxcr = zero
      dfxcde = zero
      x1eq = zero
      x2eq = zero
      dx1de = zero
      plastic = zero
c
      ixvf = 1
      checkrx = .True.
      do while((ixvf.lt.nrrx).and.(checkrx))
      lstskip = 4*(ixvf-1)+3
      xeplp = statev(lstskip+1)
      x1p = max(statev(lstskip+2),x10)
      x2p = max(statev(lstskip+3),x20)
      fxnp = max(statev(lstskip+4),fxn0)
      xi = (/x1p,x2p/)
      xj = (/x1p,x2p/)
      fxinfo = (/fxc,fxcr,fxnp/)
      call getR(R,dRdX,fxnvec,Reps,xj,xi,fxinfo,
     &  depl,dtime,temp,props,nprops)
      fx = dsqrt(R(1)*R(1) + R(2)*R(2))
      fxd = dRdX(1,1)*dRdX(2,2)-dRdX(2,1)*dRdX(1,2)
      icount = 0
      newtmax=15
      if(xi(1).eq.(one))then
      newtmax = 50
      endif
      do while((icount.lt.newtmax).and.(dabs(fx).ge.toler))
      icount = icount+1
      if(dabs(fxd).gt.zero)then
      dXdR = reshape((/dRdX(2,2),-dRdX(2,1),
     &      -dRdX(1,2),dRdX(1,1)/),(/2,2/))/fxd
      xjupd = reshape(matmul(dXdR,reshape(R,(/2,1/))),(/2/))
      xj=xj-xjupd
      xj = (/max(dabs(xj(1)),x10),max(dabs(xj(2)),x20)/)
      fxinfo = (/fxc,fxcr,fxnp/)
      call getR(R,dRdX,fxnvec,Reps,xj,xi,fxinfo,
     &  depl,dtime,temp,props,nprops)
      fx = dsqrt(R(1)*R(1) + R(2)*R(2))
      fxd = dRdX(1,1)*dRdX(2,2)-dRdX(2,1)*dRdX(1,2)
      else
      xj = (/x1p,x2p/)
      fx = zero
      endif
      enddo
      x1 = max(xj(1),x10)
      x2 = max(xj(2),x20)
      fxn = min(dabs(fxnvec(1)),one)
      if(fxn.le.(two*fxn0))then
      checkrx = .False.
      fxn = zero
      endif
      fxnr = fxnvec(2)
      dfxndx1 = fxnvec(3)
      dfxndx2 = fxnvec(4)
      xepl = xeplp*fxcp/fxc+depl
      if ((fxc-fxn).gt.(two*fxn0))then
      tempstatev(lstskip+1) = xepl
      tempstatev(lstskip+2) = x1
      tempstatev(lstskip+3) = x2
      tempstatev(lstskip+4) = fxn
c
      x1eq = x1eq + x1*(fxc-fxn)
      plastic = plastic + xepl*(fxc-fxn)
      x2eq = x2eq + x2*(fxc-fxn)
      dxeqI = dxeqI + (fxc-fxn)/dxc
      dxeqI2 = dxeqI2 + (fxc-fxn)/dxc2
      
      if(dabs(fxd).gt.0)then
      dx1de = dXdR(1,1)*Reps(1)+dXdR(1,2)*Reps(2)
      dx2de = dXdR(2,1)*Reps(1)+dXdR(2,2)*Reps(2)
      dfxnde = dfxndx1*dx1de+dfxndx2*dx2de
      dx1de = dx1de + dx1de*(fxc-fxn) +
     &      x1*(dfxcde - dfxnde)
      dfxcde = dfxnde
      fxc = fxn
      fxcp = fxnp
      fxcr = fxnr
      endif
      endif
      ixvf = ixvf+1
      end do
      
      
      tempstatev(1) = plastic
      tempstatev(2) = x1eq
      tempstatev(3) = x2eq
c
      sqx1 = dsqrt(x1eq)
      sec = sig0*sqx1
      sy = siga + emusf*sfe*sec
c partial derivatives
c     d(sec)/d(epl)
      dsecdepl = half*sig0*dx1de/sqx1
c     d(sfe)/d(epl)
      dsfedepl = (sfe0*(one-sfel**(one/qe))**(one/pe-one)*
     &     sfel**(one/qe-one)/(pe*qe*rate))/dtime
c     Total
      dsy = emusf*(sfe*dsecdepl+dsfedepl*sec)   
      return
      end
      