      subroutine scawf2(point, npoint, dmat, nbasis, eta, ngmx, nfirst,
     &                  nlast, ntype, vlist, ncmx)
C
C     GENERATE THE ELECTROSTATIC POTENTIAL DIRECTLY FROM QM WAVE FUNCTION
C
C     INPUT: POINT - Flattened array of size (npoint x 8). 3 first columns
C                    contain points (x,y,z). 5-th column is filled with 
C                    potential values
C
C
      implicit double precision (a-h,o-z) 
      double precision point(npoint*8)
      integer npoint, ngmx, ncmx 
      integer nfirst(nbasis), nlast(nbasis), ntype(nbasis)
      double precision eta(ngmx,5), vlist(ncmx,4) 
      double precision dmat(nbasis*nbasis)

      double precision Airu(10), Ajsv(10), Aktw(10)                    
      double precision p(3), sf(10,3), tf(20)
      double precision fact(20), g(50)
      double precision gtoC(1000), dfact(9)
      integer nr(20,3)
      data zero,one,two,half,quart /0.0d00,1.0d00,2.0d00,0.5d00,0.25d00/
      data pi,pitern,onep5/3.141592653589d00,5.568327997d00,1.5d00/
      data fact/1.0D00,1.0D00,2.0D00,6.0D00,24.0D00,120.0D00,
     &          720.0D00,5040.0D00,40320.0D00,362880.0D00,3628800.0D00,
     &          39916800.0D00,479001600.0D00,6227020800.0D00,6*0.0D00/
      data dfact/1.0D00,1.0D00,3.0D00,15.0D00,105.0D00,945.0D00,
     &           10395.0D00,135135.0D00,2027025.0D00/
      data nr /
     &     0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,1,0,1,0,1,
     &     0,0,1,0,0,2,0,1,0,1,0,3,0,1,0,2,2,0,1,1,
     &     0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,0,1,2,2,1/
Cf2py INTENT(In,OUT) point
C
C     NORMALIZE THE PRIMITIVES
C
      do i=1,ngmx
         gtoC(i) = eta(i,5)
      end do
      do j = 1, nbasis
         jtyp = ntype(j); js = nfirst(j); jf = nlast(j)
         l = nr(jtyp,1); m = nr(jtyp,2); n = nr(jtyp,3)
c         write(*,*) j, l, m, n
         do i = js, jf
            alpha = eta(i,4); SOO = pitern*(half/alpha)**onep5
            t1 = dfact(l+1)/alpha**l
            t2 = dfact(m+1)/alpha**m
            t3 = dfact(n+1)/alpha**n
            eta(i,5) = one/dsqrt(SOO*t1*t2*t3)
         end do
      end do
      do j = 1, nbasis
         jtyp = ntype(j); js = nfirst(j); jf = nlast(j)
         l = nr(jtyp,1); m = nr(jtyp,2); n = nr(jtyp,3)
         suma = zero
         do ii = js, jf
            do jj = js, jf
               t = one/(eta(ii,4)+eta(jj,4))
               SOO = pitern*(t**onep5)*eta(ii,5)*eta(jj,5)
               t = half*t
               t1 = dfact(l+1)*t**l
               t2 = dfact(m+1)*t**m
               t3 = dfact(n+1)*t**n
               suma = suma + gtoC(ii)*gtoC(jj)*SOO*t1*t2*t3
            end do
         end do
         suma = one / dsqrt(suma)
         do ii = js, jf
            gtoC(ii) = gtoC(ii) * suma
         end do
      end do
 
      do ii = 1, ngmx
         eta(ii,5) = eta(ii,5) * gtoC(ii)
      end do
C
C     ITERATE OVER BASIS FUNCTIONS
C
      do iii=1,nbasis
      do jjj=1,iii

         ij = nbasis*(iii-1) + jjj
         dmatij = dmat(ij)

         ityp=ntype(iii)  ; jtyp = ntype(jjj)                          
         l1 = nr(ityp,1); m1 = nr(ityp,2); n1 = nr(ityp,3)
         l2 = nr(jtyp,1); m2 = nr(jtyp,2); n2 = nr(jtyp,3)
         imax = l1+l2+1 ; jmax = m1+m2+1 ; kmax = n1+n2+1
         maxall = imax
         if ( maxall.LT.jmax ) maxall = jmax
         if ( maxall.LT.kmax ) maxall = kmax
         if ( maxall.LT.2 )  maxall = 2 
         iss = nfirst(iii); il = nlast(iii)
         jss = nfirst(jjj); jl = nlast(jjj)
                                                                       
         rAB = (eta(iss,1)-eta(jss,1))**2 + 
     &         (eta(iss,2)-eta(jss,2))**2 + 
     &         (eta(iss,3)-eta(jss,3))**2 
                                                                       
         do irun = iss, il                      
            do jrun = jss, jl                      
               aexp = eta(irun,4); anorm = eta(irun,5) 
               bexp = eta(jrun,4); bnorm = eta(jrun,5)
                                                                      
               t1 = aexp + bexp; deleft = one/t1                   
                                                                 
               p(1) = (aexp*eta(irun,1)+bexp*eta(jrun,1))*deleft
               p(2) = (aexp*eta(irun,2)+bexp*eta(jrun,2))*deleft
               p(3) = (aexp*eta(irun,3)+bexp*eta(jrun,3))*deleft
                                                                 
               pax = p(1) - eta(irun,1)
               pay = p(2) - eta(irun,2)
               paz = p(3) - eta(irun,3)
                                                                 
               pbx = p(1) - eta(jrun,1)
               pby = p(2) - eta(jrun,2)
               pbz = p(3) - eta(jrun,3)
                                                                 
               prefa = dexp(-aexp*bexp*rAB/t1)*pi*anorm*bnorm/t1
                                                                   
               prefa = two*prefa                       
               expab = dexp(-aexp*bexp*rAB/t1)         
               s00=(pi/t1)**1.5*expab                  
               dum = one;  tf(1) = one; del = half/t1  
               do n = 2, maxall                        
                  tf(n) = tf(n-1)*dum*del              
                  dum = dum + two                      
               end do                                  
C                                                                         
C              START ACCUMULATE POTENTIAL INTEGRAL
C                                           
               m = imax + jmax + kmax -2                 
               do n = 1, imax
                  sf(n,1) = fj(l1,l2,n-1,pax,pbx)      
               end do
               do n = 1, jmax
                  sf(n,2) = fj(m1,m2,n-1,pay,pby)      
               end do
               do n = 1, kmax
                  sf(n,3) = fj(n1,n2,n-1,paz,pbz)      
               end do
                                                                    
               do iin=1,npoint
                  ipx = 8*(iin-1) + 1
                  ipy = ipx+1
                  ipz = ipy+1
                  ipv = ipz+2
                                                                  
                  xp = point(ipx); yp = point(ipy); zp = point(ipz)
                                                          
                  pn = zero                                           
                  cpx = p(1) - xp
                  cpy = p(2) - yp
                  cpz = p(3) - zp
                  pcsq = cpx*cpx + cpy*cpy + cpz*cpz
                  t = t1*pcsq
                                                                      
                  call auxg(m,t,g)       
                                                                      
                  epsi = quart/t1                         
                  do ii = 1, 10
                     Airu(ii) = zero
                     Ajsv(ii) = zero
                     Aktw(ii) = zero
                  end do
                                                             
                  call aform(imax,sf,fact,cpx,epsi,Airu,1)
                  call aform(jmax,sf,fact,cpy,epsi,Ajsv,2)
                  call aform(kmax,sf,fact,cpz,epsi,Aktw,3)
                                                                      
                  do ii = 1, imax                   
                     do jj = 1, jmax
                        do kk = 1, kmax
                           nu = ii + jj + kk - 2
                           pn = pn + Airu(ii)*Ajsv(jj)*Aktw(kk)*g(nu)
                        end do
                     end do
                  end do
                  vij = prefa*pn*dmatij
                  if (iii.NE.jjj) vij = vij * two
                                                         
                  point(ipv) = point(ipv) - vij
               end do
            end do                 
         end do                

      end do
      end do
C
C     NUCLEAR CONTRIBUTION
C
      do iin=1,npoint
         ipx = 8*(iin-1) + 1
         ipy = ipx+1
         ipz = ipy+1
         ipv = ipz+2
                                                         
         xp = point(ipx); yp = point(ipy); zp = point(ipz)

         v_nuc = zero

         do inuc=1, ncmx
            xi = vlist(inuc,1) - xp
            yi = vlist(inuc,2) - yp
            zi = vlist(inuc,3) - zp
            r  = dsqrt(xi*xi + yi*yi + zi*zi)
            if (r.GT.zero) v_nuc = v_nuc + vlist(inuc,4)/r
         end do
         point(ipv) = point(ipv) + v_nuc
      end do

      return
      end


      subroutine scawfn(point, npoint, dmat, nbasis, eta, ngmx, nfirst,
     &                  nlast, ntype, vlist, ncmx)
cf2py intent(in,out) point
      double precision point(npoint*8)
      double precision dmat(nbasis*nbasis), eta(ngmx,5), vlist(ncmx,4)
      integer nr(20,3), nfirst(nbasis), nlast(nbasis)
      integer ntype(nbasis)
      integer nbasis, ngmx, ncmx, npoint

      double precision v, totnai, oeitot, ovltot, kintot, x, y, z
      double precision xi, yi, zi, r, d, dd, vij, dij
      double precision alpha, pitern, SOO, onep5, half, one, two, zero
      double precision t, t1, t2, t3, suma
      double precision gtoC(1000), dfact(9)
      integer i, j, ii, jj, ij, ipx, ipy, ipz, ipv, ipd, l, m, n
      integer jtyp, js, jf
      data nr /
     &     0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,1,0,1,0,1,
     &     0,0,1,0,0,2,0,1,0,1,0,3,0,1,0,2,2,0,1,1,
     &     0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,0,1,2,2,1/
      data zero,one,two,half,onep5,pitern/0.0D+00,1.0D+00,2.0D+00,
     &                            0.5D+00,1.5D+00,
     &                            5.568327997d00/
      data dfact/1.0D00,1.0D00,3.0D00,15.0D00,105.0D00,945.0D00,
     &           10395.0D00,135135.0D00,2027025.0D00/
C
C     NORMALIZE THE PRIMITIVES
C
      do i=1,ngmx
         gtoC(i) = eta(i,5)
      end do
      do j = 1, nbasis
         jtyp = ntype(j); js = nfirst(j); jf = nlast(j)
         l = nr(jtyp,1); m = nr(jtyp,2); n = nr(jtyp,3)
c         write(*,*) j, l, m, n
         do i = js, jf
            alpha = eta(i,4); SOO = pitern*(half/alpha)**onep5
            t1 = dfact(l+1)/alpha**l
            t2 = dfact(m+1)/alpha**m
            t3 = dfact(n+1)/alpha**n
            eta(i,5) = one/dsqrt(SOO*t1*t2*t3)
         end do
      end do
      do j = 1, nbasis
         jtyp = ntype(j); js = nfirst(j); jf = nlast(j)
         l = nr(jtyp,1); m = nr(jtyp,2); n = nr(jtyp,3)
         suma = zero
         do ii = js, jf
            do jj = js, jf
               t = one/(eta(ii,4)+eta(jj,4))
               SOO = pitern*(t**onep5)*eta(ii,5)*eta(jj,5)
               t = half*t
               t1 = dfact(l+1)*t**l
               t2 = dfact(m+1)*t**m
               t3 = dfact(n+1)*t**n
               suma = suma + gtoC(ii)*gtoC(jj)*SOO*t1*t2*t3
            end do
         end do
         suma = one / dsqrt(suma)
         do ii = js, jf
            gtoC(ii) = gtoC(ii) * suma
         end do
      end do
 
      do ii = 1, ngmx
         eta(ii,5) = eta(ii,5) * gtoC(ii)
      end do
C
C     EVALUATE THE DENSITY AND POTENTIAL AT EACH POINT
C     
      do i=1,npoint
c      do i=1,nx
c      do j=1,ny
c         ipx = ny*8*(i-1) + 8*(j-1) + 1
         ipx = 8*(i-1) + 1
         ipy = ipx + 1
         ipz = ipy + 1
         ipd = ipz + 1
         ipv = ipd + 1
         x = point(ipx); y = point(ipy); z = point(ipz)
         v = zero
         d = zero
c        nuclear contribution
         do ii=1,ncmx
            xi = vlist(ii,1) - x
            yi = vlist(ii,2) - y
            zi = vlist(ii,3) - z
            r  = dsqrt(xi*xi + yi*yi + zi*zi)
            if (r.GT.zero) v = v + vlist(ii,4)/r
         end do
c        electronic contribution
         do ii=1,nbasis
         do jj=1,ii
            ij = nbasis*(ii-1)+jj
            dd = dmat(ij)
            call genoei(ii, jj, x, y, z, eta, ngmx, nfirst, nlast, 
     &                          ntype, nr, 20, vlist, ncmx, ncmx, 
     &                          ovltot, kintot, totnai, oeitot)
c            write(*,*) ii, jj, ovltot
c            write(*,*) totnai
            dij = ovltot*dd
            vij = totnai*dd
            if (ii.NE.jj) then
                dij = dij * two
                vij = vij * two
            end if
            d = d + dij
            v = v + vij
         end do
         end do
         point(ipd) = d
         point(ipv) = v
c         write(*,*) d
c         STOP
      end do
c      end do

      return
      end
C --------------------------------- HELPER SUBROUTINES AND FUNCTIONS ----------------------

      subroutine genoei(i, j, xp, yp, zp, eta, ngmx, nfirst,
     &      nlast, ntype, nr, ntmx, vlist, noc, ncmx, ovltot, kintot,
     &      totnai, oeitot)
      implicit double precision (a-h,o-z) 
      integer i, j, ngmx, ncmx, noc, ntmx
      integer nfirst(*), nlast(*), ntype(*), nr(ntmx,3) 
      double precision ovltot, kintot 
      double precision eta(ngmx,5), vlist(ncmx,4) 

      double precision Airu(10), Ajsv(10), Aktw(10)                    
      double precision p(3), sf(10,3), tf(20)
      double precision fact(20), g(50)
      double precision kin
      data zero,one,two,half,quart /0.0d00,1.0d00,2.0d00,0.5d00,0.25d00/
      data pi/3.141592653589d00/
      data fact/1.0D00,1.0D00,2.0D00,6.0D00,24.0D00,120.0D00,
     &          720.0D00,5040.0D00,40320.0D00,362880.0D00,3628800.0D00,
     &          39916800.0D00,479001600.0D00,6227020800.0D00,6*0.0D00/

      ityp=ntype(i)  ; jtyp = ntype(j)                            
      l1 = nr(ityp,1); m1 = nr(ityp,2); n1 = nr(ityp,3)
      l2 = nr(jtyp,1); m2 = nr(jtyp,2); n2 = nr(jtyp,3)
      imax = l1+l2+1 ; jmax = m1+m2+1 ; kmax = n1+n2+1
      maxall = imax
      if ( maxall.LT.jmax ) maxall = jmax
      if ( maxall.LT.kmax ) maxall = kmax
      if ( maxall.LT.2 )  maxall = 2 
      iss = nfirst(i); il = nlast(i)
      jss = nfirst(j); jl = nlast(j)

      rAB = (eta(iss,1)-eta(jss,1))**2 + 
     &      (eta(iss,2)-eta(jss,2))**2 + 
     &      (eta(iss,3)-eta(jss,3))**2 

      oeitot = zero 
      totnai = zero
      kintot = zero
      ovltot = zero

      do irun = iss, il                      
         do jrun = jss, jl                      
            aexp = eta(irun,4); anorm = eta(irun,5) 
            bexp = eta(jrun,4); bnorm = eta(jrun,5)
     
            t1 = aexp + bexp; deleft = one/t1                   
                                                              
            p(1) = (aexp*eta(irun,1)+bexp*eta(jrun,1))*deleft
            p(2) = (aexp*eta(irun,2)+bexp*eta(jrun,2))*deleft
            p(3) = (aexp*eta(irun,3)+bexp*eta(jrun,3))*deleft
                                                              
            pax = p(1) - eta(irun,1)
            pay = p(2) - eta(irun,2)
            paz = p(3) - eta(irun,3)
                                                              
            pbx = p(1) - eta(jrun,1)
            pby = p(2) - eta(jrun,2)
            pbz = p(3) - eta(jrun,3)
                                                              
            prefa = dexp(-aexp*bexp*rAB/t1)*pi*anorm*bnorm/t1

            prefa = two*prefa                        
            expab = dexp(-aexp*bexp*rAB/t1)         
            s00=(pi/t1)**1.5*expab                  
            dum = one;  tf(1) = one; del = half/t1  
            do n = 2, maxall                        
               tf(n) = tf(n-1)*dum*del              
               dum = dum + two                      
            end do                                  
                                                    
c            ox0  = ovrlap(l1, l2,   pax, pbx, tf)   
c            oy0  = ovrlap(m1, m2,   pay, pby, tf)   
c            oz0  = ovrlap(n1, n2,   paz, pbz, tf)   
c            ox2  = ovrlap(l1, l2+2, pax, pbx, tf)   
c            oxm2 = ovrlap(l1, l2-2, pax, pbx, tf)   
c            oy2  = ovrlap(m1, m2+2, pay, pby, tf)   
c            oym2 = ovrlap(m1, m2-2, pay, pby, tf)   
c            oz2  = ovrlap(n1, n2+2, paz, pbz, tf)   
c            ozm2 = ovrlap(n1, n2-2, paz, pbz, tf)   
c            ov0  = ox0*oy0*oz0; ovl = ov0*s00       
c            ov1  = ox2*oy0*oz0; ov4 = oxm2*oy0*oz0  
c            ov2  = ox0*oy2*oz0; ov5 = ox0*oym2*oz0  
c            ov3  = ox0*oy0*oz2; ov6 = ox0*oy0*ozm2  
c
c            ovltot = ovltot + anorm*bnorm*ovl
c            
c            xl=dfloat(l2*(l2-1));   xm=dfloat(m2*(m2-1))               
c            xn=dfloat(n2*(n2-1));   xj=dfloat(2*(l2+m2+n2)+3)          
c            kin=s00*(bexp*(xj*ov0-two*bexp*(ov1+ov2+ov3))-
c     &          half*(xl*ov4+xm*ov5 +xn*ov6))   
c
c            kintot = kintot + anorm*bnorm*kin
c
            tnai = zero                                   
            m = imax + jmax + kmax -2                 
            do n = 1, imax
               sf(n,1) = fj(l1,l2,n-1,pax,pbx)      
            end do
            do n = 1, jmax
               sf(n,2) = fj(m1,m2,n-1,pay,pby)      
            end do
            do n = 1, kmax
               sf(n,3) = fj(n1,n2,n-1,paz,pbz)      
            end do

               pn = zero   
               cpx = p(1) - xp
               cpy = p(2) - yp
               cpz = p(3) - zp
               pcsq = cpx*cpx + cpy*cpy + cpz*cpz
               t = t1*pcsq

               call auxg(m,t,g)       

               epsi = quart/t1                         
               do ii = 1, 10
                  Airu(ii) = zero
                  Ajsv(ii) = zero
                  Aktw(ii) = zero
               end do
                                                          
               call aform(imax,sf,fact,cpx,epsi,Airu,1)
               call aform(jmax,sf,fact,cpy,epsi,Ajsv,2)
               call aform(kmax,sf,fact,cpz,epsi,Aktw,3)

               do ii = 1, imax                   
                  do jj = 1, jmax
                     do kk = 1, kmax
                        nu = ii + jj + kk - 2
                        pn = pn + Airu(ii)*Ajsv(jj)*Aktw(kk)*g(nu)
                     end do
                  end do
               end do
                                                      
               tnai = tnai - pn
c               write(*,*) tnai
                                                           
            totnai = totnai + prefa*tnai
         end do                 
      end do                

c      oeitot = totnai + kintot
      return
      end

      double precision function ovrlap(l1,l2,pax,pbx,tf)
      implicit double precision (a-h,o-z)
      integer l1, l2
      double precision pax,pbx
      double precision tf(*)
      
      double precision zero, one, dum
      data zero,one/0.0d00,1.0d00/
                                                               
      if( (l1.LT.0) .OR. (l2.LT.0) ) then 
         ovrlap = zero
         return
      end if
                                                             
      if ( (l1.EQ.0) .AND. (l2.EQ.0) ) then
         ovrlap = one
         return
      end if
                                                                   
      dum = zero;  maxkk = (l1 + l2)/2 + 1
                                                              
      do kk = 1, maxkk
         dum = dum + tf(kk)*fj(l1,l2,2*kk-2,pax,pbx)
      end do
                                                      
      ovrlap = dum
                                             
      return
      end

      double precision function fj(l, m, j, a, b)           
                                                            
      implicit double precision (a-h,o-z)
      integer l, m, j
      double precision a,b
                                                            
      double precision sum, term, aa, bb
      integer i, imax, imin
      double precision fact(20)
      data fact/1.0D00,1.0D00,2.0D00,6.0D00,24.0D00,120.0D00,
     &          720.0D00,5040.0D00,40320.0D00,362880.0D00,3628800.0D00,
     &          39916800.0D00,479001600.0D00,6227020800.0D00,6*0.0D00/
                                                           
      imax = min(j, l)
      imin = max(0, j-m)
                                                            
      sum = 0.0D00
      do i = imin, imax
                                                            
         term = fact(l+1)*fact(m+1)/(fact(i+1)*fact(j-i+1))
         term = term/(fact(l-i+1)*fact(m-j+i+1))
         aa = 1.0D00; bb = 1.0D00
         if ( (l-i) .NE. 0 ) aa = a**(l-i)
                                                            
         if ( (m+i-j) .NE. 0 ) bb = b**(m+i-j)
                                                            
         term= term*aa*bb
         sum = sum +term
                                                            
      end do
                                                            
      fj = sum
                                                            
      return
      end

      double precision function fmch(nu, x, y)
      implicit double precision (a-h,o-z)
      double precision x, y
      integer nu
      double precision ten, half, one, zero, rootpi4, xd
      double precision term, partialsum
      integer m, i, numberofterms, maxone, maxtwo
      data zero,half,one,rootpi4,ten
     & /0.0D00,0.5D00,1.0D00,0.88622692D00,10.0D00/
      data maxone/50/, maxtwo/200/

      m = nu
      a = dfloat(m)
      if (x.LE.ten) then
          a = a + half                             
          term = one/a
          partialsum = term
          do i = 2, maxone
             a = a + one
             term = term * x/a
             partialsum = partialsum + term
             if (term/partialsum < 1.0D-08) exit
          end do
                                                   
          if (i.EQ.maxone) then
             write(6,200)
200          format('i > 50 in fmch')
             STOP
          end if
          fmch = half*partialsum*y
          return 
      else
          b = a + half                                           
          a = a - half
          xd = one/x
          approx = rootpi4*dsqrt(xd)*xd**m 
          if (m.GT.0) then
              do i=1, m
                 b = b - one
                 approx = approx * b
              end do
          end if
          fimult = half*y*xd
          partialsum = zero
          
          if (fimult.EQ.zero) then
              fmch = approx
              return
          end if
                                                                 
          fiprop = fimult / approx
          term = one
          partialsum = term
          numberofterms = maxtwo
          do i=2, numberofterms
             term = term * a * xd
             partialsum = partialsum + term
             if (dabs(term*fiprop/partialsum).LE.1.0D-08) then
                 fmch = approx - fimult*partialsum
                 return
             end if
             a = a - one
          end do
          write(6,201)
 201      format(' numberofterms reached in fmch')
          STOP
      end if
      end

      subroutine auxg(mmax, x, g)                       
                                                        
      implicit double precision (a-h,o-z)
      integer mmax
      double precision x, g(*)
      double precision fmch
      double precision two, y
      integer mp1mx, mp1, md, mdm
      data two/2.0D00/
      y = dexp(-x)
      mp1mx = mmax+1
      g(mp1mx) = fmch(mmax, x, y)
      if ( mmax .LT. 1 ) go to 303
      do mp1 = 1, mmax                                
         md  = mp1mx - mp1
         mdm = md - 1
         g(md) = (two*x*g(md+1) + y )/dfloat(2*mdm+1)
      end do

 303  return
      end

      subroutine aform(imax, sf, fact, cpx, epsi, Airu, xyorz)
      implicit double precision (a-h,o-z)
      integer imax, xyorz
      double precision Airu(*), fact(*), sf(10,*)
      
      double precision one
      data one/1.0d00/

      do i = 1, imax
         ai = (-one)**(i-1)*sf(i,xyorz)*fact(i)
         irmax = (i-1)/2 + 1
         do ir = 1, irmax
            irumax = irmax -ir +1
            do iru = 1, irumax
               iq = ir + iru -2
               ip = i  -2*iq -1
               at5 = one
               if ( ip .GT. 0 ) at5 = cpx**ip
               tiru=ai*(-one)**(iru-1)*at5*epsi**iq
     &             /(fact(ir)*fact(iru)*fact(ip+1))
               nux = ip + iru
               Airu(nux) = Airu(nux) + tiru
            end do
         end do
      end do
                                
      return
      end

