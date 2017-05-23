c Changes to marsgo are made to continue dismantling the arrays used in the routine
c in particular the mm array is being split into mm & mmWORK. mm does not change
c in marsgo & any routine that it calls.  mmWork however appears to change in the 
c call to csp.  The strange thing is that while mm is dimensioned as nXp the call
c to csp is done using mm(1, p+1) which is out of range; fortran is a piece of shit
c to allow this.  The only thing that changes is the addition of mmWork as an argument
c to marsgo & csp.  There was a dummy array added to mars1 & cvmars of length 10;
c only to satisfy the compiler. I have no intention to use these subroutines.

c extracted from mars_nolog.f
      subroutine marsgo (n,p,x,y,w,nk,ms,df,fv,mi,lx,it,xm,xs,az,tb,cm,s
     1c,db,d,mm,mmwork)
      implicit none
      integer n,p,nk,ms,mi,it,mm(n,p),lx(p),mmwork(*)
      logical elg,newbf
      double precision x(n,p),y(n),w(n),xm(p),xs(p),tb(5,nk),cm(*),sc(n,
     1*),vcst(3),tx(5),df,fv,az
      double precision db(n,*),d(nk,*)
      double precision yb,yv,sw,s,t,u,v,we,sy,a,b,xb,xx,xd,ssq,alr
      double precision dx,wn,se,tt,txt,xt,st,su,yc,eps,rsq,dy,dv,asq0
      character*28 hol
      integer ix,nmin,k,jas,nnl,nst,nnr,j0,mj,mp,nli,kl,ll,ja,jp,
     1mm1,mn,me,mel,jd1,jd2,jft,nopt,kcp,l,nnord,nnt,lbf,ict,nc,nop,
     2k1,newb,nal,jn,mk,nep,j,mkp1,mkp2,i,kcp0,m,mtot,kr,jq
      double precision big,fln,alf,xa,h,sj,fjn,fkr,gcv,asm,tcsts,val,df1
     1,cst,cfac,txi,txl,fvr,tx1,tcst,tcmx,txm,xk
      double precision phi
      integer jf,ibfext
      data ix,alr,eps,big,fln,nmin,alf,vcst  /0,1.d-7,1.d-4,9.9e30,-1.0,
     15,.05,1.0,.666667,.333333/
c      if(it.gt.0) write(it,97)
      mk=nk
      df1=0.0
      nep=0
      t=0.d0
      u=0.d0
      v=0.d0
      we=0.d0
      sy=0.d0
      ssq=0.d0
      se=0.d0
      txt=0.d0
      xt=0.d0
      dy=0.d0
      kl=0
      kcp0=0
      tx1=0.0
      do 1 j=1,p
      if(lx(j).eq.0) go to 1
      if(x(mm(1,j),j).ge.x(mm(n,j),j)) go to 1
      nep=nep+1
      cst=vcst(iabs(lx(j)))
      if(mi.eq.1) cst=dmin1(cst,vcst(2))
      df1=df1+cst
    1 continue
      if(nep .ne. 0) go to 2
c     if(it.gt.0) write(it,'('' no predictor variables.'')')
      stop
    2 if(nep.eq.1) df1=vcst(3)
      cfac=df1/nep
      df1=df*cfac
      mkp1=mk+1
      mkp2=mk+2
      sw=0.d0
      wn=sw
      yb=wn
      s=yb
      do 3 i=1,n
      sw=sw+w(i)
      wn=wn+w(i)**2
      yb=yb+w(i)*y(i)
    3 continue
      yb=yb/sw
      wn=sw**2/wn
      do 4 i=1,n
      s=s+w(i)*(y(i)-yb)**2
    4 continue
      yv=s/sw
      tcst=1.0
      tcmx=wn-df1*vcst(1)-2.0
      if(cm(1) .le. 0.0) go to 7
      i=2
      go to 6
    5 i=i+(2)
    6 if((2)*((i)-(2*p)).gt.0) go to 7
      if(cm(i).gt.0.0) kcp0=cm(i+1)+.1
      go to 5
    7 m=0
      mtot=m
      txm=yv/(1.d0-1.d0/wn)**2
      rsq=yv*sw
      kr=0
      nopt=0
c      if(it.gt.0) write(it,98) m,txm,0.0,1.0
      if(fln.lt.0.0) fln=1.0+4.0/wn
      call addpar(0)
    8 if(m.ge.mk.or.tcst.ge.tcmx) go to 69
      nopt=nopt+1
      call itrpar(nopt)
      mm1=m
      m=m+1
      txi=big
      kcp=kcp0
      asq0=rsq/sw
    9 call nxtpar(l,jq)
      if(l.lt.0) go to 53
      txl=big
      if(nnord(l,tb) .lt. mi) go to 10
      call updpar(0,-1.d0)
      go to 9
   10 call blf0(l,0,n,x,w,cm,sc,nnt,sc(1,mkp1))
      lbf=0
      if(nnt .gt. nmin) go to 11
      call updpar(0,-1.d0)
      go to 9
   11 nep=0
      do 12 jp=1,p
      if(x(mm(1,jp),jp).ge.x(mm(n,jp),jp)) go to 12
      if(jf(l,jp,tb).ne.0) go to 12
      call isfac(l,jp,mm1,tb,cm,ja)
      if(ja.lt.0) go to 12
      if(.not.elg(jp,l,lx,tb,cm)) go to 12
      nep=nep+1
   12 continue
      if(nep .ne. 0) go to 13
      call updpar(0,-1.d0)
      go to 9
   13 call mnspan(ms,alf,nep,nnt,mn,me,mel)
      if(nnt .gt. max0(me,mel)) go to 14
      call updpar(0,-1.d0)
      go to 9
  14  continue  

C beginning of marsgo_block in python
      call mrsgo1(
     1n, p, x, y, w, nk, lx, tb, cm, sc, db, d, m, mm, mmwork, mm1,
     1ms, fv, fln, tx, txl, tx1, txi, txt, df,
     1mk, mkp1, mkp2, nep, yb, sw, s, t, u, v, we, sy, ssq, rsq, se, 
     1xt, dy, kcp0, eps, big, nmin, alf, l, me, mel, mp, nop,
     1jq, jp, ja, jas, lbf, nnt, kr, mn,
     1nc, kcp, k1, k, j, i, b, a)
C end of marsgo_block in python

      call updpar(jq,asq0-tx1)
      go to 9
   53 continue

      jp=tx(2)+.1
      call selpar(int(tx(4)+.1))
      if(cm(2*jp) .le. 0.) go to 54
      nc=int(cm(2*jp+1)+.1)-int(cm(2*jp)+.1)+1
      kcp0=kcp0+nc
   54 continue
      if(jas .le. 0) go to 60
      call getnst(jas,cm,jn,kcp,cm(kcp0+1))
      tb(2,m)=jn
      tb(3,m)=kcp0
      kcp0=kcp0+kcp
      tb(4,m)=tx(4)
      k1=kr
      call blf(int(tx(4)+.1),n,sc,sc(1,mkp1))
      tx(4)=m
      call update(2,n,m,kr,x,y,w,sw,yb,tb,cm,sc,sc(1,mkp1),db,d,d(1,3))
      if(kr.gt.k1) rsq=rsq-d(kr,1)**2
      call addpar(m)
      if(m .ge. mk) go to 58
      m=m+1
      tb(2,m)=-tb(2,m-1)
      do 55 i=3,4
      tb(i,m)=tb(i,m-1)
   55 continue
      if(ibfext(m,tb,cm) .eq. 0) go to 56
      m=m-1
      go to 58
   56 do 57 i=1,n
      sc(i,m)=phi(m,i,n,x,tb,cm)
   57 continue
      call addpar(m)
   58 if(it .le. 0) go to 59
      mp=m-1
      tcst=(nopt-1)*df1+kr+1.0
      fjn=jn
      fkr=kr
      gcv=(rsq/sw)/(1.d0-tcst/wn)**2
      call holl(jn,cm,tb(3,m),hol)
c      if(m.eq.mtot+1) write(it,100) m,gcv,fkr,tcst,fjn,hol,tb(4,m)
c      if(m.eq.mtot+2) write(it,99) m,mp,gcv,fkr,tcst,fjn,hol,tb(4,m)
   59 mtot=m
      m=m+1
      if(m.gt.mk) go to 69
   60 do 61 i=1,5
      tb(i,m)=tx(i)
   61 continue
      k1=kr
      call blf(int(tx(4)+.1),n,sc,sc(1,mkp1))
      call update(2,n,m,kr,x,y,w,sw,yb,tb,cm,sc,sc(1,mkp1),db,d,d(1,3))
      if(kr.gt.k1) rsq=rsq-d(kr,1)**2
      call addpar(m)
      if(m .ge. mk .or. (cm(2*jp) .le. 0.0) .and. (tx(3) .le. x(mm(1,jp)
     1,jp))) go to 66
      m=m+1
      do 62 i=1,4
      tb(i,m)=tx(i)
   62 continue
      tb(2,m)=-tb(2,m)
      if(cm(2*jp) .le. 0.0) go to 64
      do 63 i=1,n
      sc(i,m)=phi(m,i,n,x,tb,cm)
   63 continue
      go to 65
   64 k1=kr
      call update(2,n,m,kr,x,y,w,sw,yb,tb,cm,sc,sc(1,mkp1),db,d,d(1,3))
      if(kr.gt.k1) rsq=rsq-d(kr,1)**2
   65 call addpar(m)
   66 tcst=nopt*df1+kr+1.0
      if(it .le. 0) go to 68
      mp=m-1
      jp=dabs(tx(2))+.1
      fkr=kr
      gcv=(rsq/sw)/(1.d0-tcst/wn)**2
      if(cm(2*jp) .le. 0.0) go to 67
      call holl(jp,cm,tx(3),hol)
c      if(m.eq.mtot+1) write(it,100) m,gcv,fkr,tcst,tx(2),hol,tx(4)
c      if(m.eq.mtot+2) write(it,99) m,mp,gcv,fkr,tcst,tx(2),hol,tx(4)
      go to 68
   67 xk=xm(jp)+xs(jp)*tx(3)
c      if(m.eq.mtot+1) write(it,93) m,gcv,fkr,tcst,tx(2),xk,tx(4)
c      if(m.eq.mtot+2) write(it,94) m,mp,gcv,fkr,tcst,tx(2),xk,tx(4)
   68 mtot=m
      go to 8
   69 mk=min0(m,mk)
      m=mk+1
      k=m
      go to 71
   70 k=k+1
   71 if((k).gt.(nk)) go to 72
      tb(1,k)=0.0
      go to 70
   72 call sscp(n,m,sc,y,w,yb,yv,sw,db,d)
      call lsf1(db,m,d,yb,alr,b,d(1,2),a,d(1,3))
      nli=0
      do 73 k=1,mk
      if(d(k,2).ne.0.d0) nli=nli+1
   73 continue
      df1=df1*nopt+nli
      tcst=df1+1.0
      df1=df1/nli
      do 74 k=1,nk
      tb(5,k)=df1
   74 continue
      asm=(b/sw)/(1.d0-tcst/wn)**2
      tcsts=tcst
      az=a
      do 75 k=1,mk
      tb(1,k)=0.0
      if(d(k,2).ne.0.d0) tb(1,k)=d(k,2)
   75 continue
      if(ix .eq. 0) go to 81
      sc(1,1)=(cfac*nopt)/nli
      sc(2,1)=wn
      sc(3,1)=yv
      sc(4,1)=yb
      do 80 k=nli,nk
      call array(k+4,n,i,j)
      sc(i,j)=b/sw
      k1=k*(nk+1)+3
      l=0
      go to 77
   76 l=l+1
   77 if((l).gt.(nk)) go to 80
      k1=k1+1
      call array(k1,n,i,j)
      if(l .ne. 0) go to 78
      sc(i,j)=a
      go to 76
   78 if(l .le. mk) go to 79
      sc(i,j)=0.0
      go to 76
   79 sc(i,j)=d(l,2)
      go to 76
   80 continue
      call array((nk+1)**2+4,n,i,j)
      sc(i,j)=mk
      kl=nli
   81 do 88 ll=2,nli
      call bkstp(db,m,d,yb,alr,b,d(1,2),a,k,d(1,3))
      if(k.eq.0) go to 89
      if(ix .eq. 0) go to 86
      call array(kl+3,n,i,j)
      sc(i,j)=b/sw
      kl=kl-1
      k1=kl*(nk+1)+3
      l=0
      go to 83
   82 l=l+1
   83 if((l).gt.(nk)) go to 86
      k1=k1+1
      call array(k1,n,i,j)
      if(l .ne. 0) go to 84
      sc(i,j)=a
      go to 82
   84 if(l .le. mk) go to 85
      sc(i,j)=0.0
      go to 82
   85 sc(i,j)=d(l,2)
      go to 82
   86 tcst=tcst-df1
      b=(b/sw)/(1.d0-tcst/wn)**2
      if(b.ge.asm) go to 88
      asm=b
      tcsts=tcst
      az=a
      do 87 i=1,mk
      tb(1,i)=0.0
      if(d(i,2).ne.0.d0) tb(1,i)=d(i,2)
   87 continue
   88 continue
   89 if(txm .gt. asm) go to 91
      asm=txm
      tcsts=1.0
      az=yb
      do 90 i=1,nk
      tb(1,i)=0.0
   90 continue
   91 if(it .le. 0) go to 92
c      write(it,95)
      call coefpr(it,mk,az,tb,cm,xs)
c      write(it,96) asm,tcsts
   92 return
      entry setfln(val)
      fln=val
      return
      entry setalf(val)
      alf=val
      return
      entry setmin(nal)
      nmin=nal
      return
      entry setcta(val)
      vcst(2)=val
      return
      entry setctl(val)
      vcst(3)=val
      return
      entry xvmrgo(nal)
      ix=nal
      return
      entry setalr(val)
      alr=val
      return
   93 format('   ',i3,'    ',g12.4,2('   ',d5.1),'       ',d3.0,  '
     1 ',g12.4,'      ',d8.0)
   94 format(' ',i3,' ',i3,'  ',g12.4,2('   ',d5.1),'       ',d3.0,  '
     1    ',g12.4,'      ',d8.0)
   95 format(/,' final model after backward stepwise elimination:')
   96 format(/,'   (piecewise linear) gcv = ',g12.4,'   #efprms = ',d5.1
     1)
   97 format(//,' forward stepwise knot placement:',//  '  basfn(s)    g
     1cv      #indbsfns  #efprms',  '   variable      knot            pa
     1rent')
   98 format('   ',i3,'    ',g12.4,2('   ',d5.1))
   99 format(' ',i3,' ',i3,'  ',g12.4,2('   ',d5.1),'       ',d3.0,a28,d
     13.0)
  100 format('   ',i3,'    ',g12.4,2('   ',d5.1),'       ',d3.0,a28,d3.0
     1)
      end
