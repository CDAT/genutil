C this routine was created by identifying a block of code from marsgo
C that could be isolated.  A fortran/python interface will be 
C constructed so it  functions properly.
      subroutine mrsgo1(
     1n, p, x, y, w, nk, lx, tb, cm, sc, db, d, m, mm, mmwork, mm1,
     1ms, fv, fln, tx, txl, tx1, txi, txt, df,
     1mk, mkp1, mkp2, nep, yb, sw, s, t, u, v, we, sy, ssq, rsq, se, 
     1xt, dy, kcp0, eps, big, nmin, alf, l, me, mel, mp, nop, 
     1jq, jp, ja, jas, lbf, nnt, kr, mn,
     1nc, kcp, k1, k, j, i, b, a)
     
      implicit none
      integer n, p, nk, ms, mm(n,p), mmwork(*), lx(p)
      logical elg, newbf
      double precision x(n,p), y(n), w(n), tb(5,nk), cm(*)
      double precision sc(n,*), tx(5), df, fv
      double precision db(n,*), d(nk,*)
      double precision yb, sw, s, t, u, v, we, sy, a, b
      double precision se, txt, xt, eps, rsq, ssq, dy
      integer nmin, k, jas, mp, ja, jp
      integer mm1, mn, me, mel, kcp, l, nnt, lbf, nc, nop
      integer k1, mk, nep, j, mkp1, mkp2, i, kcp0, m, kr, jq
      double precision big, fln, alf, txi, txl, tx1
      integer jf
C local parameters
      integer jd1, jd2, ict, newb, nnl, nst, nnr, j0, mj, h, sj, jft
      double precision dx, tt, st, su, yc, dv, xa, xb, xx, xd, fvr

C beginning of marsgo_block in python
      if(jq .ne. 0) go to 15
      jd1=1
      jd2=p
      go to 16
   15 jd1=jq
      jd2=jd1
      print *, 'in mrsgo1 jd1,jd2=', jd1,jd2
   16 do 52 jp=jd1,jd2
      print *, 'in mrsgo1 jp, m=', jp, m
      if(x(mm(1,jp),jp).ge.x(mm(n,jp),jp)) go to 52
      if(jf(l,jp,tb).ne.0) go to 52
      call isfac(l,jp,mm1,tb,cm,ja)
      if(ja.lt.0) go to 52
      if(.not.elg(jp,l,lx,tb,cm)) go to 52
      if(ja .ne. 0) go to 18
      if(lbf .eq. 0) go to 19
      call blf0(l,0,n,x,w,cm,sc,nnt,sc(1,mkp1))
      lbf=0
      call mnspan(ms,alf,nep,nnt,mn,me,mel)
      go to 19
   18 call blf0(l,ja,n,x,w,cm,sc,nnt,sc(1,mkp1))
      lbf=1
      if(nnt.le.nmin) go to 52
      call mnspan(ms,alf,nep,nnt,mn,me,mel)
      if(nnt.le.max0(me,mel)) go to 52
   19 fvr=1.0
      if(jft(mm1,jp,tb).eq.0) fvr=1.0+fv
      ict=0
      if(lx(jp) .ge. 0) go to 20
      ict=1
      nc=int(cm(2*jp+1)+.1)-int(cm(2*jp)+.1)+1
c formerly
c     call csp(jp,nc,m,n,x,y,w,nk,tb,cm,kcp,yb,d,kr,nnt,  sw,me,mkp2,nop
c    1,sc(1,mkp1),db,d(1,3), mm(1, p+1)) 
      call csp(jp,nc,m,n,x,y,w,nk,tb,cm,kcp,yb,d,kr,nnt,  sw,me,mkp2,nop
     1,sc(1,mkp1),db,d(1,3),mmwork) 
      print *, 'in mrsgo1 nop, m, ict=', nop, m, ict
      if(nop.eq.0) go to 52
      go to 45
c start of mrsgo2
   20 tb(2,m)=jp
      tb(3,m)=x(mm(1,jp),jp)
      tb(4,m)=l
      k1=kr
      ssq=rsq
      call update(1,n,m,kr,x,y,w,sw,yb,tb,cm,sc,sc(1,mkp1),db,d,d(1,3))
      if(kr .le. k1) go to 21
      rsq=rsq-d(kr,1)**2
      tb(1,m)=rsq/sw
      go to 22
   21 tb(1,m)=big
   22 continue
      print *, lx(jp), 3, m, mk, nnt , me+mel
      if((lx(jp) .ne. 3) .and. ((m .lt. mk) .and. (nnt .gt. me+mel))) go
     1 to 26
      tb(1,m)=rsq/sw
      newbf=newb(m,tb).eq.0
      if(fvr*tb(1,m) .gt. txl .or. .not.(newbf)) go to 23
      txl=fvr*tb(1,m)
      tx1=tb(1,m)
      jq=jp
      print *, fvr,tb(1,m),txi, .not.(newbf)
   23 if(fvr*tb(1,m) .gt. txi .or. .not.(newbf)) go to 25
      txi=fvr*tb(1,m)
      tx(1)=tb(1,m)
      do 24 i=2,4
      tx(i)=tb(i,m)
   24 continue
      jas=ja
   25 kr=k1
      rsq=ssq
      print *, 'rsq=', rsq
      go to 52
   26 mm1=m
      m=m+1
      tb(1,m)=big
      xa=0.0
      j=n
      nnl=nnt
      nst=0
      nnr=-1
   27 j0=j
   28 mj=mm(j,jp)
      h=sc(mj,mkp1)
      if(w(mj) .le. 0.0 .or. h .le. 0.0) go to 29
      nst=nst+1
      nnl=nnl-1
      nnr=nnr+1
   29 if(x(mm(j-1,jp),jp).lt.x(mm(j,jp),jp) .and.nst.ge.mn.and.nnl.ge.me
     1l.and.nnr.ge.me) go to 30
      j=j-1
      if(j.le.1) go to 30
      go to 28
   30 if(j.le.1) go to 45
      nst=0
      xb=xa
      xa=x(mm(j,jp),jp)
      if(j0 .ne. n) go to 34
      v=0.d0
      u=v
      t=u
      we=t
      se=we
      sy=se
      dy=sy
      i=1
      go to 32
   31 i=i+1
   32 if((i).gt.(kr)) go to 33
      d(i,2)=0.d0
      d(i,3)=d(i,2)
      go to 31
   33 txt=x(mm(1,jp),jp)+x(mm(n,jp),jp)
      xt=0.5*txt
      go to 37
   34 dx=xb-xa
      dy=dy+dx*sy
      we=we+dx*se
      v=v+dx*(2.d0*u-(xb+xa-txt)*t)
      i=1
      go to 36
   35 i=i+1
   36 if((i).gt.(kr)) go to 37
      d(i,2)=d(i,2)+dx*d(i,3)
      go to 35
   37 do 40 k=j,j0
      mj=mm(k,jp)
      h=sc(mj,mkp1)
      if(w(mj).le.0.0.or.h.le.0.0) go to 40
      xx=x(mj,jp)
      xd=xx-xa
      su=w(mj)*h
      st=su*xd
      yc=y(mj)-yb
      dy=dy+st*yc
      sy=sy+su*yc
      we=we+st
      se=se+su
      sj=w(mj)*h**2
      v=v+sj*xd**2
      t=t+sj
      u=u+sj*(xx-xt)
      i=1
      go to 39
   38 i=i+1
   39 if((i).gt.(kr)) go to 40
      tt=db(mj,i)
      d(i,2)=d(i,2)+st*tt
      d(i,3)=d(i,3)+su*tt
      go to 38
   40 continue
      dv=v-we**2/sw
      if(dv .le. 0.d0) go to 44
      a=0.d0
      b=a
      i=1
      go to 42
   41 i=i+1
   42 if((i).gt.(kr)) go to 43
      s=d(i,2)
      a=a+s*d(i,1)
      b=b+s**2
      go to 41
   43 b=dv-b
      if(b .le. eps*dv) go to 44
      b=-(dy-a)**2/b
      if(b .ge. tb(1,m)) go to 44
      tb(1,m)=b
      tb(3,m)=xa
   44 j=j-1
      if(j.le.1) go to 45
      go to 27
c start mrsgo3
   45 tb(2,m)=jp
      tb(4,m)=l
      tb(1,m)=(rsq+tb(1,m))/sw
      if(ict .ne. 0 .or. tb(1,mm1) .gt. fln*tb(1,m)) go to 46
      mp=mm1
      go to 47
   46 mp=m
   47 newbf=newb(mp,tb).eq.0
      if(fvr*tb(1,mp) .ge. txl .or. .not.(newbf)) go to 48
      txl=fvr*tb(1,mp)
      tx1=tb(1,mp)
      jq=jp
   48 if(fvr*tb(1,mp) .ge. txi .or. .not.(newbf)) go to 51
      txi=fvr*tb(1,mp)
      tx(1)=tb(1,mp)
      do 49 i=2,4
      tx(i)=tb(i,mp)
   49 continue
      jas=ja
      if(ict .eq. 0) go to 51
      do 50 i=1,nc
      cm(kcp0+i)=cm(kcp+i)
   50 continue
      kcp=kcp0+nc
      tx(3)=kcp0
      print *, 'ict,m=', ict,m
   51 if(ict .ne. 0) go to 52
      m=mm1
      mm1=m-1
      kr=k1
      rsq=ssq
      print *, 'm=', m
   52 continue
C end of marsgo_block in python

      return
      end
