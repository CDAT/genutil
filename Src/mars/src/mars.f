c
c Multivariate Adaptive Regression Splines (MARS modeling, version 3.6).
c
c
c Coded and copyright (c) by Jerome H. Friedman (3/25/93).
c
c
c
c                         A Micro User's Guide
c                                  to
c                               MARS 3.6
c
c                          Jerome H. Friedman
c                 ex         Stanford University
c
c    MARS 3.6 is a collection of subroutines that implement the multivariate
c adaptive regression spline strategy for data fitting and function
c approximation described in Friedman (1991a, 1991b, 1993). It is a general-
c ization of MARS 1.0 described in Friedman (1988) and MARS 2.5 described in
c Friedman (1991a). It implements as a special case a palindromically
c invariant version of the TURBO fitting technique for smoothing and additive
c modeling described in Friedman and Silverman (1989).
c
c    These subroutines represent a set of tools that can be invoked from a user
c coded program to perform various analyses. The user routine is responsible for
c reading the data into memory and passing it to the MARS subroutines, along
c with the various parameter settings, as arguments. This set of subroutines

c can also form the basis for incorporating this methodology into a
c statistical language or package.
c
c    The user interface subroutines are:
c
c MARS: computes the mars model from the data and provides summary
c     information for interpreting it.
c
c MISS: sets up mars to handle missing values.
c
c PLOT: constructs graphical output, useful for interpreting the
c     continuous part of the mars model, in a format that can be plotted
c     with a local graphics package.
c
c CATPRT: prints tables useful for interpreting the purely categorical parts
c     of the mars model.
c
c SLICE: produces user selected lower dimensional representations of the
c    mars model to aid in interpretation.
c
c FMOD: computes response estimates for given predictor variable vectors
c    given the mars model.
c
c   It should be noted that this is a new methodology with which there is,
c as of this time, little collective experience.
c
c
c References:
c
c [1] Friedman, J. H. (1988). Fitting functions to noisy data in high
c     dimensions. Proc., Twentieth Symposium on the Interface, Wegman, Gantz,
c     and Miller, eds. American Statistical Association, Alexandria, VA. 3-43.
c
c [2] Friedman, J. H. (1991a). Multivariate adaptive regression splines
c     (with discussion).  Annals of Statistics, 19, 1-141 (March).
c
c [3] Friedman, J. H. (1991b). Estimating functions of mixed ordinal and
c     categorical variables using adaptive splines. Department of Statistics,
c     Stanford University, Tech. Report LCS108.
c
c [4] Friedman, J. H. (1993). Fast MARS. Department of Statistics,
c     Stanford University, Tech. Report LCS110.
c
c [5] Friedman, J. H. and Silverman, B. W. (1989). Flexible parsimonious
c     smoothing and additive modeling (with discussion). TECHNOMETRICS, 31,
c     3-39 (Feburary).
c
c
c
c User interface subroutines:
c
c       All arguments in the calling sequence of user called MARS
c       subroutines must be of the same type in the calling program
c       as indicated in the subroutine's documentation below. All
c       calling sequence arrays must be dimensioned in the calling
c       program as indicated in the documentation below. This includes
c       all workspace arrays. The leading dimensions of multidimensional
c       arrays must match, whereas the last dimension and singly 
c       dimensioned arrays can be as large or larger than that indicated.
c
c       More detailed explanations for some of the quantities below can
c       be found in the indicated sections of references [2], [3], and/or
c       [4] above.
c
c
c
c call mars (n,p,x,y,w,nk,mi,lx,fm,im,sp,dp,mm):
c
c input:
c n = number of observations.
c p = number of predictor variables per observation (integer).
c x(n,p) = predictor variable data matrix.
c y(n) = response value for each observation.
c w(n) = weight (mass) for each observation.
c nk = maximum number of basis functions.(Ref[2] Sec. 3.6, Ref[3] Sec. 2.3)
c mi = maximum number of variables per basis function (interaction level).
c      mi=1 => additive modeling (main effects only);
c      mi>1 => up to mi-variable interactions allowed.
c lx(p) = predictor variable flags; lx(i) corresponds to the ith variable:
c    lx(i) =  0 : exclude variable from model.
c             1 : ordinal variable - no restriction.
c             2 : ordinal variable that can only enter additively;
c                 no interactions with other variables.
c             3 : ordinal variable that can enter only linearly.
c            -1 : categorical variable - no restriction.
c            -2 : categorical variable that can only enter additively;
c                 no interactions with other variables.
c
c output:
c fm(3+nk*(5*mi+nmcv+6)+2*p+ntcv), im(21+nk*(3*mi+8)) = mars model.
c    (nmcv = maximum number of distinct values for any categorical variable;
c     ntcv = total number of distinct values over all categorical variables.)
c
c    note: upon return im(1) and im(2) contain the lengths of the fm and im
c          arrays (respectively) actually used by the program.
c
c workspace:
c sp(n*(max(nk+1,2)+3)+max(3*n+5*nk+p,2*p,4*n)+2*p+4*nk) : double precision.
c dp(max(n*nk,(nk+1)*(nk+1))+max((nk+2)*(nmcv+3),4*nk)) : double precision.
c mm(n*p+2*max(mi,nmcv)) : integer.
c
c
c defaults:
c    the following quanities are set to default values in mars. their values
c    can be changed by executing any of the following statements before the
c    call to mars, with the new value as the argument.
c
c call speed(is):
c is = speed acceleration factor (1-5).
c    larger values progressively sacrifice optimization thoroughness for
c    computational speed advantage. This usually results in marked decrease
c    in computing time with little or no effect on resulting approximation
c    accuracy (especially useful for exploratory work).
c    is = 1 => no acceleration.
c    is = 5 => maximum speed advantage.
c    (default: is=4) (Ref [4] Secs. 3.0 - 4.0)
c
c call logit(il):
c il = ordinary/logistic regression flag.
c    il=0 => ordinary least-squares regression.
c    il=1 => logistic regression.
c    (default: il=0). If logistic regression is selected (il=1) then
c    the response variable is assumed to take on only the values 0/1 and
c    the mars model is for the log-odds: f(x) = log (Pr(Y=1:x)/Pr(Y=0:x)).
c    (Ref[2] Sec. 4.5)
c
c call setdf(df):
c df = number of degrees-of-freedom charged for (unrestricted)
c      knot optimization. (default: df=3.0)
c      (Ref[2] Sec. 3.6, Ref[3] Sec. 2.3)
c
c call xvalid(ix):
c ix = control parameter for sample reuse technique used to automatically
c    estimate smoothing parameter df (see above) from the data.
c    ix = 0 => no effect (default). value used for df is set by user if
c            setdf(df) is called (see above), otherwise default value 
c            (df=3.0) is used.
c    ix > 0 => ix - fold cross-validation.
c    ix < 0 => single validation pass using every (-ix)th (randomly selected)
c            observation as an independent test set.
c    if ix.ne.0 then call setdf(df) (see above) has no effect. if ix > 0,
c    computation increases roughly by a factor of ix over that for ix = 0.
c    for ix < 0 computation increases approximately by a factor of two.
c    (Ref[3] Sec. 2.3)
c 
c call stseed(is):
c is = seed for internal random number generator used to group observation
c      subsets for validation (ix.ne.0). (default: is=987654321).
c
c call print(it):
c it = fortran file number for printed output. (it.le.0 => no printed output.)
c      note that this controls printed output for all user called mars
c      routines. (default: it=6).
c
c call setfv(fv):
c fv = (fractional) incremental penalty for increasing the number of variables
c      in the mars model. sometimes useful with highly collinear designs as it
c      may produce nearly equivalent models with fewer predictor variables,
c      aiding in interpretation. (fv .ge. 0)
c    fv=0.0  => no penalty (default).
c    fv=0.05 => moderate penalty.
c    fv=0.1  => heavy penality.
c    the best value depends on the specific situation and some user
c    experimentation using different values is usually required. this option
c    should be used with some care. (Ref[2] Sec. 5.3)
c
c call setic(ic):
c ic = flag restricting categorical - ordinal interactions.
c    ic=0 => no effect (default).
c    ic=1 => interactions between categorical and ordinal variables prohibited.
c    ic=2 => maximum number of ordinal variables participating in any
c            interaction is restricted to two. categorical interactions are
c            unrestricted.
c    the restrictions associated with a value of ic are imposed in addition
c    to those that are controlled by the mi and lx flags (see above).
c
c call setint(i,j,k):
c i,j = predictor variable numbers.
c   k = interaction flag:
c      k.eq.0 => interactions between variables i and j are prohibited.
c      k.ne.0 => interactions between variables i and j are permitted
c                if allowed by the mi, lx, and ic parameter values (see above).
c                (default)
c   i.eq.0 .or. j.eq.0 => reset to defaults for all predictor variables.
c
c   this call can be executed repeatedly to invoke (or remove) multiple
c   constraints.
c
c call nest (n,i,j,nv,vals);
c nests predictor variable i within categorical predictor variable j.
c links the existance of a value for var(i) to a subset of values for var(j).
c variable j must be categorical (lx(j) < 0, see above). (Ref[3] Sec. 3.3)
c
c n = same as in mars (see above).
c i = variable to be nested (ordinal or categorical).
c j = variable to which i is nested (categorical).
c nv = size of corresponding subset of values for variable j.
c vals(nv) = specific values of variable j for which values of variable i
c            are defined to exist.
c
c setting nv = 0 removes the nesting associated with i and j previously set.
c setting i = 0 and/or j = 0 removes all nesting (default).
c
c this call can be executed repeatedly to invoke (or remove) nestings.
c (recursive nesting is not implemented. j may not itself be nested to
c another variable.)
c
c note: variable nesting does NOT override the interaction constraints set by
c       calling setic or setint (see above). the mi and lx flags (see above)
c       do not limit interactions on variables to which others are nested.
c       they control nested and other nonnested variables as described above
c       except that interactions with variables to which others are nested
c       are ignored in applying the constraints.
c
c
c call setms(ms): 
c ms = minimum span (minimum number of observations between each knot).
c      ms .le. 0 => default value (depending on n and p) is used.
c      (default: ms=0). (Ref[2] Sec. 3.8)
c
c
c the following three routines (miss, mkmiss, xmiss) are used to enable
c mars to deal with various aspects of missing predictor values in the
c training data and/or future data to be predicted. for problems with
c no such missing values, these three routines are of no concern.
c
c
c call miss (n,p,x,lx,xm,flg,pn,xn,lxn,xs,xp);
c
c called  (optionally) before mars (see above) to indicate the presence of
c missing values in the predictor variable data matrix.
c
c sets up mars to accomodate missing values.
c produces as output transformations of some of the original mars input
c quanities defining the problem. the new transformed quanities define
c the same problem - but with missing values - and replace the corresponding
c original quanities in the call to mars (see above) and plot and slice
c (see below). in particular, a new predictor data matrix is created
c with extra (dummy) variables - one for each original variable with
c missing values - to indicate observations for which the corresponding
c original variable is not missing. each such original variable is
c automatically nested (see above) to this (corresponding) dummy variable.
c also produces as output a slicing vector to be used as input to
c slice (see below) to produce the mars model corresponding to no missing
c values on any of the original predictor variables, for interpretation.
c (Ref[3] Sec. 3.4)
c
c input:
c n,p,x,lx = original intended input to mars (see above).
c xm(p) = vector giving missing value flag for each original variable.
c         x(i,j) = xm(j) => x(i,j) is missing.
c flg = input to slice (see below).
c
c output:
c pn,xn,lxn = corresponding transformed quanities to be used as input to mars.
c pn = number of predictor variables in transformed data matrix (integer,
c      .le. 2*p)
c xn(n,pn) = transformed data matrix.
c lxn(pn) = predictor variable flags for transformed data matrix.
c xs(pn) = input for slice (see below). slicing vector that produces a slice
c          of the augmented predictor variable space to produce the nonmissing
c          mars model.
c xp(2*p+1) = input for xmiss (see below).
c
c notes:
c    (1) the value of the output quanity pn is less than or equal to 2*p.
c        the other output quanities (arrays) should be dimensioned
c        xn(n,2*p), lxn(2*p), xs(2*p) in the calling program to be safe.
c    (2) if there are no missing values then the output quanities 
c        pn,xn,lxn will be identical to p,x,lx respectively, and
c        xs(j)=flg, j=1,p.
c    (3) the corresponding quanities for input and output can be the same in
c        the calling program. however, they should have the corresponding
c        dimension equal to 2*p, and they will be altered.
c    (4) dimensions of all relevant workspace arrays in mars and other
c        user callable routines must be large enough to accomodate the
c        the increased number of variables pn (.le. 2*p) produced.
c
c
c call mkmiss (n,p,x,y,w,xm,pm,nnx,nn,xnn,yn,wn,sc);
c
c called (optionally) before miss (see above) to generate additional data
c with missing values by resampling from original (input) data.
c
c used to train mars for missing predictor values when they are under
c represented in the training data. takes as input original data and
c produces as output a new (larger) data set containing both the original
c data and additional data resampled from it with predictor variable values
c flagged as missing. (Ref[3] Sec. 3.4)
c
c input:
c n,p,x,y,w = original data set (see mars above).
c xm(p) = same as in miss (see above).
c pm(p) = vector of fractions of missing values in output sample:
c    pm(j) = fraction of the total output data with jth predictor variable
c            value missing.
c nnx = maximum sample size of new (output) data.
c
c output:
c nn = sample size of generated output data set (input to miss, mars, and
c      other user called routines, in place of n).
c xnn(nn,p) = new generated predictor data matrix to be used as input to
c             miss (see above) in place of x.
c yn(nn),wn(nn) = new generated data to be used as input to mars and other
c                 user called routines in place of y and w.
c
c workspace:
c sc(p*nnx): double precision.
c
c notes:
c    (1) the output arrays should be dimensioned xnn(nnx,p),yn(nnx),wn(nnx)
c        in the calling program to be safe.
c    (2) if much larger fractions of missing values are requested than
c        exist in the orginal training data then the size of the output
c        (resampled) data set can become very large (nn = big).
c    (3) if the size of the output (resampled) data set reaches the
c        specified maximum (nnx) then the actual fractions of missing
c        values for each variable will be less than those specified in the
c        input vector pm(p).
c    (4) dimensions of all relevant workspace arrays in mars and other
c        user callable routines must be large enough to accomodate the
c        the increased number of observations nn (.le. nnx) produced.
c
c
c call xmiss (n,x,xm,xp,xn);
c
c must be called before fmod (see below) if miss (see above) was called
c before mars (see above) to handle missing values. produces a new set
c of covariate vectors, from the original set intended for fmod, to replace
c the original set in the call to fmod. (Ref[3] Sec. 3.4)
c
c input:
c n = number of covariate vectors (see fmod below).
c x(n,p) = original covariate vectors intended for fmod.
c xm, xp = same as in miss (see above).
c
c output:
c xn(n,pn) = new set of covariate vectors to be used as input to fmod.
c
c notes:
c    (1) the value of the quanity pn is less than or equal to 2*p (p = the
c        number of original predictor variables). the output quanity (array)
c        should be dimensioned xn(n,2*p) in the calling program to be safe.
c    (2) the corresponding quanities x and xn can be the same in
c        the calling program. however, it should have the second
c        dimension equal to 2*p, and it will be altered.
c
c
c
c
c The following subroutines can be called only after mars.
c
c
c
c call plot (m,n,p,x,fm,im,ngc,ngs,icx,nk,nc,crv,ns,srf,sp,mm):
c
c computes plots of all purely additive and all purely bivariate ordinal
c contributions to the mars model,in a form suitable for displaying
c with a computer graphics package. If there are no interactions between
c categorical and ordinal variables present in the mars model, this
c subroutine returns plots of the (compound) anova functions (Ref[2] Sec.3.5,
c eqns (25) and (27)) for those ordinal variables involved in at most two
c (ordinal) variable interactions. If categorical-ordinal interactions are
c present, it returns the corresponding plots for each ordinal contribution
c in the categorical-ordinal decomposition (Ref[3] Sec. 3.2).  since the
c locations of the plotted functions are arbitrary, they are all translated
c to have zero minimum value.
c
c input:
c   m = model flag:
c     = 1 => plot piecewise-linear mars model.
c     = 2 => plot piecewise-cubic mars model. (Ref[2] Sec. 3.7)
c   n = number of observations.
c   p = number of predictor variables per observation (integer).
c   x,fm,im = same as in mars (see above).
c   ngc = number of raster points for computing curve estimates.
c   ngs = number of raster points on each axis for computing surface estimates.
c   icx = convex hull flag:
c       = 0 => plot surface estimates over entire range of argument limits.
c       > 0 => plot surface estimates only inside the convex hull of the
c             (bivariate) point set.
c   nk = maximum number of basis functions.(Ref[2] Sec. 3.6, Ref[3] Sec. 2.3)
c
c output:
c    nc = number of curves (purely additive ordinal anova functions).
c    crv(ngc,2,nc) = additive ordinal anova functions in anova
c                    decomposition order.
c       crv(.,1,m) = ordered abscissa values for mth anova function.
c       crv(.,2,m) = corresponding ordinate values.
c    ns = number of surfaces (purely bivariate ordinal anova functions).
c    srf(ngs,ngs,ns) = bivariate (plus associated univariate) ordinal
c                      anova functions in anova decomposition order.
c       srf(.,.,m) = total contribution (bivariate + univariate) of ordinal
c                  variables associated with mth bivariate anova function of the
c                  mars model, in a square raster format (ngs x ngs) over the
c                  ranges of the two variables. the first and second indicies
c                  correspond to the first and second variables respectively.
c
c workspace:
c    sp(max(4*ngs*ngs,ngc,2*n)) : double precision.
c    mm(max(2*(mi+1),nmcv)) : integer.
c
c
c
c call catprt (m,nk,fm,im,sp,mm):
c
c prints all univariate and bivariate contributions to the purely categorical
c part of the mars model that are not involved in higher order (categorical)

c interactions. These are the (compound) anova functions (Ref[2] Sec. 3.5,
c eqns (25) and (27)) for the purely categorical part of the categorical-
c ordinal decomposition (Ref[3] Sec. 3.2, eqns (32b) and (38b)). since
c the locations of the printed functions are arbitrary they are all
c translated to have zero minimum value. the functions are printed as tables
c of integers in the range [0,99] or [0,9]. the function value for each
c entry is the product of the corresponding integer and the scale factor
c printed above the table.
c
c input:
c m = model flag:
c   = 1 => piecewise-linear model.
c   = 2 => piecewise-cubic  model. (Ref[2] Sec. 3.7)
c nk = maximum number of basis functions.(Ref[2] Sec. 3.6, Ref[3] Sec. 2.3)
c fm,im = same as in mars (see above).
c
c workspace:
c sp(nmcv*nmcv) : double precision.
c mm(2*nk+nmcv) : integer.
c
c
c
c call slice (flg,xs,n,p,x,nk,mi,fm,im,fmn,imn,sp,mm):
c
c computes the mars model within a lower dimensional subspace defined by an
c axis oriented slice of the predictor variable space. the slice is selected
c by assigning specific values to a subset of the predictor variables. the
c returned model is a function of the variables complement to the selected
c subset, and represents the mars model conditioned on the specified values
c for the selected variables. this new lower dimensional sliced model can
c be input to plot and/or catprt for interpretation in the same manner as the
c original mars model. (Ref[2] Sec. 4.7, Ref[3] Sec. 2.4)
c
c input:
c flg = flag for indicating that a predictor variable is not in the subset
c    defining the slice. its value should be outside the range of values for
c    all predictor variables.
c xs(p) = vector defining the slice:
c    xs(i).eq.flg => do not condition on ith variable.
c    xs(i).ne.flg => condition on ith variable at value stored in xs(i).
c n = number of observations.
c p = number of predictor variables per observation (integer).
c x = same as in mars (see above).
c nk = maximum number of basis functions.(Ref[2] Sec. 3.6, Ref[3] Sec. 2.3)
c mi = maximum number of variables per basis function (interaction level).
c fm,im = arrays defining mars model (output from mars - see above).
c
c output:
c fmn,imn = corresponding arrays defining the sliced model
c             (dimensioned same as in mars - see above).
c
c workspace:
c sp(2*nk+2*p+max(nk,3*p)) : double precision.
c mm(2*p) : integer.
c
c
c
c call fmod (m,n,x,fm,im,f,sp):
c
c calculates mars model response estimates for sets of covariate vectors.
c
c input:
c m = model flag:
c   = 1 => piecewise-linear mars model.
c   = 2 => piecewise-cubic mars model. (Ref[2] Sec. 3.7)
c n = number of covariate vectors.
c x(n,p) = covariate vectors.
c fm,im = same as in mars (see above).
c
c output:
c f(n) = value of the mars model estimate for each covariate vector.
c
c workspace:
c sp(n,2) : double precision.
c
c
c
c call cvinfo (dfs,pse,nbf):
c 
c returns results of sample reuse procedure for estimating optimal smoothing
c parameter df (see above). can only be called if xvalid(ix) was called
c before mars with ix.ne.0 (see above). (Ref[2] Sec 3.6, Ref[3] Sec. 2.3)
c
c output:
c dfs = optimal smoothing parameter estimate.
c pse = estimate of corresponding predictive-squared-error.
c nbf = estimate of associated number of (nonconstant) basis functions.
c
c 
c
      subroutine mars (n,p,x,y,w,nk,mi,lx,fm,im,sp,dp,mm)
      implicit none
      integer n,p,nk,mi,lx(*),im(*),mm(*)
      integer lcm
      double precision x(*),y(*),w(*),fm(*),sp(*)
      real tstart,ttaken
      double precision dp(*)
      character*120 logfile, fn
c      data logfile /'mars.log            '/
      data logfile /'mars.log'/
      tstart = (secnds(0.0))
      open(99, file=logfile)
c      print *,logfile
      im(3)=n
      im(4)=p
      im(5)=nk
      im(6)=mi
      im(7)=16
c     the following was changed from the original: 5*nk replace by 5*(nk+1)
c     this fixes the bug in purcat which surfaced when debugging the python
c     
      im(8)=im(7)+5*(nk+1)
      im(9)=im(8)+2*nk*mi
      im(10)=im(9)+3*(nk+2)
      im(2)=im(10)+nk*mi-1
      im(11)=1
      im(12)=2
      im(13)=im(12)+5*nk
      im(14)=im(13)+1
      im(15)=im(14)+nk*(5*mi+1)
      call mars1(n,p,x,y,w,nk,mi,lx,fm(im(11)),fm(im(12)),fm(im(15)),im(
     1im(7)),  im(im(8)),im(im(9)),im(im(10)),fm(im(13)),fm(im(14)),sp,d
     1p,mm)
      im(1)=im(15)+lcm(p,nk,fm(im(12)),fm(im(15)))-1
      ttaken = secnds(tstart)
c      write(99,100) ttaken
      print 100, ttaken
      close(99)
      return
  100 format(/,' MARS run time = ',f13.6, ' seconds')
      entry setlog(fn)
c      print *, fn
      logfile = fn
      return
      end
      subroutine plot (m,n,p,x,fm,im,ngc,ngs,icx,nk,nc,crv,ns,srf,sp,mm)
      implicit none
      integer m,n,p,ngc,ngs,icx,nk,nc,ns,im(*),mm(*)
      double precision x(n,p),fm(*),crv(*),srf(*),sp(*)
      if(m .ne. 1) go to 1
      call plotl(im(3),im(4),x,im(5),im(im(7)),im(im(8)),im(im(9)),im(im
     1(10)),  fm(im(12)),fm(im(15)),ngc,ngs,icx,nc,crv,ns,srf,sp,mm)
      return
    1 call plotc(im(3),im(4),x,im(5),im(im(7)),im(im(8)),im(im(9)),im(im
     1(10)),  fm(im(14)),fm(im(15)),ngc,ngs,icx,nc,crv,ns,srf,sp,mm)
      return
      end
      subroutine catprt (m,nk,fm,im,sp,mm)
      implicit none
      integer m,im(*),mm(*),nk
      double precision fm(*),sp(*)
      call ctprt1(m,im(5),im(im(7)),im(im(8)),fm(im(12)),fm(im(15)),fm(i
     1m(14)),sp,mm)
      return
      end
      subroutine slice (flg,xs,n,p,x,nk,mi,fm,im,fmn,imn,sp,mm)
      implicit none
      integer n,p,nk,mi,im(*),imn(*),mm(*)
      double precision flg,xs(*),x(*),fm(*),fmn(*),sp(*)
      integer i
      do 1 i=1,15
      imn(i)=im(i)
    1 continue
      i=im(15)
      go to 3
    2 i=i+1
    3 if((i).gt.(im(1))) go to 4
      fmn(i)=fm(i)
      go to 2
    4 call slice1(flg,xs,im(3),im(4),x,im(5),fm(im(11)),fm(im(12)),fm(im
     1(15)),  im(im(7)),im(im(8)),im(im(9)),im(im(10)),fm(im(13)),fm(im(
     114)),  fmn(im(11)),fmn(im(12)),imn(im(7)),imn(im(8)),imn(im(9)),im
     1n(im(10)),  fmn(im(13)),fmn(im(14)),sp,mm)
      return
      end
      subroutine fmod (m,n,p,x,fm,im,f,sp)      
      implicit none
      integer m,n,im(*), p
      double precision x(*),fm(*),f(*),sp(*)
      if(m .ne. 1) go to 1
      call fmrs(n,x,im(5),fm(im(11)),fm(im(12)),fm(im(15)),f)
      return
    1 call cmrs(n,x,fm(im(15)),im(im(7)),im(im(8)),im(im(9)),im(im(10)),
     1  fm(im(13)),fm(im(14)),f,sp)
      return
      end
      subroutine print(it)
      implicit none
      integer it
      call printm(it)
      call printg(it)
      call printc(it)
      call prtslc(it)
      return
      end
      subroutine setint(i,j,k)
      implicit none
      integer i,j,k
      integer il,m1,m2,l,ig,ll,it
      integer mlist
      parameter(mlist=1000)
      integer m(2,mlist)
      save m
      data il /0/
      if((i .ne. 0) .and. (j .ne. 0)) go to 1
      il=0
      return
    1 if(i.eq.j) return
      m1=min0(i,j)
      m2=max0(i,j)
      if(k .ne. 0) go to 6
      l=1
      go to 3
    2 l=l+1
    3 if((l).gt.(il)) go to 4
      if(m1.eq.m(1,l).and.m2.eq.m(2,l)) return
      go to 2
    4 il=il+1
      if(il .le. mlist) go to 5
      write(6,  '('' increase parameter mlist in subroutine setint to gr
     1eater than'',           i5,/,'' and recompile.'')') il
      stop
    5 m(1,il)=m1
      m(2,il)=m2
      return
    6 ig=0
      l=1
      go to 8
    7 l=l+1
    8 if((l).gt.(il)) go to 10
      if(m1 .ne. m(1,l) .or. m2 .ne. m(2,l)) go to 7
      ig=1
   10 if(ig.eq.0) return
      il=il-1
      ll=l
      go to 12
   11 ll=ll+1
   12 if((ll).gt.(il)) go to 13
      m(1,ll)=m(1,ll+1)
      m(2,ll)=m(2,ll+1)
      go to 11
   13 return
      entry intlst(it)
      if(it.le.0) return
      if(il.eq.0) return
      write(it,'(/,'' interactions prohibited between:'')')
c     do 14 l=1,il
c     write(it,'(''    var('',i3,'')  and  var('',i3,'')'')') m(1,l),m(2
c    1,l)
c  14 continue
      return
      entry intalw(i,j,k)
      k=1
      m1=min0(i,j)
      m2=max0(i,j)
      l=1
      go to 16
   15 l=l+1
   16 if((l).gt.(il)) go to 18
      if(m1 .ne. m(1,l) .or. m2 .ne. m(2,l)) go to 15
      k=0
   18 return
      end
      subroutine mars1 (n,p,x,y,w,nk,mi,lx,az,tb,cm,kp,kv,lp,lv,bz,tc,sp
     1,dp,mm)
      implicit none
      integer n,p,nk,mi,kp(5,*),kv(2,*),lp(3,*),lv(*),mm(n,*),lx(p)
      double precision az,bz,x(n,p),y(n),w(n),tb(5,nk),cm(*),tc(*),sp(*)
      double precision dp(*)
      integer ms,il,it,ic,ix,i1,im,is,i2,mal,i,j,k, ii, jj
      double precision df,fv,sw,wn,ef,s,yh,gcv,val,t
      logical stelg,z00001
      data ms,df,il,fv,it,ic,ix /0,3.0,0,0.0,99,0,0/

      if(it.gt.0) write(it,11)
      if(it.gt.0) write(it,10) n,p,nk,ms,mi,df,il,fv,ic
      if(it.gt.0) write(it,12)
      if(it.gt.0) write(it,'('' var: '',5('' '',20i3,/))') (i,i=1,p)
      if(it.gt.0) write(it,'('' flag:'',5('' '',20i3,/))') (lx(i),i=1,p)
c     print *, ' '
c     do 321 i = 1, n
c        print *,'M1 ',x(i,1),' ',x(i,2),' ',x(i,3),' ',x(i,4),
c    1         ' ',y(i)
c321  continue
      call intlst(it)
      call nstlst(it)
      i1=max0(n*(nk+1),2*n)+1
      im=i1+n+max0(3*n+5*nk,2*p,4*n,2*n+5*nk+p)
      is=im+p
      i2=max0(n*nk,(nk+1)*(nk+1))+1
      call rspnpr(it,il,n,y,w,mm)
      do 2 j=1,p
      do 1 i=1,n
      mm(i,j)=i
    1 continue
      call psort(x(1,j),mm(1,j),1,n)
    2 continue
      call ordpr(it,n,p,x,lx,mm)
      call atoscl (n,p,w,x,lx,mm,sp(im),sp(is),cm,x)
      call catpr(it,n,p,x,cm,mm(1,p+1))
      call oknest(it,p,lx,cm)
      if(ix.ne.0) call cvmars (ix,n,p,x,y,w,nk,ms,df,fv,mi,lx,it,sp(im),
     1sp(is),tb,cm,sp,dp,dp(i2),mm, sp(is+p),sp(is+p+2*n))

      call marsgo  (n,p,x,y,w,nk,ms,df,fv,mi,lx,it,sp(im),sp(is),az,tb,c
     1m,sp,dp,dp(i2),mm)
      if(il .le. 0) go to 6
      call logitl(n,x,y,w,nk,il,az,tb,cm,sp,dp)
      if(it .le. 0) go to 6
      sw=0.0
      wn=sw
      do 3 i=1,n
      sw=sw+w(i)
      wn=wn+w(i)**2
    3 continue
      wn=sw**2/wn
      ef=1.0
      do 4 k=1,nk
      if(tb(1,k).ne.0.0) ef=ef+tb(5,k)
    4 continue
      ef=1.0/(1.0-ef/wn)**2
      s=0.0
      t=s
      call fmrs(n,x,nk,az,tb,cm,sp)
      do 5 i=1,n
      yh=1.0/(1.0+exp(-sp(i)))
      gcv=ef*(y(i)-yh)**2
      s=s+w(i)*gcv
      t=t+w(i)*yh*(1.0-yh)
    5 continue
      s=s/sw
      t=t/sw
      write(it,13) s,t
    6 if(it .le. 0) go to 7

      if(il.eq.0) call anova (n,x,y,w,nk,it,tb,cm,lp,lv,sp,dp)
      if(il.gt.0) call anoval(n,x,y,w,nk,il,it,az,tb,cm,lp,lv,sp,dp)
    7 call ccoll (nk,tb,cm,kp,kv,lp,lv,mm)
      call cubic (n,p,x,y,w,nk,it,tb,cm,kp,kv,lp,lv,bz,tc,sp,sp(i1),sp(i
     11+2*p),mm,dp)
      if(il .le. 0) go to 9
      call logitc(n,x,y,w,nk,il,cm,tb,kp,kv,lp,lv,bz,tc,sp,sp(i1+4*n),
     1dp)
      if(it .le. 0) go to 9
      call cmrs(n,x,cm,kp,kv,lp,lv,bz,tc,sp,sp(n+1))
      s=0.0
      t=s
      do 8 i=1,n
      yh=1.0/(1.0+exp(-sp(i)))
      gcv=ef*(y(i)-yh)**2
      s=s+w(i)*gcv
      t=t+w(i)*yh*(1.0-yh)
    8 continue
      s=s/sw
      t=t/sw
      write(it,14) s,t
    9 if(it.gt.0) call varimp (n,p,x,y,w,nk,il,it,az,tb,cm,sp,sp(p+1),dp
     1)
      call orgpl(sp(im),sp(is),nk,tb,cm)
      call orgpc(sp(im),sp(is),lp,lv,tc)
      call sclato(n,p,x,sp(im),sp(is),cm,x)
      return
      entry setms(mal)
      ms=mal
      return
      entry setdf(val)
      df=val
      return
      entry printm(mal)
      it=mal
      return
      entry logit(mal)
      il=mal
      return
      entry setfv(val)
      fv=val
      return
      entry setic(mal)
      ic=mal
      z00001=stelg(ic)
      return
      entry xvalid(mal)
      ix=mal
      call xvmrgo(ix)
      return
   10 format(/' input parameters (see doc.):',/,  '    n     p    nk
     1ms    mi     df    il    fv     ic',/,  ' ',i5,i5,i6,i6,i6,f8.3,i5
     1,f7.3,i6)
   11 format(/,' MARS modeling, version 3.6 (3/25/93)',/)
   12 format(/' predictor variable flags:')
   13 format(/' piecewise-linear logistic gcv =',g12.4,'   ave var =',g1
     12.4)
   14 format(/' piecewise-cubic logistic gcv =',g12.4,'   ave var =',g12
     1.4)
      end
      subroutine plotc (n,p,x,nk,kp,kv,lp,lv,tc,cm,ngc,ngs,icx,nc,crv,ns
     1,srf,sp,mm)
      implicit none
      integer n,p,nk,ngc,ngs,icx,nc,ns,kp(5,*),kv(2,*),lp(3,*),lv(*),
     1mm(*)
      double precision x(n,p),tb(5,nk),tc(*),cm(*),crv(ngc,2,*),srf(ngs,
     1ngs,*),sp(*),zl(2),zu(2)
      integer it,jnt,iz,ll,nf,i,j,k,k1,k2,k4,ko,l,ncx,jj,jv,m,l1,l2,ne,
     1nh,nal,ngsq,nxs
      double precision big,d,dc,r,dl,fx,d1,d2
      data big,it /1.e30,6/
      if(it.gt.0) write(it,'(/'' mars graphics (piecewise-cubic):'',/)')
      jnt=2
      go to 1
      entry plotl (n,p,x,nk,kp,kv,lp,lv,tb,cm,ngc,ngs,icx,nc,crv,ns,srf,
     1sp,mm)
      if(it.gt.0) write(it,'(/'' mars graphics (piecewise-linear):'',/)'
     1)
      jnt=1
    1 ngsq=ngs**2
      iz=2*ngsq
      d=1.0/(ngs-1)
      dc=1.0/(ngc-1)
      ll=1
      nc=0
      ns=nc
    2 if(kp(1,ll).lt.0) go to 36
      if(kp(3,ll) .gt. 0) go to 3
      ll=ll+1
      go to 2
    3 nf=kp(3,ll)
      k4=kp(4,ll)-1
      k1=kp(1,ll)
      k2=kp(2,ll)
      if(it .le. 0) go to 7
      if(k1 .ne. 0) go to 4
      write(it,'('' pure ordinal contribution:'')')
      go to 7
    4 continue
      write(it,'('' categorical - ordinal interaction:'')')
      do 6 i=1,k1
      jj=kv(1,k2+i-1)
      j=iabs(jj)
      k=kv(2,k2+i-1)
      ncx=int(cm(2*j+1)+.1)-int(cm(2*j)+.1)+1
      do 5 l=1,ncx
      mm(l)=cm(k+l)+.1
      if(jj.lt.0) mm(l)=mod(mm(l)+1,2)
    5 continue
      write(it,'('' x('',i3,'') ='',70i1/80i1)') j,(mm(l),l=1,ncx)
    6 continue
    7 do 35 k=1,nf
      l=lp(1,k+k4)
      if(l.gt.2) go to 35
      ko=lp(2,k+k4)
      if(l .ne. 1) go to 17
      j=0
      jv=lv(ko)
      do 9 m=k,nf
      l1=lp(1,m+k4)
      if(l1.eq.1) go to 9
      l2=lp(2,m+k4)-1
      do 8 i=1,l1
      if(jv.eq.lv(l2+i)) j=1
    8 continue
      if(j.eq.1) go to 10
    9 continue
   10 if(j.eq.1) go to 35
      nc=nc+1
      zl(1)=big
      zu(1)=-big
      do 11 i=1,n
      r=x(i,jv)
      zl(1)=dmin1(zl(1),r)
      zu(1)=dmax1(zu(1),r)
   11 continue
      dl=(zu(1)-zl(1))*dc
      do 12 i=1,ngc
      crv(i,1,nc)=zl(1)+dl*(i-1)
   12 continue
      if(jnt .ne. 1) go to 13
      call fun(l,jv,ngc,crv(1,1,nc),nk,tb,cm,k1,kv(1,k2),crv(1,2,nc),mm)
      go to 14
   13 call cfun (l,jv,ngc,crv(1,1,nc),nf,lp(1,k4+1),lv,tc(kp(5,ll)),  cr
     1v(1,2,nc),sp,mm)
   14 dl=big
      do 15 i=1,ngc
      dl=dmin1(dl,crv(i,2,nc))
   15 continue
      fx=0.0
      do 16 i=1,ngc
      crv(i,2,nc)=crv(i,2,nc)-dl
      fx=dmax1(fx,crv(i,2,nc))
   16 continue
      if(it.gt.0) write(it,39) nc,jv,fx
      go to 35
   17 j=0
      mm(1)=lv(ko)
      mm(2)=lv(ko+1)
      do 19 m=k,nf
      l1=lp(1,m+k4)
      if(l1.le.2) go to 19
      l2=lp(2,m+k4)-1
      do 18 i=1,l1
      if(mm(1).eq.lv(l2+i).or.mm(2).eq.lv(l2+i)) j=1
   18 continue
      if(j.eq.1) go to 20
   19 continue
   20 if(j.eq.1) go to 35
      ns=ns+1
      zl(1)=big
      zl(2)=zl(1)
      zu(1)=-big
      zu(2)=zu(1)
      do 22 j=1,2
      do 21 i=1,n
      r=x(i,mm(j))
      zl(j)=dmin1(zl(j),r)
      zu(j)=dmax1(zu(j),r)
   21 continue
   22 continue
      do 23 j=1,2
      dl=(zu(j)-zl(j))/(ngs-3)
      zu(j)=zu(j)+dl
      zl(j)=zl(j)-dl
   23 continue
      ne=0
      d1=d*(zu(1)-zl(1))
      d2=d*(zu(2)-zl(2))
      do 25 j=1,ngs
      do 24 i=1,ngs
      ne=ne+1
      sp(iz+ne)=zl(1)+d1*(i-1)
      sp(iz+ngsq+ne)=zl(2)+d2*(j-1)
   24 continue
   25 continue
      dl=big
      if(jnt .ne. 1) go to 26
      call pair(mm,ngsq,sp(iz+1),nk,tb,cm,k1,kv(1,k2),  srf(1,1,ns),sp,m
     1m(3))
      go to 27
   26 call cpair(mm,ngsq,sp(iz+1),nf,lp(1,k4+1),lv,  tc(kp(5,ll)),srf(1,
     11,ns),sp)
   27 if(icx .le. 0) go to 29
      call cvxhul(n,x(1,mm(1)),x(1,mm(2)),big,nh,sp)
      if(it .le. 0 .or. 3*nh .lt. iz) go to 28
      nxs=dsqrt(dfloat(3*nh)*0.5D0)+1.1D0
c     write(it,38) nxs
   28 call hulset(ngsq,sp(iz+1),big,nh,sp,srf(1,1,ns))
   29 do 31 j=1,ngs
      do 30 i=1,ngs
      if(i.eq.1.or.j.eq.1.or.i.eq.ngs.or.j.eq.ngs.or.srf(i,j,ns).ge.big)
     1 go to 30
      dl=dmin1(dl,srf(i,j,ns))
   30 continue
   31 continue
      fx=0.0
      do 34 j=1,ngs
      do 33 i=1,ngs
      if((i .ne. 1) .and. ((j .ne. 1) .and. ((i .ne. ngs) .and. ((j .ne.
     1 ngs) .and. (srf(i,j,ns) .lt. big))))) go to 32
      srf(i,j,ns)=0.0
      go to 33
   32 srf(i,j,ns)=srf(i,j,ns)-dl
      fx=dmax1(fx,srf(i,j,ns))
   33 continue
   34 continue
      if(it.gt.0) write(it,40) ns,mm(1),mm(2),fx
   35 continue
      ll=ll+1
      go to 2
   36 continue
      if(it.gt.0) write(it,37) nc,ns
      return
      entry printg(nal)
      it=nal
      return
   37 format(/,' ',i3,' curves and',i3,' surfaces.'/)
   38 format(' plot: convex hull too large. increase ngs to',i6)
   39 format('   crv',i3,':  x(',i2,').  max =',g12.4)
   40 format('   srf',i3,':  x(',i2,'), x(',i2,').  max =',g12.4)
      end
      subroutine ctprt1 (m,nk,kp,kv,tb,cm,tc,sc,js)
      implicit none
      integer m,nk,kp(5,*),kv(2,*),js(*)
      double precision tb,cm(*),tc(*),sc(*)
      integer it,nc,j,jj,nl,i,k,ncat,nv,j1,j2,n1,n2,na,nb,ja,jb,nal
      double precision big,xm,xx,px,rx,rxp,s1,x
      double precision cvlv
      data big,it /9.9e30,6/
      if(it.le.0) return
      nc=ncat(kp)
      if(nc.eq.0) return
      write(it,'(/,'' there are'',i3,'' purely categorical basis functio
     1ns.'')') nc
      write(it,'('' purely additive and bivariate contributions follow''
     1)')
      if(m .ne. 1) go to 1
      write(it,'('' (piecewise-linear fit):'')')
      go to 2
    1 continue
      write(it,'('' (piecewise-cubic fit):'')')
    2 call catv(1,kp,kv,nv,js)
      do 8 jj=1,nv
      j=js(jj)
      xm=big
      xx=-big
      nl=int(cm(2*j+1)+.1)-int(cm(2*j)+.1)+1
      do 3 i=1,nl
      sc(i)=cvlv(m,1,j,i,nk,kp,kv,tb,cm,tc)
      xm=dmin1(xm,sc(i))
      xx=dmax1(xx,sc(i))
    3 continue
      px=99.0
      if(nl.gt.26) px=9.0
      rx=xx-xm
      rxp=rx/px
      write(it,'(/,'' f( x('',i3,'') ) : scale ='',g12.4)') j,rxp
      if(rxp.le.0.0) go to 8
      do 4 i=1,nl
      js(i+nv)=(sc(i)-xm)/rxp+.5
    4 continue
      if(nl .gt. 26) go to 5
      write(it,28) (i,i=1,nl)
      write(it,28) (js(i+nv),i=1,nl)
      go to 8
    5 if(nl .gt. 38) go to 6
      write(it,29) (i,i=1,nl)
      write(it,29) (js(i+nv),i=1,nl)
      go to 8
    6 if(nl .gt. 78) go to 7
      write(it,30) (mod(i,10),i=1,nl)
      write(it,30) (js(i+nv),i=1,nl)
      go to 8
    7 continue
      write(it,37) 78
    8 continue
      call catv(2,kp,kv,nv,js)
      do 27 jj=1,nv
      j1=js(2*jj-1)
      j2=js(2*jj)
      xm=big
      xx=-big
      n1=int(cm(2*j1+1)+.1)-int(cm(2*j1)+.1)+1
      n2=int(cm(2*j2+1)+.1)-int(cm(2*j2)+.1)+1
      k=0
      do 10 i=1,n1
      s1=cvlv(m,1,j1,i,nk,kp,kv,tb,cm,tc)
      js(2*nv+1)=i
      do 9 j=1,n2
      js(2*nv+2)=j
      k=k+1
      sc(k)=s1+cvlv(m,2,js(2*jj-1),js(2*nv+1),nk,kp,kv,tb,cm,tc)
    9 continue
   10 continue
      do 12 j=1,n2
      s1=cvlv(m,1,j2,j,nk,kp,kv,tb,cm,tc)
      do 11 i=1,n1
      k=j+n2*(i-1)
      sc(k)=sc(k)+s1
   11 continue
   12 continue
      k=0
      do 14 i=1,n1
      do 13 j=1,n2
      k=k+1
      x=sc(k)
      xx=dmax1(xx,x)
      xm=dmin1(xm,x)
   13 continue
   14 continue
      na=min0(n1,n2)
      nb=max0(n1,n2)
      if(na .ne. n1) go to 15
      ja=j1
      jb=j2
      go to 16
   15 ja=j2
      jb=j1
   16 px=99.0
      if(na.gt.25) px=9.0
      rx=xx-xm
      rxp=rx/px
      write(it,'(/,'' f( x('',i3,''), x('',i3,'') ) : scale ='',g12.4)')
     1ja,jb,rxp
      if(rxp.le.0.0) go to 27
      if(na .le. 75) go to 17
      write(it,37) 75
      go to 27
   17 if(na .gt. 25) go to 18
      write(it,34) (i,i=1,na)
      go to 20
   18 if(na .gt. 37) go to 19
      write(it,35) (i,i=1,na)
      go to 20
   19 continue
      write(it,36) (mod(i,10),i=1,na)
   20 do 26 j=1,nb
      do 23 i=1,na
      if(na .ne. n1) go to 21
      k=j+n2*(i-1)
      go to 22
   21 k=i+n2*(j-1)
   22 js(i+2*nv)=(sc(k)-xm)/rxp+.5
   23 continue
      if(na .gt. 25) go to 24
      write(it,31) j,(js(i+2*nv),i=1,na)
      go to 26
   24 if(na .gt. 37) go to 25
      write(it,32) j,(js(i+2*nv),i=1,na)
      go to 26
   25 continue
      write(it,33) j,(js(i+2*nv),i=1,na)
   26 continue
   27 continue
      return
      entry printc(nal)
      it=nal
      return
   28 format(' ',26i3)
   29 format(' ',38i2)
   30 format(' ',78i1)
   31 format(' ',i3,' ',25i3)
   32 format(' ',i3,' ',37i2)
   33 format(' ',i3,' ',75i1)
   34 format('     ',25i3)
   35 format('     ',37i2)
   36 format('     ',75i1)
   37 format(' function not printed (more than',i3,' categorical values)
     1.')
      end
      subroutine slice1 (flg,xs,n,p,x,nk,az,tb,cm,kp,kv,lp,lv,bz,tc,azn,
     1tbn,kpn,kvn,  lpn,lvn,bzn,tcn,sp,mm)
      implicit none
      integer n,p,nk,kp(5,*),kv(2,*),lp(3,*),lv(*),kpn(5,*),kvn(2,*),
     1lpn(3,*),lvn(*),mm(*)
      double precision flg,az,xs(p),x(n,p),tb(5,nk),cm(*),tc(*),tbn(5,nk
     1),tcn(*),sp(*),bz,azn,bzn
      integer it,ni,m,i,j,ig,i1,i2,i3
      double precision big,xl,xr
      data it,big /6,9.9e30/
      ni=0
      do 1 m=1,nk
      if(tb(1,m).ne.0.0) ni=ni+1
    1 continue
      if(ni .ne. 0) go to 3
      kpn(1,1)=-1
      lpn(1,1)=0
      do 2 m=1,nk
      tbn(1,m)=0.0
    2 continue
      azn=0.0
      bzn=azn
c     if(it.gt.0) write(it,'('' slice: original mars model = constant.''
c    1)')
      return
    3 if(it .le. 0) go to 5
c     write(it,'(/,'' sliced mars model: flag ='',g12.4)') flg
c     write(it,'(/,'' slice:'')')
      do 4 j=1,p
      if(xs(j).eq.flg) go to 4
c     write(it,'('' x('',i3,'') ='',g12.4)') j,xs(j)
    4 continue
    5 i1=2*nk+1
      i2=i1+2*p
      i3=max0(i2+p,i1+nk)
      do 7 j=1,p
      xl=big
      xr=-xl
      do 6 i=1,n
      xl=dmin1(xl,x(i,j))
      xr=dmax1(xr,x(i,j))
    6 continue
      sp(j+i3-1)=xr-xl
      sp(j+i3-1+p)=xl
    7 continue
      call reducq(flg,xs,nk,tb,cm,tc,kp,kv,lp,lv,sp(i3),sp,sp(i1),sp(i2)
     1)
      call reducl(flg,xs,nk,az,tb,cm,bz,sp,sp(i3),azn,tbn,bzn,sp(i1))
      ni=0
      do 8 m=1,nk
      if(tbn(1,m).ne.0.0) ni=ni+1
    8 continue
      if(ni .ne. 0) go to 10
      kpn(1,1)=-1
      lpn(1,1)=0
      do 9 m=1,nk
      tbn(1,m)=0.0
    9 continue
      azn=0.0
      bzn=azn
      if(it.gt.0) write(it,'('' sliced mars model = constant.'')')
      return
   10 if(it.gt.0) call slova(nk,it,tbn,ni,lpn,lvn)
      call ccoll(nk,tbn,cm,kpn,kvn,lpn,lvn,mm)
      call qslice(p,nk,tbn,cm,sp,kpn,kvn,lpn,lvn,tcn,sp(i3),sp(i1),mm)
      return
      entry prtslc(ig)
      it=ig
      return
      end
      subroutine cmrs (n,x,cm,kp,kv,lp,lv,bz,tc,y,sc)
      implicit none
      integer n,kp(5,*),kv(2,*),lp(3,*),lv(*)
      double precision bz,x(n,*),tc(*),cm(*),y(n),sc(n,2)
      integer ifg,la,ll,i,j,jl,jj,k,kk,ic,il,kp3,m,l,nt,lb,i1,l1
      integer icat
      data ifg /0/
      do 1 i=1,n
      y(i)=bz
    1 continue
      ll=1
      la=ll
      l1=la
    2 if(kp(1,ll).lt.0) go to 19
      do 3 i=1,n
      sc(i,1)=1.0
    3 continue
      if(kp(1,ll) .le. 0) go to 11
      jl=kp(1,ll)
      do 10 il=1,jl
      k=kp(2,ll)+il-1
      jj=kv(1,k)
      j=iabs(jj)
      kk=kv(2,k)
      do 9 i=1,n
      if(sc(i,1).eq.0.0) go to 9
      if(ifg .ne. 0) go to 4
      ic=icat(x(i,j),j,cm)
      go to 5
    4 ic=x(i,j)+.1
    5 if(ic .ne. 0) go to 6
      sc(i,1)=0.0
      go to 7
    6 sc(i,1)=cm(ic+kk)
    7 if(jj .ge. 0) go to 9
      if(sc(i,1) .ne. 0.0) go to 8
      sc(i,1)=1.0
      go to 9
    8 sc(i,1)=0.0
    9 continue
   10 continue
      go to 12
   11 if(kp(3,ll) .gt. 0) go to 12
      ll=ll+1
      go to 2
   12 if(kp(3,ll) .ge. 0) go to 14
      k=-kp(3,ll)
      ll=ll+1
      do 13 i=1,n
      if(sc(i,1).eq.0.0) go to 13
      y(i)=y(i)+tc(k)
   13 continue
      go to 2
   14 kp3=kp(3,ll)
      do 18 m=1,kp3
      l=lp(1,l1)
      nt=lp(3,l1)
      lb=la+5*l*nt-1
      do 17 j=1,nt
      do 15 i=1,n
      sc(i,2)=sc(i,1)
   15 continue
      call que(j,l,nt,lv(lp(2,l1)),n,x,tc(la),sc(1,2))
      do 16 i=1,n
      y(i)=y(i)+tc(lb+j)*sc(i,2)
   16 continue
   17 continue
      la=lb+nt+1
      l1=l1+1
   18 continue
      ll=ll+1
      go to 2
   19 return
      entry stcmrs(i1)
      ifg=i1
      return
      end
      subroutine fmrs (n,x,nk,az,tb,cm,y)
      implicit none
      integer n,nk
      double precision az,x(n,*),tb(5,nk),cm(*),y(n)
      double precision s
      integer ifg,i,m,ip,j,k,i1
      double precision phi,t,u
      integer icat
      data ifg /0/
      do 13 i=1,n
      s=az
      do 12 m=1,nk
      if(tb(1,m).eq.0.0) go to 12
      phi=1.0
      ip=m
    1 if(ip.le.0) go to 11
      t=tb(2,ip)
      j=abs(t)+.1
      if(cm(2*j) .le. 0.0) go to 8
      if(ifg .ne. 0) go to 2
      k=icat(x(i,j),j,cm)
      go to 3
    2 k=x(i,j)+.1
    3 if(k .ne. 0) go to 4
      u=0.0
      go to 5
    4 u=cm(k+int(tb(3,ip)+.1))
    5 if(t .ge. 0.0) go to 9
      if(u .ne. 0.0) go to 6
      u=1.0
      go to 9
    6 u=0.0
      go to 9
    8 u=dmax1(0.0D0,sign(1.0D0,t)*(x(i,j)-tb(3,ip)))
    9 if(u .ne. 0.0) go to 10
      phi=0.0
      go to 11
   10 phi=phi*u
      ip=tb(4,ip)+.1
      go to 1
   11 s=s+tb(1,m)*phi
   12 continue
      y(i)=s
   13 continue
      return
      entry stfmrs(i1)
      ifg=i1
      return
      end
      subroutine marsgo (n,p,x,y,w,nk,ms,df,fv,mi,lx,it,xm,xs,az,tb,cm,s
     1c,db,d,mm)
      implicit none
      integer n,p,nk,ms,mi,it,mm(n,*),lx(p), ii, jj
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

      if(it.gt.0) write(it,97)
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
      if(it.gt.0) write(it,'('' no predictor variables.'')')
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
      if(it.gt.0) write(it,98) m,txm,0.0,1.0
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
    9 continue
      call nxtpar(l,jq)
c      print *, 'parent, jq=', l, jq
      if(l.lt.0) go to 53
c      print *, 'before marsgo1 parent, nextVar=', l, jq
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
   14 continue
c     print *, 'begin mrsgo1'
      if(jq .ne. 0) go to 15
      jd1=1
      jd2=p
      go to 16
   15 jd1=jq
      jd2=jd1
   16 continue
c      print *, 'marsgo1 jd1, jd2=',jd1, jd2
c      print *, 'entering marsgo1'
      do 52 jp=jd1,jd2
cc      print *, 'marsgo1(findBestChildVariable) variable=',jp
      if(x(mm(1,jp),jp).ge.x(mm(n,jp),jp)) go to 52
c      print *,'jf, l, jp=', jf(l,jp,tb), l, jp
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
      call csp(jp,nc,m,n,x,y,w,nk,tb,cm,kcp,yb,d,kr,nnt,  sw,me,mkp2,nop
     1,sc(1,mkp1),db,d(1,3),mm(1,p+1))
c      print *, 'after csp, kcp, nc=', kcp, nc
c      print *, (cm(ii), ii=1,kcp+nc)
      if(nop.eq.0) go to 52
      go to 45
   20 continue
      tb(2,m)=jp
      tb(3,m)=x(mm(1,jp),jp)
      tb(4,m)=l
      k1=kr
      ssq=rsq
      call update(1,n,m,kr,x,y,w,sw,yb,tb,cm,sc,sc(1,mkp1),db,d,d(1,3))
c      print *, 'k1, kr, rsq, DY[kr]=', k1, kr, rsq, d(kr,1)
      if(kr .le. k1) go to 21
c      print *, 'k1, kr, rsq, DY[kr]=', k1, kr, rsq, d(kr,1)
      rsq=rsq-d(kr,1)**2
      tb(1,m)=rsq/sw
c      print *, 'rsq=', rsq
      go to 22
   21 tb(1,m)=big
   22 if((lx(jp) .ne. 3) .and. ((m .lt. mk) .and. (nnt .gt. me+mel))) go
     1 to 26
      tb(1,m)=rsq/sw
      newbf=newb(m,tb).eq.0
      if(fvr*tb(1,m) .gt. txl .or. .not.(newbf)) go to 23
      txl=fvr*tb(1,m)
      tx1=tb(1,m)
      jq=jp
   23 if(fvr*tb(1,m) .gt. txi .or. .not.(newbf)) go to 25
      txi=fvr*tb(1,m)
      tx(1)=tb(1,m)
      do 24 i=2,4
      tx(i)=tb(i,m)
   24 continue
      jas=ja
   25 kr=k1
      rsq=ssq
      go to 52
   26 mm1=m
c      print *, 'after 26, mm1=', mm1
      m=m+1
      tb(1,m)=big
      xa=0.0
c getPotentialKnots
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
cc      print *, 'b, coef =',b, tb(1,m)
c     include a small amount for error
      if(b+1.e-10 .ge. tb(1,m)) go to 44
      tb(1,m)=b
      tb(3,m)=xa
cc      print *,'j, j0, coef, knot=', j, j0, b, xa
   44 j=j-1
      if(j.le.1) go to 45
      go to 27
c end findBestKnot
   45 continue
      tb(2,m)=jp
      tb(4,m)=l
      tb(1,m)=(rsq+tb(1,m))/sw
cc      print *, 'proposedHS= ', (tb(jj,m), jj=1,4)
c      print *, 'mm1=', mm1, m, nop
      if (ict .ne. 0) go to 46
      if (tb(1,mm1) .gt. fln*tb(1,m)) go to 46
      mp=mm1
      go to 47
   46 mp=m
   47 newbf=newb(mp,tb).eq.0
c      print *, 'in decide'
c      print *, 'hockeyStick=', (tb(jj,mp),jj=1,4)
c      print *, 'fvr, tb(1,mp), txl=', fvr, tb(1,mp), txl
      if(fvr*tb(1,mp) + 1.e-10 .ge. txl .or. .not.(newbf)) go to 48
c      print *, 'pass 1'
      txl=fvr*tb(1,mp)
      tx1=tb(1,mp)
      jq=jp
   48 continue
c      print *, 'fvr, tb(1,mp), txi=', fvr, tb(1,mp), txi
      if(fvr*tb(1,mp) + 1.e-10 .ge. txi .or. .not.(newbf)) go to 51
c      print *, 'pass 2'
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
      print *, 'in decideNextVariable:kcp, kcp0, nc =', kcp, kcp0, nc
      print *, (cm(kcp0+ii), ii=1,nc)
      print *, 'leaving decide'
c     end decideNextVariable   
   51 if(ict .ne. 0) go to 52
      m=mm1
      mm1=m-1
      kr=k1
      rsq=ssq
   52 continue
      call updpar(jq,asq0-tx1)
c      print *, 'after mrsgo1 bfIndex, rsq=', m, rsq
c      print *, 'after marsgo1 nextVar, asq0-tx1=', jq, asq0-tx1
cc      print *, 'after marsgo1 newBF=', (tx(jj), jj=1,5)
c      print *,'after mrsgo1(findBestChildVar) tb=',(tx(jj),jj=1,4)
      go to 9
   53 continue
c      print *,'after mrsgo1(findBestChildVa) tb=',(tb(3,jj),jj=1,nk)
c      print *,'after mrsgo1(findBestChildVa) tb=',(tb(jj,m),jj=1,5)
      jp=tx(2)+.1
c      print *, 'jp=', jp
c      print *,'after mrsgo1 tb=',(tb(1,jj),jj=1,nk)
      call selpar(int(tx(4)+.1))
c      print *, '    jp, cm(2*jp)=', jp, cm(2*jp)
      if(cm(2*jp) .le. 0.) go to 54
      nc=int(cm(2*jp+1)+.1)-int(cm(2*jp)+.1)+1
      kcp0=kcp0+nc
   54 if(jas .le. 0) go to 60
      print *, 'nested data'
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
      if(m.eq.mtot+1) write(it,100) m,gcv,fkr,tcst,fjn,hol,tb(4,m)
      if(m.eq.mtot+2) write(it,99) m,mp,gcv,fkr,tcst,fjn,hol,tb(4,m)
   59 mtot=m
      m=m+1
      if(m.gt.mk) go to 69
   60 do 61 i=1,5
      tb(i,m)=tx(i)
   61 continue
c      print *, 'after addBasisFunction, tb=', (tb(ii,m), ii=1,4)
      k1=kr
      call blf(int(tx(4)+.1),n,sc,sc(1,mkp1))
      call update(2,n,m,kr,x,y,w,sw,yb,tb,cm,sc,sc(1,mkp1),db,d,d(1,3))
      if(kr.gt.k1) rsq=rsq-d(kr,1)**2
      call addpar(m)
c      print *, 'm, mk, kr, kl,tx[3], x[mm[1,jp],jp]=',m, mk, kr, 
c     1kl,tx(3), x(mm(1,jp),jp)
      if(m .ge. mk .or. (cm(2*jp) .le. 0.0) .and. (tx(3) .le. x(mm(1,jp)
     1,jp))) go to 66
      m=m+1
      do 62 i=1,4
      tb(i,m)=tx(i)
   62 continue
      tb(2,m)=-tb(2,m)
c      print *,'after addReflectedBasisFunction, tb=', (tb(ii,m), ii=1,4)
      if(cm(2*jp) .le. 0.0) go to 64
      do 63 i=1,n
      sc(i,m)=phi(m,i,n,x,tb,cm)
   63 continue
      go to 65
   64 k1=kr
      call update(2,n,m,kr,x,y,w,sw,yb,tb,cm,sc,sc(1,mkp1),db,d,d(1,3))
      if(kr.gt.k1) rsq=rsq-d(kr,1)**2
   65 call addpar(m)
   66 continue
cc      print *,'after addBasisFunction knots=',(tb(3,jj),jj=1,nk)
cc      print *,'after addBasisFunction bf=',(tb(jj,m),jj=1,5)
      tcst=nopt*df1+kr+1.0
      if(it .le. 0) go to 68
      mp=m-1
      jp=abs(tx(2))+.1
      fkr=kr
      gcv=(rsq/sw)/(1.d0-tcst/wn)**2
c      print *,'after marsgo1 m, rsq, sw, tcst, wn, gcv=',m, rsq, sw, 
c     1tcst, wn, gcv
      if(cm(2*jp) .le. 0.0) go to 67
      call holl(jp,cm,tx(3),hol)
      print *, 'hol=', hol
      if(m.eq.mtot+1) write(it,100) m,gcv,idint(fkr+.5),idint(tcst+.5),
     1idint(tx(2)+.5),hol,idint(tx(4)+.5)
      if(m.eq.mtot+2) write(it,99) m,mp,gcv,idint(fkr+.5),idint(tcst+.5)
     1,idint(tx(2)+.5),hol,idint(tx(4)+.5)
      go to 68
   67 xk=xm(jp)+xs(jp)*tx(3)
c      print *, 'xk, xm, xs, tx, m', xk, xm(jp), xs(jp), tx(3), m
      if(m.eq.mtot+1) write(it,93) m,gcv,idint(fkr+.5),idint(tcst+.5),
     1idint(tx(2)+.5),xk,idint(tx(4)+.5)
      if(m.eq.mtot+2) write(it,94) m,mp,gcv,idint(fkr+.5),idint(tcst+.5)
     1,idint(tx(2)+.5),xk,idint(tx(4)+.5)
   68 mtot=m
      go to 8
   69 continue
      mk=min0(m,mk)
      m=mk+1
      k=m
      go to 71
   70 k=k+1
   71 if((k).gt.(nk)) go to 72
      tb(1,k)=0.0
      go to 70
   72 call sscp(n,m,sc,y,w,yb,yv,sw,db,d)
c      print *, 'd='
c      do 998 ii=1,mk
c 998  print *, (d(jj,ii), jj=1,mk)
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
   81 continue
      do 88 ll=2,nli
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
      write(it,95)
c      print *, 'xs=', xs
c      do 999 ii=1,mk
c 999  print *, (tb(jj,ii), jj=1,4)
      call coefpr(it,mk,az,tb,cm,xs)
      write(it,96) asm,tcsts
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
   93 format('   ',i3,'    ',g12.4,2('   ',i7),'       ',i7,  '  ',
     1g12.4,'  ',i7)
   94 format(' ',i3,' ',i3,'  ',g12.4,2('   ',i7),'       ',i7,  '  ',
     1g12.4,'  ',i7)
   95 format(/,' final model after backward stepwise elimination:')
   96 format(/,'   (piecewise linear) gcv = ',g12.4,'   #efprms = ',f5.1
     1)
   97 format(//,' forward stepwise knot placement:',//  '  basfn(s)    g
     1cv      #indbsfns  #efprms',  '   variable      knot            pa
     1rent')
   98 format('   ',i3,'    ',g12.4,2('   ',f5.1))
   99 format(' ',i3,' ',i3,'  ',g12.4,2('   ',i7),'       ',i7,a28,i3)
  100 format('   ',i3,'    ',g12.4,2('   ',i7),'       ',i7,a28,i7)
      end
c     end marsgo
      subroutine addpar (ib)
      implicit none
      integer ib
      integer maxdph
      parameter(maxdph=2000)
      double precision que(2,maxdph),sp(maxdph)
      integer m(maxdph),n(maxdph),jp(2,maxdph)
      double precision val
      integer lq,kp,itr,ktr,i,j,k,l,iarg,jq,mpr,mtr,jj
      double precision big,arg,beta
      save mpr,mtr,ktr,big,beta,lq,kp,itr,jp,que,n,m
      data big,mpr,mtr,beta /9.9e30,10,5,1.0/
      if(ib .ne. 0) go to 1
      lq=1
      que(1,1)=big
      que(2,1)=0.0
      m(1)=1
      kp=0
      itr=kp
      n(1)=itr
      jp(1,1)=n(1)
      jp(2,1)=jp(1,1)
      ktr=0.5*(mpr-1)+.1
      return
    1 continue
c      print *,'   begin addpar m=', (m(jj), jj=1,lq)
      if(que(1,lq).ge.-0.5) go to 2
      lq=lq-1
      go to 1
    2 i=1
      go to 4
    3 i=i+1
    4 if((i).gt.(lq)) go to 7
      if(que(1,i).ge.-0.5) go to 3
      lq=lq-1
      do 6 j=i,lq
      n(j)=n(j+1)
      do 5 k=1,2
      jp(k,j)=jp(k,j+1)
      que(k,j)=que(k,j+1)
    5 continue
    6 continue
      i=i-1
      go to 3
    7 lq=lq+1
      if(lq .le. maxdph) go to 8
      write(6, '('' increase parameter maxdph in subroutine addpar to ''
     1)')
      write(6,'('' '',i10,''  or larger, and recompile mars.'')') lq
      stop
    8 que(1,lq)=big
      que(2,lq)=0.0
      n(lq)=ib
      jp(1,lq)=0
      jp(2,lq)=jp(1,lq)
      do 9 i=1,lq
      m(i)=i
      sp(i)=que(1,i)
    9 continue
      call psort(sp,m,1,lq)
c      print *,'   middle addpar m=', (m(jj), jj=1,lq)
      do 10 i=1,lq
      j=m(i)
      sp(j)=i+beta*(itr-que(2,j))
   10 continue
c      print *,'   middle addpar sp=', (sp(jj), jj=1,lq)
      call psort(sp,m,1,lq)
      kp=max0(0,lq-mpr)
c      print *,'   end addpar ib, kp=', ib, kp
c      print *,'   end addpar lq=', lq
c      print *,'   end addpar que=', (que(1,jj), jj=1,lq)
c      print *,'   end addpar m=', (m(jj), jj=1,lq)
      return
      entry nxtpar (l,jq)
      kp=kp+1
      if(kp .le. lq) go to 11
      l=-1
c      print *,'   in nxtpar parent jq, kp, lq=', l, jq, kp, lq
      return
   11 l=n(m(kp))
      if(itr-jp(2,m(kp)).gt.mtr.or.itr.le.ktr) jp(1,m(kp))=0
      jq=jp(1,m(kp))
c      print *, '   in nxtpar kp, lq=', kp, lq
c      print *, '   in nxtpar m=', (m(jj), jj=1,kp)
c      print *, '   in nxtpar n=', (n(jj), jj=1,m(kp))
c      print *, '   in nxtpar jp=', (jp(1,jj), jj=1,m(kp))
c      print *, '   in nxtpar parent jq=', l, jq
      return
      entry updpar (jq,val)
      que(1,m(kp))=val
      que(2,m(kp))=itr
      if(jp(1,m(kp)) .ne. 0) go to 12
      jp(1,m(kp))=jq
      jp(2,m(kp))=itr
   12 continue
c      print *, '   in updpar jq, val=', jq, val
c      print *, '   in updpar kp, m=', kp, (m(jj), jj=1,kp)
c      print *, '   in updpar jp=', (jp(1,jj), jj=1,m(kp))
      return
      entry selpar(ib)
c      print *,'    in selpar, ib, lq', ib, lq
      do 13 i=lq,1,-1
      if(n(i).ne.ib) go to 13
      jp(1,i)=0
      go to 14
   13 continue
   14 return
      entry itrpar(iarg)
      itr=iarg
      return
      entry setmpr(iarg)
      mpr=iarg
      return
      entry setbta(arg)
      beta=arg
      return
      entry setfrq(arg)
      mtr=1.0/dmax1(arg,0.01D0)+.1
      return
      end
      subroutine speed(is)
      implicit none
      integer is
      integer j
      integer lque(5)
      double precision freq(5)
      save lque,freq
      data lque /9999,20,20,10,5/
      data freq /9.e30,9.e30,0.2,0.2,0.2/
      j=is
      if(is.lt.1) j=1
      if(is.gt.5) j=5
      call setmpr(lque(j))
      call setfrq(freq(j))
      return
      end
      subroutine atoscl(n,p,w,x,lx,mm,xm,xs,cm,z)
      implicit none
      integer n,p,lx(p),mm(n,p)
      double precision w(n),x(n,p),z(n,p),xm(p),xs(p),cm(*)
      integer i,j,k,ip,nct,nc,j0,n2,n2p1
      double precision s,t,sw
      sw=0.d0
      do 1 j=1,n
      sw=sw+w(j)
    1 continue
      ip=0
      nct=ip
      do 12 i=1,p
      if(lx(i) .ne. 0) go to 2
      xm(i)=0.0
      xs(i)=xm(i)
      go to 12
    2 if(lx(i) .ge. 0) go to 8
      nc=0
      xm(i)=ip
      j=1
      nct=nct+1
    3 j0=j
      if(j .ge. n) go to 5
    4 if(x(mm(j+1,i),i).gt.x(mm(j,i),i)) go to 5
      j=j+1
      if(j.ge.n) go to 5
      go to 4
    5 ip=ip+1
      cm(ip)=x(mm(j,i),i)
      nc=nc+1
      do 6 k=j0,j
      z(mm(k,i),i)=nc
    6 continue
      j=j+1
      if(j.gt.n) go to 7
      go to 3
    7 xs(i)=nc
      go to 12
    8 s=0.d0
      t=s
      do 9 j=1,n
      s=s+w(j)*x(j,i)
    9 continue
      s=s/sw
      xm(i)=s
      do 10 j=1,n
      z(j,i)=x(j,i)-s
      t=t+w(j)*z(j,i)**2
   10 continue
      xs(i)=1.0
      if(t.le.0.d0) go to 12
      t=dsqrt(t/sw)
      xs(i)=t
      t=1.d0/t
      do 11 j=1,n
      z(j,i)=t*z(j,i)
   11 continue
   12 continue
      n2=2*p+1
      if(nct .ne. 0) go to 14
      do 13 i=1,n2
      cm(i)=0.0
   13 continue
      return
   14 n2p1=n2+1
      i=ip
      go to 16
   15 i=i+(-1)
   16 if((-1)*((i)-(1)).gt.0) go to 17
      cm(i+n2)=cm(i)
      go to 15
   17 j=0
      i=2
      go to 19
   18 i=i+(2)
   19 if((2)*((i)-(n2)).gt.0) go to 22
      j=j+1
      if(lx(j) .ge. 0) go to 20
      cm(i)=xm(j)+n2p1
      cm(i+1)=cm(i)+xs(j)-1.0
      go to 18
   20 cm(i)=0.0
      cm(i+1)=cm(i)
      go to 18
   22 cm(1)=nct
      call stfmrs(1)
      call stcmrs(1)
      return
      end
      subroutine orgpl (xm,xs,nk,tb,cm)
      implicit none
      integer nk
      double precision xm(*),xs(*),tb(5,nk),cm(*)
      integer m,j,ip
      double precision scl
      do 1 m=1,nk
      j=abs(tb(2,m))+.1
      if(cm(2*j).gt.0.0) go to 1
      tb(3,m)=xm(j)+xs(j)*tb(3,m)
    1 continue
      do 4 m=1,nk
      if(tb(1,m).eq.0.0) go to 4
      scl=1.0
      ip=m
    2 if(ip.le.0) go to 3
      j=abs(tb(2,ip))+.1
      if(cm(2*j).eq.0.0) scl=scl*xs(j)
      ip=tb(4,ip)+.1
      go to 2
    3 tb(1,m)=tb(1,m)/scl
    4 continue
      return
      end
      subroutine anova (n,x,y,w,nk,it,tb,cm,lp,lv,t,d)
      implicit none
      integer n,nk,it,lp(3,*),lv(*)
      double precision x(n,*),y(n),w(n),tb(5,nk),cm(*),t(n,nk)
      double precision d(nk,*),s,u,sw,yv,wn,yb
      integer nkp1,nkp2,nkp3,lm,i,j,k,m,ni,na,nim1,k2,
     1i2,l,ll,im,np,jf
      double precision eft,efm,efp
      integer nord
      double precision phi,varf
      if(it.le.0) return
      nkp1=nk+1
      nkp2=nkp1+1
      nkp3=nkp2+1
      lm=nkp3+nk
      sw=0.d0
      wn=sw
      s=wn
      u=s
      do 1 i=1,n
      sw=sw+w(i)
      wn=wn+w(i)**2
      s=s+w(i)*y(i)
    1 continue
      s=s/sw
      yb=s
      wn=sw**2/wn
      do 2 i=1,n
      u=u+w(i)*(y(i)-s)**2
    2 continue
      yv=u/sw
      eft=1.0
      do 3 m=1,nk
      if(tb(1,m).ne.0.0) eft=eft+tb(5,m)
    3 continue
      ni=0
      do 9 m=1,nk
      if(tb(1,m).eq.0.0) go to 9
      ni=ni+1
      s=0.d0
      do 4 j=1,n
      t(j,ni)=phi(m,j,n,x,tb,cm)
      s=s+w(j)*t(j,ni)
    4 continue
      s=s/sw
      do 5 j=1,n
      t(j,ni)=t(j,ni)-s
    5 continue
      do 7 i=1,ni
      s=0.d0
      do 6 j=1,n
      s=s+w(j)*t(j,i)*t(j,ni)
    6 continue
      d(i,ni)=s
    7 continue
      s=0.d0
      do 8 j=1,n
      s=s+w(j)*t(j,ni)*(y(j)-yb)
    8 continue
      d(ni,nkp1)=s
      d(ni,nkp2)=tb(1,m)
      lv(ni)=m
    9 continue
      if(ni .ne. 0) go to 10
      write(it,26)
      return
   10 do 11 m=1,ni
      t(m,1)=lv(m)
   11 continue
      write(it,24) ni
      call coll(nk,tb,lp,lv,lp(1,nkp1))
      m=1
   12 if(lp(1,m).eq.0) go to 13
      m=m+1
      go to 12
   13 na=m-1
      m=1
      nim1=ni-1
      if(na .ne. 1) go to 14
      k2=lp(2,m)
      i2=lp(1,m)+k2-1
      efm=eft-1.0
      u=yv/(1.d0-1.d0/wn)**2
      s=dsqrt(varf(nk,d,d(1,nkp2),sw,1,ni))
      write(it,25) m,s,u,lp(3,m),efm,(lv(i),i=k2,i2)
      return
   14 do 23 m=1,na
      k2=lp(2,m)
      l=lp(1,m)
      i2=l+k2-1
      ll=k2-1
      np=ni
      do 19 im=1,ni
      i=t(im,1)+.1
      if(nord(i,tb) .eq. l) go to 15
      t(im,2)=0.0
      go to 19
   15 k=0
      do 16 j=1,l
      if(jf(i,lv(ll+j),tb).eq.1) go to 16
      k=1
      go to 17
   16 continue
   17 continue
      if(k .ne. 1) go to 18
      t(im,2)=0.0
      go to 19
   18 t(im,2)=1.0
      np=np-1
   19 continue
   20 k=0
      do 21 i=1,nim1
      if(t(i,2) .le. t(i+1,2)) go to 21
      k=1
      if (ni.lt.nkp2) then
        call exch(nk,nkp2,i,d,t,t(1,2))
      else
        call exch(nk,ni,i,d,t,t(1,2))
      end if
   21 continue
      if(k.eq.0) go to 22
      go to 20
   22 call lsf(nk,np,nkp1,0.d0,d,d(1,lm),s,u,d(1,nkp3),1)
      efm=efp(lp(1,m),lv(lp(2,m)),nk,tb)
      u=(u/sw+yv)/(1.d0-(eft-efm)/wn)**2
      s=dsqrt(varf(nk,d,d(1,nkp2),sw,np+1,ni))
      write(it,25) m,s,u,lp(3,m),efm,(lv(i),i=k2,i2)
   23 continue
      return
   24 format(/' anova decomposition on',i3,' basis functions:',/  '  fun
     1. std. dev.     -gcv    #bsfns  #efprms  variable(s)')
   25 format(' ',i3,' ',2g12.4,'  ',i2,'      ',f4.1,'  ',20i4)
   26 format(/' estimated optimal model = response mean.')
      end
      subroutine anoval (n,x,y,w,nk,il,it,az,tb,cm,lp,lv,sc,d)
      implicit none
      integer n,nk,il,it,lp(3,*),lv(*)
      double precision az,x(n,*),y(n),w(n),tb(5,nk),cm(*),sc(n,*)
      double precision d(nk,*),sw,yv,wn,yb
      integer i,ni,m,i2,k,k2,ip,na,l,ll,nord,j,jf, ii
      double precision eft,efm,u,a0,efp
      if(it.le.0) return
      sw=0.d0
      yb=sw
      yv=yb
      wn=yv
      do 1 i=1,n
      sw=sw+w(i)
      wn=wn+w(i)**2
      yb=yb+w(i)*y(i)
    1 continue
      yb=yb/sw
      wn=sw**2/wn
      do 2 i=1,n
      yv=yv+w(i)*(y(i)-yb)**2
    2 continue
      yv=yv/sw
      eft=1.0
      ni=0
      do 3 m=1,nk
      if(tb(1,m).eq.0.0) go to 3
      ni=ni+1
      eft=eft+tb(5,m)
    3 continue
      if(ni .ne. 0) go to 4
      write(it,14)
      return
    4 continue
      write(it,12) ni
      call coll(nk,tb,lp,lv,lp(1,nk+1))
      m=1
    5 if(lp(1,m).eq.0) go to 6
      m=m+1
      go to 5
    6 na=m-1
      m=1
      if(na .ne. 1) go to 7
      k2=lp(2,m)
      i2=lp(1,m)+k2-1
      efm=eft-1.0
      u=yv/(1.d0-1.d0/wn)**2
      write(it,13) m,u,lp(3,m),efm,(lv(i),i=k2,i2)
      return
    7 ip=nk+4
      do 11 m=1,na
      k2=lp(2,m)
      l=lp(1,m)
      i2=l+k2-1
      ll=k2-1
      call cptb(nk,tb,sc(1,ip))
      do 10 i=1,nk
      if(tb(1,i).eq.0.0) go to 10
      if(nord(i,tb).ne.l) go to 10
      k=0
      do 8 j=1,l
      if(jf(i,lv(ll+j),tb).eq.1) go to 8
      k=1
      go to 9
    8 continue
    9 if(k.eq.1) go to 10
      call setz(i,sc(1,ip))
   10 continue
      a0=az
      call vp(n,x,y,w,nk,il,yb,sw,a0,sc(1,ip),cm,u,sc,d)
      efm=efp(lp(1,m),lv(lp(2,m)),nk,tb)
      u=u/(1.d0-(eft-efm)/wn)**2
      write(it,13) m,u,lp(3,m),efm,(lv(i),i=k2,i2)
   11 continue
      return
   12 format(/' logit anova decomposition on',i3,' basis functions:',/
     1'  fun.    -gcv    #bsfns  #efprms  variable(s)')
   13 format(' ',i3,' ',g12.4,'   ',i2,'     ',f4.1,'    ',20i4)
   14 format(/' estimated optimal model = response mean.')
      end
      subroutine cptb(nk,tb,ub)
      implicit none
      integer nk
      double precision tb(5,nk),ub(5,*)
      integer k,m,l
      do 2 m=1,nk
      do 1 k=1,5
      ub(k,m)=tb(k,m)
    1 continue
    2 continue
      return
      entry setz(l,ub)
      ub(1,l)=0.0
      return
      end
      subroutine fun (l,jv,n,x,nk,tb,cm,jl,kv,t,js)
      implicit none
      integer l,n,nk,jl,jv(l),kv(2,jl),js(*)
      double precision x(n,l),tb(5,nk),cm(*),t(n)
      double precision s
      integer i,j,k,m,ip
      double precision phi,u
      integer icf,nordc,jf
      do 8 i=1,n
      s=0.d0
      do 7 m=1,nk
      if(icf(m,tb,cm,jl,kv,js).eq.0) go to 7
      if(nordc(1,m,tb,cm).ne.l) go to 7
      k=0
      do 1 j=1,l
      if(jf(m,jv(j),tb).eq.1) go to 1
      k=1
      go to 2
    1 continue
    2 if(k.eq.1) go to 7
      phi=1.0
      ip=m
    3 if(ip.le.0) go to 6
      u=tb(2,ip)
      j=abs(u)+.1
      if(cm(2*j) .eq. 0.0) go to 4
      ip=tb(4,ip)+.1
      go to 3
    4 do 5 k=1,l
      if(j.eq.jv(k)) j=k
    5 continue
      phi=phi*dmax1(0.0D0,sign(1.0D0,u)*(x(i,j)-tb(3,ip)))
      ip=tb(4,ip)+.1
      go to 3
    6 s=s+tb(1,m)*phi
    7 continue
      t(i)=s
    8 continue
      return
      end
      subroutine cubic (n,p,x,y,w,nk,it,tb,cm,kp,kv,lp,lv,bz,tc,t,z,sc,j
     1s,d)
      implicit none
      integer n,p,nk,it,kp(5,*),kv(2,*),lp(3,*),lv(*),js(*), ii
      double precision bz,x(n,p),y(n),w(n),tb(5,nk),cm(*),tc(*),t(n,nk),
     1z(2,p),sc(n)
      double precision d(nk,*),s,u,sw,yb,wn,yv
      integer i,j,k,ll,la,l1,lt,jl,il,jj,kk,ic,ni,m,nkp1,nkp2,nkp3,lm,
     1kp3,l,nt,jp,le
      double precision big,xl,xr,eft
      data big /9.9e30/
      yb=0.d0
      sw=yb
      wn=sw
      yv=wn
      do 1 i=1,n
      sw=sw+w(i)
      wn=wn+w(i)**2
      yb=yb+w(i)*y(i)
    1 continue
      yb=yb/sw
      wn=sw**2/wn
      do 2 i=1,n
      yv=yv+w(i)*(y(i)-yb)**2
    2 continue
      yv=yv/sw
      ni=0
      do 3 m=1,nk
      if(tb(1,m).ne.0.0) ni=ni+1
    3 continue
      if(ni .ne. 0) go to 4
      bz=yb
      u=yv/(1.0-1.0/wn)**2
c     if(it.gt.0) write(it,34) ni,u
      return
    4 nkp1=nk+1
      nkp2=nk+2
      nkp3=nk+3
      lm=nkp3+nk
      do 6 i=1,p
      xl=big
      xr=-xl
      do 5 j=1,n
      xl=dmin1(xl,x(j,i))
      xr=dmax1(xr,x(j,i))
    5 continue
      z(1,i)=xl
      z(2,i)=xr
    6 continue
      ll=1
      la=ll
      l1=la
      lt=0
    7 if(kp(1,ll).lt.0) go to 20
      do 8 i=1,n
      sc(i)=1.0
    8 continue
c      print *, 'in cubic  ll, FLAG1, FLAG2=', ll, kp(1,ll), kp(3,ll)
      if(kp(1,ll) .le. 0) go to 12
      jl=kp(1,ll)
      do 11 il=1,jl
      k=kp(2,ll)+il-1
      jj=kv(1,k)
      j=iabs(jj)
      kk=kv(2,k)
      do 10 i=1,n
      if(sc(i).eq.0.0) go to 10
      ic=x(i,j)+.1
      sc(i)=cm(ic+kk)
      if(jj .ge. 0) go to 10
      if(sc(i) .ne. 0.0) go to 9
      sc(i)=1.0
      go to 10
    9 sc(i)=0.0
   10 continue
   11 continue
      go to 13
   12 if(kp(3,ll) .gt. 0) go to 13
      ll=ll+1
      go to 7
   13 if(kp(3,ll) .gt. 0) go to 15
      lt=lt+1
      kp(5,ll)=0
      do 14 i=1,n
      t(i,lt)=sc(i)
   14 continue
      go to 19
   15 kp3=kp(3,ll)
      kp(5,ll)=la
      do 18 m=1,kp3
      l=lp(1,l1)
      nt=lp(3,l1)
      call knts(l,nt,lv(lp(2,l1)),kp(1,ll),kv(1,kp(2,ll)),  nk,tb,cm,tc(
     1la),js)
      call side(l,nt,lv(lp(2,l1)),z,tc(la))
      do 17 jp=1,nt
      lt=lt+1
      do 16 i=1,n
      t(i,lt)=sc(i)
   16 continue
      call que(jp,l,nt,lv(lp(2,l1)),n,x,tc(la),t(1,lt))
   17 continue
      l1=l1+1
      la=la+nt*(5*l+1)
   18 continue
   19 ll=ll+1
      go to 7
   20 continue
c      print *, 'lt=', lt
      do 26 j=1,lt
      s=0.d0
      u=s
      do 21 i=1,n
      s=s+w(i)*t(i,j)
   21 continue
      s=s/sw
      d(j,nkp2)=s
      do 22 i=1,n
      t(i,j)=t(i,j)-s
   22 continue
      s=0.d0
      do 23 i=1,n
      s=s+w(i)*(y(i)-yb)*t(i,j)
   23 continue
      d(j,nkp1)=s
      do 25 k=1,j
      s=0.d0
      do 24 i=1,n
      s=s+w(i)*t(i,k)*t(i,j)
   24 continue
      d(k,j)=s
   25 continue
c      print *, 'in cubic   j, d=', j, (d(ii,j), ii=1,j)
   26 continue
      call lsf(nk,lt,nkp1,yb,d,d(1,lm),s,u,d(1,nkp3),1)
      eft=1.0
      do 27 i=1,nk
      if(tb(1,i).ne.0.0) eft=eft+tb(5,i)
   27 continue
      u=(u/sw+yv)/(1.0-eft/wn)**2
      bz=s
      ll=1
      l1=ll
      le=la-1
      la=0
      lt=la
   28 if(kp(1,ll).lt.0) go to 33
      if(kp(1,ll) .ne. 0 .or. kp(3,ll) .gt. 0) go to 29
      ll=ll+1
      go to 28
   29 if(kp(3,ll) .gt. 0) go to 30
      le=le+1
      kp(3,ll)=-le
      lt=lt+1
      tc(le)=d(lt,lm)
      ll=ll+1
      go to 28
   30 kp3=kp(3,ll)
      do 32 m=1,kp3
      nt=lp(3,l1)
      la=la+5*lp(1,l1)*nt
      do 31 i=1,nt
      lt=lt+1
      tc(i+la)=d(lt,lm)
   31 continue
      la=la+nt
      l1=l1+1
   32 continue
      ll=ll+1
      go to 28
   33 continue
      if(it.gt.0) write(it,34) lt,u
      return
   34 format(/' piecewise cubic fit on',i3,' basis functions, gcv =',g12
     1.4)
      end
      subroutine cfun (l,jv,n,x,nf,lp,lv,tc,t,sc,jw)
      implicit none
      integer l,n,nf,jv(l),lp(3,*),lv(*),jw(l)
      double precision x(n,l),tc(*),t(n),sc(n)
      integer i,la,l2,j,m,k,nt,lb,l1
      do 1 i=1,n
      t(i)=0.0
    1 continue
      la=1
      do 10 l1=1,nf
      if(lp(1,l1).ne.l) go to 9
      l2=lp(2,l1)-1
      do 3 j=1,l
      m=0
      do 2 k=1,l
      if(jv(j).eq.lv(k+l2)) m=1
    2 continue
      if(m.eq.0) go to 9
    3 continue
      nt=lp(3,l1)
      lb=la+5*l*nt-1
      do 8 j=1,nt
      do 5 k=1,l
      do 4 i=1,l
      if(lv(k+l2).eq.jv(i)) jw(k)=i
    4 continue
    5 continue
      do 6 i=1,n
      sc(i)=1.0
    6 continue
      call que(j,l,nt,jw,n,x,tc(la),sc)
      do 7 i=1,n
      t(i)=t(i)+tc(lb+j)*sc(i)
    7 continue
    8 continue
      go to 11
    9 la=la+lp(3,l1)*(5*lp(1,l1)+1)
   10 continue
   11 return
      end
      subroutine orgpc (xm,xs,lp,lv,tc)
      implicit none
      integer lp(3,*),lv(*)
      double precision xm(*),xs(*),tc(*)
      integer la,l1,l,nt,lb,j
      la=1
      l1=la
    1 if(lp(1,l1).eq.0) go to 3
      l=lp(1,l1)
      nt=lp(3,l1)
      lb=la+5*l*nt-1
      do 2 j=1,nt
      call scpc(xm,xs,j,l,nt,lv(lp(2,l1)),tc(la),tc(lb+j))
    2 continue
      la=lb+nt+1
      l1=l1+1
      go to 1
    3 return
      end
      subroutine pair (jv,n,x,nk,tb,cm,jl,kv,f,sc,js)
      implicit none
      integer n,nk,jl,jv(2),kv(2,jl),js(*)
      double precision x(n,*),tb(5,nk),cm(*),f(n),sc(n)
      integer i,k
      call fun(2,jv,n,x,nk,tb,cm,jl,kv,f,js)
      do 2 k=1,2
      call fun(1,jv(k),n,x(1,k),nk,tb,cm,jl,kv,sc,js)
      do 1 i=1,n
      f(i)=f(i)+sc(i)
    1 continue
    2 continue
      return
      end
      subroutine cpair (jv,n,x,nf,lp,lv,tc,f,sc)
      implicit none
      integer n,nf,jv(2),lp(3,*),lv(*),jw(2)
      double precision x(n,*),tc(*),f(n),sc(n,2)
      integer i,k
      call cfun(2,jv,n,x,nf,lp,lv,tc,f,sc,jw)
      do 2 k=1,2
      call cfun(1,jv(k),n,x(1,k),nf,lp,lv,tc,sc,sc(1,2),jw)
      do 1 i=1,n
      f(i)=f(i)+sc(i,1)
    1 continue
    2 continue
      return
      end
      subroutine logitl (n,x,y,w,nk,il,az,tb,cm,sc,d)
      implicit none
      integer n,nk,il,kp(5,*),kv(2,*),lp(3,*),lv(*)
      double precision az,x(n,*),y(n),w(n),tb(5,nk),cm(*),sc(n,*),tc(*),
     1ss(n)
      double precision d(nk,*),a,b,s,sw,yb
      integer i,k,m,ll,la,jnt,niter,l1,jj,j,kk,ic,kp3,
     1l,nt,jp,mkp1,mkp2,mkp3,mkp4,iter,mm1,mk,lt,jl
      double precision wm,thr,bz,pp,ww
      double precision phi
      data niter,wm,thr /25,0.0001,0.0001/
      k=0
      do 2 i=1,n
      k=0
      do 1 m=1,nk
      if(tb(1,m).eq.0.0) go to 1
      k=k+1
      sc(i,k)=phi(m,i,n,x,tb,cm)
    1 continue
    2 continue
      if(k .ne. 0) go to 3
      az=dlog(az/(1.0-az))
      return
    3 mk=k
      a=az
      jnt=1
      go to 19
      entry logitc (n,x,y,w,nk,il,cm,tb,kp,kv,lp,lv,bz,tc,sc,ss,d)
      ll=1
      la=ll
      l1=la
      lt=0
    4 if(kp(1,ll).lt.0) go to 17
      do 5 i=1,n
      ss(i)=1.0
    5 continue
      if(kp(1,ll) .le. 0) go to 9
      jl=kp(1,ll)
      do 8 il=1,jl
      k=kp(2,ll)+il-1
      jj=kv(1,k)
      j=iabs(jj)
      kk=kv(2,k)
      do 7 i=1,n
      if(ss(i).eq.0.0) go to 7
      ic=x(i,j)+.1
      ss(i)=cm(ic+kk)
      if(jj .ge. 0) go to 7
      if(ss(i) .ne. 0.0) go to 6
      ss(i)=1.0
      go to 7
    6 ss(i)=0.0
    7 continue
    8 continue
      go to 10
    9 if(kp(3,ll) .gt. 0) go to 10
      ll=ll+1
      go to 4
   10 if(kp(3,ll) .gt. 0) go to 12
      lt=lt+1
      do 11 i=1,n
      sc(i,lt)=ss(i)
   11 continue
      go to 16
   12 kp3=kp(3,ll)
      do 15 m=1,kp3
      l=lp(1,l1)
      nt=lp(3,l1)
      do 14 jp=1,nt
      lt=lt+1
      do 13 i=1,n
      sc(i,lt)=ss(i)
   13 continue
      call que(jp,l,nt,lv(lp(2,l1)),n,x,tc(la),sc(1,lt))
   14 continue
      l1=l1+1
      la=la+nt*(5*l+1)
   15 continue
   16 ll=ll+1
      go to 4
   17 if(lt .ne. 0) go to 18
      bz=dlog(bz/(1.0-bz))
      return
   18 mk=lt
      a=bz
      jnt=2
   19 mkp1=mk+1
      mkp2=mk+2
      mkp3=mk+3
      mkp4=mk+4
      iter=0
      if(jnt .ne. 1) go to 21
      k=0
      do 20 m=1,nk
      if(tb(1,m).eq.0.0) go to 20
      k=k+1
      d(k,mkp3)=tb(1,m)
   20 continue
      go to 27
   21 ll=1
      l1=ll
      la=0
      lt=la
   22 if(kp(1,ll).lt.0) go to 27
      if(kp(1,ll) .ne. 0 .or. kp(3,ll) .gt. 0) go to 23
      ll=ll+1
      go to 22
   23 if(kp(3,ll) .gt. 0) go to 24
      lt=lt+1
      d(lt,mkp3)=tc(-kp(3,ll))
      ll=ll+1
      go to 22
   24 kp3=kp(3,ll)
      do 26 m=1,kp3
      nt=lp(3,l1)
      la=la+5*lp(1,l1)*nt
      do 25 i=1,nt
      lt=lt+1
      d(lt,mkp3)=tc(i+la)
   25 continue
      la=la+nt
      l1=l1+1
   26 continue
      ll=ll+1
      go to 22
   27 iter=iter+1
      b=0.d0
      sw=b
      yb=sw
      do 29 i=1,n
      s=a
      do 28 m=1,mk
      s=s+d(m,mkp3)*sc(i,m)
   28 continue
      sc(i,mkp3)=s
      pp=1.0/(1.0+exp(-sc(i,mkp3)))
      ww=dmax1(pp*(1.0D0-pp),wm)    
      sc(i,mkp3)=sc(i,mkp3)+(y(i)-pp)/ww
      if(il.eq.2) ww=ww**2
      ww=ww*w(i)
      sc(i,mkp2)=ww
      sw=sw+ww
      yb=yb+ww*sc(i,mkp3)
      if(iter.gt.1) b=b+abs(pp-sc(i,mkp1))
      sc(i,mkp1)=pp
   29 continue
      if(iter.gt.niter.or.(iter.gt.1.and.b/n.lt.thr)) go to 37
      yb=yb/sw
      do 36 m=1,mk
      b=0.d0
      do 30 i=1,n
      b=b+sc(i,mkp2)*sc(i,m)
   30 continue
      b=b/sw
      mm1=m-1
      l=1
      go to 32
   31 l=l+1
   32 if((l).gt.(mm1)) go to 34
      s=0.d0
      do 33 i=1,n
      s=s+sc(i,mkp2)*(sc(i,m)-b)*sc(i,l)
   33 continue
      d(l,m)=s
      go to 31
   34 a=0.d0
      s=a
      do 35 i=1,n
      ww=sc(i,mkp2)
      pp=sc(i,m)-b
      s=s+ww*pp**2
      a=a+ww*pp*sc(i,mkp3)
   35 continue
      d(m,m)=s
      d(m,mkp1)=a
      d(m,mkp2)=b
   36 continue
      call lsf(nk,mk,mkp1,yb,d,d(1,mkp3),a,s,d(1,mkp4),1)
      go to 27
   37 if(jnt .ne. 1) go to 39
      az=a
      k=0
      do 38 m=1,nk
      if(tb(1,m).eq.0.0) go to 38
      k=k+1
      tb(1,m)=d(k,mkp3)
   38 continue
      go to 45
   39 bz=a
      ll=1
      l1=ll
      la=0
      lt=la
   40 if(kp(1,ll).lt.0) go to 45
      if(kp(1,ll) .ne. 0 .or. kp(3,ll) .gt. 0) go to 41
      ll=ll+1
      go to 40
   41 if(kp(3,ll) .gt. 0) go to 42
      lt=lt+1
      tc(-kp(3,ll))=d(lt,mkp3)
      ll=ll+1
      go to 40
   42 kp3=kp(3,ll)
      do 44 m=1,kp3
      nt=lp(3,l1)
      la=la+5*lp(1,l1)*nt
      do 43 i=1,nt
      lt=lt+1
      tc(i+la)=d(lt,mkp3)
   43 continue
      la=la+nt
      l1=l1+1
   44 continue
      ll=ll+1
      go to 40
   45 return
      end
      subroutine varimp (n,p,x,y,w,nk,il,it,az,tb,cm,vip,sc,d)
      implicit none
      integer n,p,nk,il,it
      double precision az,x(n,p),y(n),w(n),tb(5,nk),cm(*),vip(p),sc(n,*)
      double precision d(nk,*),sw,yb,yv,wn
      integer i,ip,j,nd
      double precision a0,cst,g,g0
      sw=0.d0
      yb=sw
      yv=yb
      wn=yv
      do 1 i=1,n
      sw=sw+w(i)
      wn=wn+w(i)**2
      yb=yb+w(i)*y(i)
    1 continue
      yb=yb/sw
      do 2 i=1,n
      yv=yv+w(i)*(y(i)-yb)**2
    2 continue
      yv=yv/sw
      wn=sw**2/wn
      ip=nk+4
      call varz(0,nk,tb,sc(1,ip),cst,nd)
      if(cst .ne. 1.0) go to 3
      g0=0.0
      if(il.gt.0) g0=yv
      go to 4
    3 a0=az
      call vp(n,x,y,w,nk,il,yb,sw,a0,sc(1,ip),cm,g0,sc,d)
    4 cst=1.d0/(1.d0-cst/wn)**2
      if(il .ne. 0) go to 5
      g0=(g0+yv)*cst
      go to 6
    5 g0=g0*cst
    6 do 12 j=1,p
      call varz(j,nk,tb,sc(1,ip),cst,nd)
      if(nd .ne. 0) go to 7
      vip(j)=g0
      go to 12
    7 if(cst .ne. 1.0) go to 8
      g=0.0
      if(il.gt.0) g=yv
      go to 9
    8 a0=az
      call vp(n,x,y,w,nk,il,yb,sw,a0,sc(1,ip),cm,g,sc,d)
    9 cst=1.d0/(1.d0-cst/wn)**2
      if(il .ne. 0) go to 10
      g=(g+yv)*cst
      go to 11
   10 g=g*cst
   11 vip(j)=g
   12 continue
      if(it .le. 0) go to 13
      write(it,17)
      call numprt(it,p,vip)
   13 a0=0.0
      do 14 j=1,p
      vip(j)=dsqrt(dmax1(0.0D0,vip(j)-g0))
      a0=dmax1(a0,vip(j))
   14 continue
      if(a0.le.0.0) return
      do 15 j=1,p
      vip(j)=100.0*vip(j)/a0
   15 continue
      if(it .le. 0) go to 16
      write(it,18)
      call numprt(it,p,vip)
   16 return
   17 format(/,' -gcv removing each variable:')
   18 format(/,' relative variable importance:')
      end
      subroutine numprt(it,n,a)
      implicit none
      integer it,n,i
      double precision a(*)
      integer i1,i2
      i2=0
    1 if(i2.ge.n) go to 2
      i1=i2+1
      i2=i2+6
      if(i2.gt.n) i2=n
      write(it,'(/,'' '',6(''    '',i4,''    ''))') (i,i=i1,i2)
      write(it,'('' '',6g12.4)') (a(i),i=i1,i2)
      go to 1
    2 return
      end
      subroutine varz(j,nk,tb,ub,cst,nd)
      implicit none
      integer j,nk,nd
      double precision tb(5,nk),ub(5,nk)
      integer m,k
      double precision cst
      integer jf
      do 2 m=1,nk
      do 1 k=1,5
      ub(k,m)=tb(k,m)
    1 continue
    2 continue
      nd=0
      if(j .le. 0) go to 4
      do 3 m=1,nk
      if(ub(1,m).eq.0.0) go to 3
      if(jf(m,j,ub) .eq. 0) go to 3
      ub(1,m)=0.0
      nd=nd+1
    3 continue
    4 cst=1.0
      do 5 m=1,nk
      if(ub(1,m).ne.0.0) cst=cst+ub(5,m)
    5 continue
      return
      end
      subroutine vp (n,x,y,w,nk,il,yb,sw,az,tb,cm,gof,sc,d)
      implicit none
      integer n,nk,il
      double precision az,gof,x(n,*),y(n),w(n),tb(5,nk),cm(*),sc(n,nk)
      double precision d(nk,*),s,t,yb,sw
      integer i,m,k
      double precision a,pp
      if(il .ne. 0) go to 1
      call lstsqr(n,x,y,w,nk,yb,sw,tb,cm,gof,sc,d)
      return
    1 call logitl(n,x,y,w,nk,il,az,tb,cm,sc,d)
      t=0.d0
      do 3 i=1,n
      s=az
      k=0
      do 2 m=1,nk
      if(tb(1,m).eq.0.0) go to 2
      k=k+1
      s=s+tb(1,m)*sc(i,k)
    2 continue
      a=s
      pp=1.0/(1.0+exp(-a))
      t=t+w(i)*(y(i)-pp)**2
    3 continue
      gof=t/sw
      return
      end
      subroutine lstsqr (n,x,y,w,nk,yb,sw,tb,cm,gof,sc,d)
      implicit none
      integer n,nk
      double precision gof,x(n,*),y(n),w(n),tb(5,nk),cm(*),sc(n,nk)
      double precision d(nk,*),a,b,s,yb,sw
      integer i,k,l,m,mk,mkp1,mkp2,mkp3,mkp4,mm1
      double precision pp,ww
      double precision phi
      k=0
      do 2 i=1,n
      k=0
      do 1 m=1,nk
      if(tb(1,m).eq.0.0) go to 1
      k=k+1
      sc(i,k)=phi(m,i,n,x,tb,cm)
    1 continue
    2 continue
      mk=k
      mkp1=mk+1
      mkp2=mk+2
      mkp3=mk+3
      mkp4=mk+4
      do 9 m=1,mk
      b=0.d0
      do 3 i=1,n
      b=b+w(i)*sc(i,m)
    3 continue
      b=b/sw
      mm1=m-1
      l=1
      go to 5
    4 l=l+1
    5 if((l).gt.(mm1)) go to 7
      s=0.d0
      do 6 i=1,n
      s=s+w(i)*(sc(i,m)-b)*sc(i,l)
    6 continue
      d(l,m)=s
      go to 4
    7 a=0.d0
      s=a
      do 8 i=1,n
      ww=w(i)
      pp=sc(i,m)-b
      s=s+ww*pp**2
      a=a+ww*pp*y(i)
    8 continue
      d(m,m)=s
      d(m,mkp1)=a
      d(m,mkp2)=b
    9 continue
      call lsf(nk,mk,mkp1,yb,d,d(1,mkp3),a,s,d(1,mkp4),1)
      gof=s/sw
      return
      end
      function efp(l,jv,nk,tb)
      implicit none
      double precision efp
      integer l,nk,jv(l)
      double precision tb(5,nk)
      integer j,k,m
      integer nord,jf
      efp=0.0
      do 3 m=1,nk
      if(tb(1,m).eq.0.0) go to 3
      if(nord(m,tb).ne.l) go to 3
      k=0
      do 1 j=1,l
      if(jf(m,jv(j),tb).eq.1) go to 1
      k=1
      go to 2
    1 continue
    2 if(k.eq.1) go to 3
      efp=efp+tb(5,m)
    3 continue
      return
      end
      function elg(jv,l,lx,tb,cm)
      implicit none
      logical elg
      integer jv,l,lx(*)
      double precision tb(5,*),cm(*)
      integer ic,kx,ip,jl,k,i1,jb
      integer nnord,nordc
      logical stelg
      data ic /0/
      elg=.false.
      kx=iabs(lx(jv))
      if(kx.eq.0) return
      if(l .ne. 0) go to 1
      elg=.true.
      return
    1 if((kx .ne. 2) .and. (kx .ne. 3)) go to 2
      if(nnord(l,tb).gt.0) return
    2 ip=l
    3 if(ip.le.0) go to 4
      jl=abs(tb(2,ip))+.1
      ip=tb(4,ip)+.1
      go to 3
    4 k=iabs(lx(jl))
      call isnstr(jl,jb)
      if((k.eq.2.or.k.eq.3).and.jb.eq.0) return
      if(ic .ne. 1) go to 5
      if(lx(jv).lt.0.and.nordc(1,l,tb,cm).gt.0) return
      if(lx(jv).gt.0.and.nordc(2,l,tb,cm).gt.0) return
      go to 6
    5 if(ic .ne. 2) go to 6
      if(lx(jv).gt.0.and.nordc(1,l,tb,cm).ge.2) return
    6 ip=l
    7 if(ip.le.0) go to 8
      jl=abs(tb(2,ip))+.1
      call intalw(jv,jl,k)
      if(k.eq.0) return
      ip=tb(4,ip)+.1
      go to 7
    8 elg=.true.
      return
      entry stelg(i1)
      ic=i1
      return
      end
      function phi(m,i,n,x,tb,cm)
      implicit none
      double precision phi
      integer m,i,n
      double precision tb(5,*),x(n,*),cm(*)
      integer ip,j
      double precision t,u
      phi=1.0
      ip=m
    1 if(ip.le.0) go to 7
      t=tb(2,ip)
      j=abs(t)+.1
      if(cm(2*j) .le. 0.0) go to 4
      u=cm(int(x(i,j)+.1)+int(tb(3,ip)+.1))
      if(t .ge. 0.0) go to 5
      if(u .ne. 0.0) go to 2
      u=1.0
      go to 5
    2 u=0.0
      go to 5
    4 u=dmax1(0.0D0,sign(1.0D0,t)*(x(i,j)-tb(3,ip)))
    5 if(u .gt. 0.0) go to 6
      phi=0.0
      return
    6 phi=phi*u
      ip=tb(4,ip)+.1
      go to 1
    7 return
      end
      function nord(m,tb)
      implicit none
      integer nord
      integer m
      double precision tb(5,*)
      integer ip
      ip=m
      nord=0
    1 if(ip.le.0) go to 2
      nord=nord+1
      ip=tb(4,ip)+.1
      go to 1
    2 return
      end
      function jf(m,j,tb)
      implicit none
      integer jf
      integer m,j
      double precision tb(5,*)
      integer ip,jp
      ip=m
      jf=0
    1 if(ip.le.0) go to 2
      jp=abs(tb(2,ip))+.1
      if(jp.eq.j) jf=1
      ip=tb(4,ip)+.1
      go to 1
    2 return
      end
      subroutine lsf(nk,m,mkp1,yb,d,a,a0,gf,dp,k1)
      implicit none
      integer nk,m,mkp1,k1
      double precision a(m),d(nk,*),dp(nk,*),eps,yb,a0,gf,s,t,sl,dps
      integer mkp2,k,i,info
      double precision big
      data big,eps /9.9e30,1.d-05/
      mkp2=mkp1+1
      gf=big
      if(d(m,m).le.0.d0) return
      sl=1.d0+eps
      do 2 k=k1,m
      do 1 i=1,k
      dp(i,k)=d(i,k)
    1 continue
      dp(k,k)=dp(k,k)*sl
    2 continue
      do 3 k=1,m
      a(k)=d(k,mkp1)
    3 continue
      info=k1
      call spofa(dp,nk,m,info)
      if(info.ne.0) return
      call sposl(dp,nk,m,a)
      s=yb
      t=0.d0
      do 4 i=1,m
      s=s-a(i)*d(i,mkp2)
      t=t-a(i)*(d(i,mkp1)+eps*d(i,i)*a(i))
    4 continue
      a0=s
      gf=t
      return
      entry seteps(dps)
      eps=dps
      return
      end
      subroutine coll(nk,tb,lp,lv,jv)
      implicit none
      integer nk,lp(3,*),lv(*),jv(*)
      integer nord
      double precision tb(5,nk)
      integer m,mo,l1,l2,mt,l10,jg,l1m1,i,j,k,
     1ig
      mo=0
      do 1 m=1,nk
      if(tb(1,m).ne.0.0) mo=max0(mo,nord(m,tb))
    1 continue
      if(mo .ne. 0) go to 2
      lp(1,1)=0
      return
    2 l1=1
      l2=l1
      do 11 mt=1,mo
      l10=l1
      do 10 m=1,nk
      if(tb(1,m).eq.0.0.or.nord(m,tb).ne.mt) go to 10
      call jfv(m,tb,jv)
      jg=0
      l1m1=l1-1
      i=l10
      go to 4
    3 i=i+1
    4 if((i).gt.(l1m1)) go to 8
      k=lp(2,i)-1
      ig=0
      do 5 j=1,mt
      if(jv(j).eq.lv(k+j)) go to 5
      ig=1
      go to 6
    5 continue
    6 if(ig .ne. 0) go to 3
      jg=1
      lp(3,i)=lp(3,i)+1
    8 if(jg .ne. 0) go to 10
      lp(1,l1)=mt
      lp(2,l1)=l2
      lp(3,l1)=1
      k=l2-1
      do 9 i=1,mt
      lv(i+k)=jv(i)
    9 continue
      l1=l1+1
      l2=l2+mt
   10 continue
   11 continue
      lp(1,l1)=0
      return
      end
      subroutine jfv(m,tb,jv)
      implicit none
      integer m,jv(*)
      double precision tb(5,*)
      integer ip,j,l,i,k
      ip=m
      j=0
    1 if(ip.le.0) go to 2
      j=j+1
      jv(j)=abs(tb(2,ip))+.1
      ip=tb(4,ip)+.1
      go to 1
    2 if(j.eq.1) return
      j=j-1
    3 l=0
      do 4 i=1,j
      if(jv(i) .le. jv(i+1)) go to 4
      k=jv(i)
      jv(i)=jv(i+1)
      jv(i+1)=k
      l=1
    4 continue
      if(l.eq.0) go to 5
      go to 3
    5 return
      end
      subroutine side(l,nt,jv,xe,x)
      implicit none
      integer l,nt,jv(l)
      double precision xe(2,*),x(nt,*)
      integer l2,l3,l4,k,j,m
      double precision xl,xr,z,dl,x1,x2,a,dx,dr
      l2=l+l
      l3=l2+l
      l4=l3+l
      do 7 k=1,l
      xl=xe(1,jv(k))
      xr=xe(2,jv(k))
      do 6 j=1,nt
      z=x(j,k)
      if(z .gt. xl) go to 1
      x(j,k+l)=xl
      x(j,k+l2)=x(j,k+l)
      x(j,k+l3)=0.0
      x(j,k+l4)=x(j,k+l3)
      go to 6
    1 dl=z-xl
      dr=xr-z
      x1=xl
      x2=xr
      do 3 m=1,nt
      a=x(m,k)
      if(a.eq.z) go to 3
      dx=a-z
      if(dx .ge. 0.0 .or. -dx .ge. dl) go to 2
      dl=-dx
      x1=a
    2 if(dx .le. 0.0 .or. dx .ge. dr) go to 3
      dr=dx
      x2=a
    3 continue
      x1=0.5*(x1+z)
      x2=0.5*(x2+z)
      if(x(j,k+l) .le. 0.0) go to 4
      x(j,k+l)=x1
      x(j,k+l2)=x2
      go to 5
    4 x(j,k+l)=x2
      x(j,k+l2)=x1
    5 call pr(x(j,k+l),x(j,k),x(j,k+l2),x(j,k+l3),x(j,k+l4))
    6 continue
    7 continue
      return
      end
      subroutine que(jp,l,nt,jv,n,x,tc,t)
      implicit none
      integer jp,l,nt,n,jv(l)
      double precision x(n,*),tc(nt,*),t(n)
      integer l2,l3,l4,i,j,k
      double precision q
      double precision cue
      l2=l+l
      l3=l2+l
      l4=l3+l
      do 3 i=1,n
      if(t(i).eq.0.0) go to 3
      q=1.0
      do 1 k=1,l
      j=jv(k)
      q=q*cue(x(i,j),tc(jp,k+l),tc(jp,k),tc(jp,k+l2),tc(jp,k+l3),tc(jp,k
     1+l4))
      if(q.eq.0.0) go to 2
    1 continue
    2 t(i)=q
    3 continue
      return
      end
      subroutine scpc(xm,xs,jp,l,nt,jv,tc,b)
      implicit none
      integer jp,l,nt,jv(l)
      double precision b,xm(*),xs(*),tc(nt,*)
      double precision g,h,q
      integer l2,l3,l4,j,k
      l2=l+l
      l3=l2+l
      l4=l3+l
      q=1.d0
      do 1 k=1,l
      j=jv(k)
      g=xm(j)
      h=xs(j)
      q=q*h
      tc(jp,k+l)=g+h*tc(jp,k+l)
      tc(jp,k)=g+h*tc(jp,k)
      tc(jp,k+l2)=g+h*tc(jp,k+l2)
      tc(jp,k+l3)=tc(jp,k+l3)/h
      tc(jp,k+l4)=tc(jp,k+l4)/h**2
    1 continue
      b=b/q
      return
      end
      subroutine update(il,n,m,kr,x,y,w,sw,yb,tb,cm,sc,bl,d,dy,db)
      implicit none
      integer il,n,m,kr
      double precision x(n,*),y(n),w(n),tb(5,*),cm(*),sc(n,*),bl(n)
      double precision d(n,*),dy(*),db(*),b,s,yb,sw,dv,eps,v,q
      integer i,j,k,kp,nw,n0
      double precision h,t,tk,u,sg
      data eps /1.d-4/
      kp=kr+1
      b=0.d0
      t=tb(2,m)
      j=abs(t)+.1
      if(il .ne. 1) go to 3
      tk=tb(3,m)
      do 2 i=1,n
      h=bl(i)
      if(h .gt. 0.0) go to 1
      sc(i,m)=0.0
      go to 2
    1 sc(i,m)=h*(x(i,j)-tk)
      b=b+w(i)*sc(i,m)
    2 continue
      go to 17
    3 if(cm(2*j) .le. 0.0) go to 12
      k=tb(3,m)+.1
      nw=0
      n0=nw
      do 11 i=1,n
      h=bl(i)
      if(h .gt. 0.0) go to 4
      sc(i,m)=0.0
      go to 11
    4 u=cm(int(x(i,j)+.1)+k)
      if(w(i) .le. 0.0) go to 5
      nw=nw+1
      if(u.eq.0.0) n0=n0+1
    5 if(t .ge. 0.0) go to 8
      if(u .ne. 0.0) go to 6
      sc(i,m)=h
      go to 10
    6 sc(i,m)=0.0
      go to 10
    8 if(u .gt. 0.0) go to 9
      sc(i,m)=0.0
      go to 10
    9 sc(i,m)=h
   10 b=b+w(i)*sc(i,m)
   11 continue
      if(n0.eq.0.or.n0.eq.nw) return
      go to 17
   12 tk=tb(3,m)
      sg=sign(1.0D0,t)
      do 16 i=1,n
      h=bl(i)
      if(h .gt. 0.0) go to 13
      sc(i,m)=0.0
      go to 16
   13 u=dmax1(0.0D0,sg*(x(i,j)-tk))
      if(u .gt. 0.0) go to 14
      sc(i,m)=0.0
      go to 15
   14 sc(i,m)=h*u
   15 b=b+w(i)*sc(i,m)
   16 continue
   17 b=b/sw
      s=0.d0
      v=s
      do 18 j=1,kr
      db(j)=0.d0
   18 continue
      do 21 i=1,n
      d(i,kp)=sc(i,m)-b
      if(sc(i,m).le.0.0) go to 21
      q=w(i)*sc(i,m)
      s=s+q*(sc(i,m)-b)
      v=v+q*(y(i)-yb)
      j=1
      go to 20
   19 j=j+1
   20 if((j).gt.(kr)) go to 21
      db(j)=db(j)+q*d(i,j)
      go to 19
   21 continue
      if(s.le.0.d0) return
      dv=s
      j=1
      go to 23
   22 j=j+1
   23 if((j).gt.(kr)) go to 24
      s=s-db(j)**2
      go to 22
   24 if(s.lt.eps*dv) return
      j=1
      go to 26
   25 j=j+1
   26 if((j).gt.(kr)) go to 28
      do 27 i=1,n
      d(i,kp)=d(i,kp)-db(j)*d(i,j)
   27 continue
      go to 25
   28 s=1.d0/dsqrt(s)
      j=1
      go to 30
   29 j=j+1
   30 if((j).gt.(kr)) go to 31
      v=v-db(j)*dy(j)
      go to 29
   31 dy(kp)=v*s
      do 32 i=1,n
      d(i,kp)=d(i,kp)*s
   32 continue
      kr=kp
      return
      end
      subroutine pr(um,u,up,p,r)
      implicit none
      double precision um,u,up,p,r
      double precision s
      s=1.0
      if(um.gt.up) s=-1.0
      p=s*(2.0*up+um-3.0*u)/(up-um)**2
      r=s*(2.0*u-up-um)/(up-um)**3
      return
      end
      function cue(x,um,u,up,p,r)
      implicit none
      double precision cue
      double precision x,um,u,up,p,r
      double precision s,y
      s=1.0
      if(um.gt.up) s=-1.0
      y=s*x
      if(y .gt. s*um) go to 1
      cue=0.0
      return
    1 if(y .lt. s*up) go to 2
      cue=y-s*u
      return
    2 cue=p*(x-um)**2+r*(x-um)**3
      return
      end
      function varf(nk,d,a,sw,k1,k2)
      implicit none
      double precision varf
      integer nk,k1,k2
      double precision d(nk,*),a(nk),sw,s,t,u
      integer i,j
      s=0.d0
      do 4 i=k1,k2
      t=0.d0
      do 3 j=k1,k2
      if(j .gt. i) go to 1
      u=d(j,i)
      go to 2
    1 u=d(i,j)
    2 t=t+a(j)*u
    3 continue
      s=s+a(i)*t
    4 continue
      varf=s/sw
      return
      end
      subroutine exch(nk,m,k,d,a,b)
      implicit none
      integer nk,m,k
      double precision a(m),b(m)
      double precision d(nk,m),t
      integer i,j,l,km1,kp2
      double precision r
      l=k+1
      km1=k-1
      kp2=k+2
      r=a(k)
      a(k)=a(l)
      a(l)=r
      r=b(k)
      b(k)=b(l)
      b(l)=r
      do 1 j=1,2
      i=nk+j
      t=d(k,i)
      d(k,i)=d(l,i)
      d(l,i)=t
    1 continue
      t=d(k,k)
      d(k,k)=d(l,l)
      d(l,l)=t
      j=1
      go to 3
    2 j=j+1
    3 if((j).gt.(km1)) go to 4
      t=d(j,k)
      d(j,k)=d(j,l)
      d(j,l)=t
      go to 2
    4 j=kp2
      go to 6
    5 j=j+1
    6 if((j).gt.(m)) go to 7
      t=d(k,j)
      d(k,j)=d(l,j)
      d(l,j)=t
      go to 5
    7 return
      end
      function jft(m,j,tb)
      implicit none
      integer jft
      integer m,j
      double precision tb(5,*)
      integer k
      double precision abs
      k=1
      go to 2
    1 k=k+1
    2 if((k).gt.(m)) go to 4
      if(int(abs(tb(2,k))+.1) .ne. j) go to 1
      jft=1
      return
    4 jft=0
      return
      end
      subroutine coefpr (it,nk,az,tb,cm,xs)
      implicit none
      integer it,nk,i
      double precision az,tb(5,*),cm(*),xs(*),a(6)
      integer i2,i1,l2
      integer min0
      i2=0
    1 if(i2.ge.nk) go to 4
      if(i2 .ne. 0) go to 2
      i1=0
      i2=min0(5,nk)
      l2=i2+1
      a(1)=az
      call org(1,i2,tb,cm,xs,a(2))
      go to 3
    2 i1=i2+1
      i2=i2+6
      if(i2.gt.nk) i2=nk
      l2=i2-i1+1
      call org(i1,i2,tb,cm,xs,a)
    3 continue
      write(it,'(/,'' bsfn:'',6(''    '',i4,''    ''))') (i,i=i1,i2)
      write(it,'('' coef:'',6g12.4)') (a(i),i=1,l2)
      go to 1
    4 return
      end
      subroutine org(m1,m2,tb,cm,xs,a)
      implicit none
      integer m1,m2
      double precision xs(*),tb(5,*),cm(*),a(*)
      integer ip,j,k,m
      double precision s
      k=0
      do 4 m=m1,m2
      k=k+1
      if(tb(1,m) .ne. 0.0) go to 1
      a(k)=0.0
      go to 4
    1 s=1.0
      ip=m
    2 if(ip.le.0) go to 3
      j=abs(tb(2,ip))+.1
      if(cm(2*j).eq.0.0) s=s*xs(j)
      ip=tb(4,ip)+.1
      go to 2
    3 a(k)=tb(1,m)/s
    4 continue
      return
      end
      subroutine hulset (n,x,big,nh,xh,y)
      implicit none
      integer n,nh
      double precision big,x(n,*),y(n),xh(3,nh)
      integer i,j,k
      double precision a,b,s,sg,x1,x2
      do 5 j=1,n
      k=0
      x1=x(j,1)
      x2=x(j,2)
      do 3 i=1,nh
      a=xh(1,i)
      b=xh(2,i)
      sg=xh(3,i)
      if(a .lt. big) go to 1
      s=x1-b
      go to 2
    1 s=x2-a*x1-b
    2 if(s*sg .ge. 0.0) go to 3
      k=1
      go to 4
    3 continue
    4 if(k.eq.1) y(j)=big
    5 continue
      return
      end
      subroutine cvxhul (n,x1,x2,big,nh,xh)
      implicit none
      integer n,nh
      double precision big,x1(n),x2(n),xh(3,*)
      integer i,iq,k,kq,k0,lq
      double precision a,am,a0,eps,x0,y0,xm,ym,xb,yb,
     1x,y,b,s
      data eps /1.e-3/
      k0=0
      x0=big
      y0=x0
      xm=-big
      ym=xm
      xb=0.0
      yb=xb
      do 1 i=1,n
      xb=xb+x1(i)
      yb=yb+x2(i)
      x0=dmin1(x0,x1(i))
      y0=dmin1(y0,x2(i))
      xm=dmax1(xm,x1(i))
      ym=dmax1(ym,x2(i))
    1 continue
      x0=x0-eps*(xm-x0)
      y0=y0-eps*(ym-y0)
      xb=xb/n
      yb=yb/n
      nh=0
      a0=0.0
      lq=1
    2 am=big
      kq=4
      k=0
      do 15 i=1,n
      x=x1(i)
      y=x2(i)
      if(x .ne. x0) go to 5
      if(y.eq.y0) go to 15
      if(y .le. y0) go to 3
      iq=2
      go to 10
    3 iq=4
      go to 10
    5 if(x .le. x0) go to 8
      if(y .lt. y0) go to 6
      iq=1
      go to 10
    6 iq=4
      go to 10
    8 if(y .le. y0) go to 9
      iq=2
      go to 10
    9 iq=3
   10 if(iq.gt.kq) go to 15
      if(iq.lt.lq) go to 15
      if((iq .ne. 1) .and. (iq .ne. 3)) go to 11
      a=abs((y-y0)/(x-x0))
      go to 12
   11 a=abs((x-x0)/(y-y0))
   12 if(iq.eq.lq.and.a.lt.a0) go to 15
      if(iq .ge. kq) go to 13
      kq=iq
      am=a
      k=i
      go to 15
   13 if(a .ge. am) go to 14
      am=a
      k=i
      go to 15
   14 if(a .ne. am) go to 15
      if((x-x0)**2+(y-y0)**2.gt.(x1(k)-x0)**2+(x2(k)-y0)**2) k=i
   15 continue
      if(k .ne. 0) go to 16
      a0=0.0
      lq=1
      go to 2
   16 if(nh .ne. 0) go to 17
      k0=k
      go to 20
   17 if(x1(k) .ne. x0) go to 18
      a=big
      b=x0
      s=xb-b
      go to 19
   18 a=(x2(k)-y0)/(x1(k)-x0)
      b=y0-a*x0
      s=yb-a*xb-b
   19 xh(1,nh)=a
      xh(2,nh)=b
      xh(3,nh)=sign(1.0D0,s)
      if(k.eq.k0) go to 21
   20 nh=nh+1
      x0=x1(k)
      y0=x2(k)
      lq=kq
      a0=am
      go to 2
   21 return
      end
      subroutine knts (l,nt,jv,jl,kv,nk,tb,cm,x,js)
      implicit none
      integer l,nt,jl,nk,jv(l),kv(2,jl),js(*)
      double precision tb(5,nk),cm(*),x(nt,*)
      integer l1,m,k,j,ip
      double precision t
      integer jf,icf,nordc
      l1=0
      do 7 m=1,nk
      if(icf(m,tb,cm,jl,kv,js).eq.0) go to 7
      if(nordc(1,m,tb,cm).ne.l) go to 7
      k=0
      do 1 j=1,l
      if(jf(m,jv(j),tb).eq.1) go to 1
      k=1
      go to 2
    1 continue
    2 if(k.eq.1) go to 7
      ip=m
      l1=l1+1
    3 if(ip.le.0) go to 7
      t=tb(2,ip)
      j=abs(t)+.1
      if(cm(2*j) .eq. 0.0) go to 4
      ip=tb(4,ip)+.1
      go to 3
    4 k=1
    5 if(jv(k).eq.j) go to 6
      k=k+1
      go to 5
    6 x(l1,k)=tb(3,ip)
      x(l1,l+k)=sign(1.0D0,t)
      ip=tb(4,ip)+.1
      go to 3
    7 continue
      return
      end
      subroutine sclato (n,p,x,xm,xs,cm,z)
      implicit none
      integer n,p
      double precision x(n,p),xm(p),xs(p),cm(*),z(n,p)
      integer j,j1,i,l
      do 4 j=1,p
      j1=cm(2*j)+.1
      if(j1 .ne. 0) go to 2
      if(xs(j).le.0.0) go to 4
      do 1 i=1,n
      z(i,j)=xs(j)*x(i,j)+xm(j)
    1 continue
      go to 4
    2 j1=j1-1
      do 3 i=1,n
      l=x(i,j)+.1
      z(i,j)=cm(l+j1)
    3 continue
    4 continue
      call stfmrs(0)
      call stcmrs(0)
      return
      end
      function ncat(kp)
      implicit none
      integer ncat
      integer kp(5,*)
      integer ll
      ncat=0
      ll=1
    1 if(kp(1,ll).lt.0) go to 2
      if(kp(1,ll).gt.0.and.kp(3,ll).le.0) ncat=ncat+1
      ll=ll+1
      go to 1
    2 return
      end
      subroutine catv (jl,kp,kv,nv,jv)
      implicit none
      integer jl,nv,kp(5,*),kv(2,*),jv(jl,*)
      integer ll,jg,j,ig,k1,l1,i,l,k
      nv=0
      ll=1
    1 if(kp(1,ll).lt.0) go to 20
      if(kp(3,ll) .le. 0) go to 2
      ll=ll+1
      go to 1
    2 if(kp(1,ll) .eq. jl) go to 3
      ll=ll+1
      go to 1
    3 jg=0
      j=1
      go to 5
    4 j=j+1
    5 if((j).gt.(nv)) go to 9
      ig=0
      do 6 i=1,jl
      if(jv(i,j).eq.iabs(kv(1,kp(2,ll)+i-1))) go to 6
      ig=1
      go to 7
    6 continue
    7 if(ig .ne. 0) go to 4
      jg=1
    9 if(jg .ne. 1) go to 10
      ll=ll+1
      go to 1
   10 l1=ll+1
      ig=0
   11 if(kp(1,l1).lt.0) go to 17
      k1=kp(1,l1)
      if((k1 .ne. jl) .and. (kp(3,l1) .le. 0)) go to 12
      l1=l1+1
      go to 11
   12 do 15 i=1,k1
      k=iabs(kv(1,kp(2,l1)+i-1))
      do 13 j=1,jl
      l=iabs(kv(1,kp(2,ll)+j-1))
      if(l .ne. k) go to 13
      ig=1
      go to 14
   13 continue
   14 if(ig.eq.1) go to 16
   15 continue
   16 if(ig.eq.1) go to 17
      l1=l1+1
      go to 11
   17 if(ig .ne. 1) go to 18
      ll=ll+1
      go to 1
   18 nv=nv+1
      do 19 i=1,jl
      jv(i,nv)=iabs(kv(1,kp(2,ll)+i-1))
   19 continue
      ll=ll+1
      go to 1
   20 return
      end
      function cvlv (m,jl,jv,lv,nk,kp,kv,tb,cm,tc)
      implicit none
      double precision cvlv
      integer m,jl,nk,jv(jl),lv(jl),kp(5,*),kv(2,*)
      double precision tb(5,nk),cm(*),tc(*)
      double precision cvll,cvlq
      if(m .ne. 1) go to 1
      cvlv=cvll(jl,jv,lv,nk,tb,cm)
      go to 2
    1 cvlv=cvlq(jl,jv,lv,kp,kv,cm,tc)
    2 return
      end
      function cvll (jl,jv,lv,nk,tb,cm)
      implicit none
      double precision cvll
      integer jl,nk,jv(jl),lv(jl),iv(2)
      double precision tb(5,nk),cm(*)
      integer m,ig,i,ip,j
      double precision t,phi,u
      integer nord,nordc
      cvll=0.0
      if(jl.gt.2) return
      do 11 m=1,nk
      if(tb(1,m).eq.0.0) go to 11
      if(nordc(1,m,tb,cm).gt.0) go to 11
      if(nord(m,tb).ne.jl) go to 11
      call jfv(m,tb,iv)
      ig=0
      do 1 i=1,jl
      if(jv(i).eq.iv(i)) go to 1
      ig=1
      go to 2
    1 continue
    2 if(ig.eq.1) go to 11
      phi=1.0
      ip=m
    3 if(ip.le.0) go to 10
      t=tb(2,ip)
      j=abs(t)+.1
      i=1
      go to 5
    4 i=i+1
    5 if((i).gt.(jl)) go to 6
      if(jv(i).eq.j) go to 6
      go to 4
    6 u=cm(lv(i)+int(tb(3,ip)+.1))
      if(t .ge. 0.0) go to 8
      if(u .ne. 0.0) go to 7
      u=1.0
      go to 8
    7 u=0.0
    8 if(u .ne. 0.0) go to 9
      phi=0.0
      go to 10
    9 ip=tb(4,ip)+.1
      go to 3
   10 if(phi.gt.0.0) cvll=cvll+tb(1,m)
   11 continue
      return
      end
      function cvlq (jl,jv,lv,kp,kv,cm,tc)
      implicit none
      double precision cvlq
      integer jl,jv(jl),lv(jl),kp(5,*),kv(2,*)
      double precision cm(*),tc(*)
      integer ll,ig,kt,i,j,jt,k
      ll=1
      cvlq=0.0
    1 if(kp(1,ll).lt.0) go to 12
      if(kp(3,ll) .le. 0) go to 2
      ll=ll+1
      go to 1
    2 if(kp(1,ll) .eq. jl) go to 3
      ll=ll+1
      go to 1
    3 ig=0
      do 4 i=1,jl
      if(jv(i).eq.iabs(kv(1,kp(2,ll)+i-1))) go to 4
      ig=1
      go to 5
    4 continue
    5 if(ig .ne. 1) go to 6
      ll=ll+1
      go to 1
    6 kt=1
      do 9 j=1,jl
      k=kp(2,ll)+j-1
      jt=cm(lv(j)+kv(2,k))+.1
      if(kv(1,k) .ge. 0) go to 8
      if(jt .ne. 0) go to 7
      jt=1
      go to 8
    7 jt=0
    8 if(jt .ne. 0) go to 9
      kt=0
      go to 10
    9 continue
   10 if(kt .ne. 1) go to 11
      cvlq=cvlq+tc(-kp(3,ll))
   11 ll=ll+1
      go to 1
   12 return
      end
      subroutine ccoll (nk,tb,cm,kp,kv,lp,lv,jv)
      implicit none
      integer nk,kp(5,*),kv(2,*),lp(3,*),lv(*),jv(*)
      double precision tb(5,*),cm(*)
      integer ll,l1,l2,li, ii
c      print *, 'before collc, kp=', (kp(1,ii), ii=1,nk)
      call collc(nk,tb,cm,kp,kv,jv)
c      print *, 'in ccoll  kv=', kv(1,1)
c      print *, 'before purcat, kp=', (kp(1,ii), ii=1,nk)
      call purcat(nk,tb,cm,kp,kv,li,jv)
c      print *, 'after purcat, kp=', (kp(1,ii), ii=1,nk)
      ll=li+1
      l1=1
      l2=l1
    1 if(kp(1,ll).lt.0) go to 2
      kp(4,ll)=l1
      call collf(nk,tb,cm,kp(1,ll),kv(1,kp(2,ll)),l1,l2,lp,lv,jv)
      kp(3,ll)=l1-kp(4,ll)
      ll=ll+1
      go to 1
    2 lp(1,l1)=0
      return
      end
      subroutine collc (nk,tb,cm,kp,kv,jv)
      implicit none
      integer nk,kp(5,*),kv(2,*),jv(*)
      double precision tb(5,nk),cm(*)
      integer l1,l2,mc,m,mt,jg,l1m1,i,k,ig,m1,
     1m2,jj,nc,kk,jk,l10,nv,j
      double precision z
      integer nordc
      kp(1,1)=0
      kp(2,1)=1
      l1=2
      l2=1
      mc=0
      do 1 m=1,nk
      if(tb(1,m).ne.0.0) mc=max0(mc,nordc(2,m,tb,cm))
    1 continue
c      print *, 'in collc  mc=', mc
      mt=1
      go to 3
    2 mt=mt+1
    3 continue
      if((mt).gt.(mc)) go to 18
      l10=l1
      do 17 m=1,nk
      if(tb(1,m).eq.0.0.or.nordc(2,m,tb,cm).ne.mt) go to 17
      call jfvc(2,m,tb,cm,nv,jv,jv(mt+1))
      jg=0
      l1m1=l1-1
      i=l10
c      print *, 'in collc  after jfvc, jv=', jv(1)
      go to 5
    4 i=i+1
    5 if((i).gt.(l1m1)) go to 15
      k=kp(2,i)-1
      ig=0
      do 6 j=1,mt
      if(iabs(jv(j)).eq.iabs(kv(1,k+j))) go to 6
      ig=1
      go to 7
    6 continue
    7 if(ig .ne. 0) go to 13
      do 12 j=1,mt
      m1=kv(2,k+j)
      m2=jv(mt+j)
      jj=iabs(jv(j))
      nc=int(cm(2*jj+1)+.1)-int(cm(2*jj)+.1)+1
      kk=jv(j)*kv(1,k+j)
      do 10 jk=1,nc
      z=cm(jk+m2)
      if(kk .ge. 0) go to 9
      if(z .ne. 0.0) go to 8
      z=1.0
      go to 9
    8 z=0.0
    9 if(cm(jk+m1).eq.z) go to 10
      ig=1
      go to 11
   10 continue
   11 if(ig.eq.1) go to 13
   12 continue
   13 if(ig .ne. 0) go to 4
      jg=1
   15 continue
c      print *, 'in collc   jg, mt=',jg, mt
      if(jg .ne. 0) go to 17
      kp(1,l1)=mt
      kp(2,l1)=l2
      k=l2-1
c      print *, 'in collc   k, jv[1]=', k, jv(1)
      do 16 i=1,mt
      kv(1,i+k)=jv(i)
      kv(2,i+k)=jv(i+mt)
   16 continue
      l1=l1+1
      l2=l2+mt
c      print *, 'in collc   kv=', kv(1,1)
   17 continue
      go to 2
   18 kp(1,l1)=-1
      return
      end
      function nordc (l,m,tb,cm)
      implicit none
      integer nordc
      integer l,m
      double precision tb(5,*),cm(*)
      integer ip,j
      ip=m
      nordc=0
    1 if(ip.le.0) go to 4
      j=abs(tb(2,ip))+.1
      if(l .ne. 1) go to 2
      if(cm(2*j).eq.0.0) nordc=nordc+1
      go to 3
    2 if(cm(2*j).gt.0.0) nordc=nordc+1
    3 ip=tb(4,ip)+.1
      go to 1
    4 return
      end
      subroutine jfvc (l,m,tb,cm,nv,jv,jp)
      implicit none
      integer l,m,nv,jv(*),jp(*)
      double precision tb(5,*),cm(*)
      integer ip,j,ll,i,k
      ip=m
      nv=0
    1 if(ip.le.0) go to 5
      j=abs(tb(2,ip))+.1
      if(l .ne. 1) go to 3
      if(cm(2*j) .le. 0.0) go to 4
      ip=tb(4,ip)+.1
      go to 1
    3 if(cm(2*j) .ne. 0.0) go to 4
      ip=tb(4,ip)+.1
      go to 1
    4 nv=nv+1
      jv(nv)=j
      if(l.ne.1.and.tb(2,ip).lt.0.0) jv(nv)=-j
      if(l.ne.1) jp(nv)=tb(3,ip)+.1
      ip=tb(4,ip)+.1
      go to 1
    5 if(nv.le.1) return
      j=nv-1
    6 ll=0
      do 7 i=1,j
      if(iabs(jv(i)) .le. iabs(jv(i+1))) go to 7
      ll=1
      k=jv(i)
      jv(i)=jv(i+1)
      jv(i+1)=k
      if(l .eq. 1) go to 7
      k=jp(i)
      jp(i)=jp(i+1)
      jp(i+1)=k
    7 continue
      if(ll.eq.0) go to 8
      go to 6
    8 return
      end
      subroutine purcat (nk,tb,cm,kp,kv,li,jv)
      implicit none
      integer nk,li,kp(5,*),kv(2,*),jv(*)
      double precision tb(5,nk),cm(*)
      integer lm,ll,jl,ifg,jfg,m,i,j
      integer icf,nord, icf0, ii

      lm=1
    1 if(kp(1,lm).lt.0) go to 2
      lm=lm+1
      go to 1
    2 ll=1
      li=0
c      print *, 'START purcat    kv=',kv(1,1)
c      print *, 'in purcat  kp=', (kp(1,ii), ii=1,nk)
    3 if(kp(1,ll).lt.0) go to 20
c      print *, 'in purcat kp[1,ll]=', kp(1,ll)
      jl=kp(1,ll)
      if(jl .gt. 0) go to 4
      ll=ll+1
      go to 3
    4 ifg=0
      jfg=ifg
c      print *, 'in purcat    kv=',kv(1,1)
      do 6 m=1,nk
c      print *, 'before icf     kv(1,1)=', kv(1,1)
      icf0 = icf(m,tb,cm,jl,kv(1,kp(2,ll)),jv)
c      print *, 'after icf     icf, kv(1,1)=', icf0, kv(1,1)
c      if(icf(m,tb,cm,jl,kv(1,kp(2,ll)),jv).eq.0) go to 6
      if (icf0 .eq. 0) go to 6
      if(nord(m,tb) .ne. jl) go to 5
      ifg=1
      go to 6
    5 jfg=1
    6 continue
      if(ifg .ne. 0) go to 9
      if(jfg .ne. 0) go to 8
      write(6,7)
    7 format (' bug in purcat - term not found.')
      stop
    8 ll=ll+1
      go to 3
    9 li=li+1
      j=lm
c      print *, 'purcat 2, lm=', lm
      go to 11
   10 j=j+(-1)
   11 if((-1)*((j)-(li)).gt.0) go to 13
c     this loop had a problem in the original code
c     specifically it access kp(5,6) which collides
c     with the address of kv.  the remedy is noted 
c     in the subroutine mars.
      do 12 i=1,5
      kp(i,j+1)=kp(i,j)
   12 continue
      go to 10
   13 lm=lm+1
      ll=ll+1
      do 14 i=1,5
      kp(i,li)=kp(i,ll)
   14 continue
      kp(3,li)=0
      kp(4,li)=1
      kp(5,li)=0
      if(jfg .ne. 1) go to 15
      ll=ll+1
      go to 3
   15 j=ll+1
      go to 17
   16 j=j+1
   17 if((j).gt.(lm)) go to 19
      do 18 i=1,5
      kp(i,j-1)=kp(i,j)
   18 continue
      go to 16
   19 lm=lm-1
c      print *, 'in purcat  kp=', (kp(1,ii), ii=1,nk)
c      print *, 'in purcat    kv=',kv(1,1)
      go to 3
   20 continue
c      print *, 'end purcat'
      return
      end
      subroutine collf (nk,tb,cm,jl,kv,l1,l2,lp,lv,jv)
      implicit none
      integer nk,jl,l1,l2,kv(2,*),lp(3,*),lv(*),jv(*)
      double precision tb(5,*),cm(*)
      integer mo,m,mt,l10,jg,l1m1,i,k,ig,j,nv
      integer icf,nordc
      mo=0
      do 1 m=1,nk
      if(icf(m,tb,cm,jl,kv,jv).ne.0) mo=max0(mo,nordc(1,m,tb,cm))
    1 continue
      if(mo.eq.0) return
      do 10 mt=1,mo
      l10=l1
      do 9 m=1,nk
      if(icf(m,tb,cm,jl,kv,jv).eq.0) go to 9
      if(nordc(1,m,tb,cm).ne.mt) go to 9
      call jfvc(1,m,tb,cm,nv,jv,jv)
      jg=0
      l1m1=l1-1
      i=l10
      go to 3
    2 i=i+1
    3 if((i).gt.(l1m1)) go to 7
      k=lp(2,i)-1
      ig=0
      do 4 j=1,mt
      if(jv(j).eq.lv(k+j)) go to 4
      ig=1
      go to 5
    4 continue
    5 if(ig .ne. 0) go to 2
      jg=1
      lp(3,i)=lp(3,i)+1
    7 if(jg .ne. 0) go to 9
      lp(1,l1)=mt
      lp(2,l1)=l2
      lp(3,l1)=1
      k=l2-1
      do 8 i=1,mt
      lv(i+k)=jv(i)
    8 continue
      l1=l1+1
      l2=l2+mt
    9 continue
   10 continue
      return
      end
      function icf (m,tb,cm,jl,kv,jv)
      implicit none
      integer icf
      integer m,jl,kv(2,jl),jv(*)
      double precision tb(5,*),cm(*)
      integer nv,j,l1,l2,k,kk,nc,i
      double precision z
      integer nordc

      icf=0
      if(tb(1,m).eq.0.0.or.nordc(2,m,tb,cm).ne.jl) return
      if(jl .ne. 0) go to 1
      icf=1
      return
    1 call jfvc(2,m,tb,cm,nv,jv,jv(jl+1))
c      print *, 'in icf jl, jv[1], kv=',jl, jv(1), kv(1,1)
      do 2 j=1,jl
      if(iabs(jv(j)).ne.iabs(kv(1,j))) return
    2 continue
c      print *, 'in icf jl=', jl
      do 6 j=1,jl
      l1=kv(2,j)
      l2=jv(jl+j)
      k=2*iabs(jv(j))
      kk=jv(j)*kv(1,j)
      nc=int(cm(k+1)+.1)-int(cm(k)+.1)+1
      do 5 i=1,nc
      z=cm(i+l2)
      if(kk .ge. 0) go to 4
      if(z .ne. 0.0) go to 3
      z=1.0
      go to 4
    3 z=0.0
    4 if(cm(i+l1).ne.z) return
    5 continue
    6 continue
      icf=1
      return
      end
      function icat (x,j,cm)
      implicit none
      integer icat
      integer j
      double precision x,cm(*)
      integer j0,j1,j2,k
      j0=cm(2*j)+.1
      j1=j0
      j2=cm(2*j+1)+.1
    1 if(j2.eq.j1+1) go to 5
      k=(j1+j2)/2
      if(cm(k) .ne. x) go to 2
      icat=k-j0+1
      return
    2 if(cm(k) .ge. x) go to 3
      j1=k
      go to 1
    3 j2=k
      go to 1
    5 if(x .ne. cm(j1)) go to 6
      icat=j1-j0+1
      go to 8
    6 if(x .ne. cm(j2)) go to 7
      icat=j2-j0+1
      go to 8
    7 icat=0
    8 return
      end
      subroutine csp  (jp,nc,m,n,x,y,w,nk,tb,cm,kcp,yb,d,kr,ntt,sw,me,mk
     1p2,nop,sc,db,sp,mm)
      implicit none
      integer jp,nc,m,n,nk,kcp,kr,ntt,me,mkp2,nop,mm(nc,2)
      double precision x(n,*),y(n),w(n),tb(5,nk),cm(*),sc(n)
      double precision yb,sw,d(nk,*),db(n,*),sp(mkp2,*),a,b,s,eps,dv,dy
      integer mk,mkp1,n1,j,i,k,ns,jj,nrt,k1,js,nr
      double precision big,h,wh,bof0,bof1
      data eps,big /1.d-4,9.9e30/
      nop=0
      if(nc .gt. 1) go to 1
      tb(1,m)=big
      return
    1 mk=mkp2-2
      mkp1=mk+1
      n1=nc+1
      do 3 j=1,n1
      do 2 i=1,mkp2
      sp(i,j)=0.d0
    2 continue
    3 continue
      do 4 j=1,nc
      mm(j,2)=0
    4 continue
      do 7 i=1,n
      h=sc(i)
      if(h.le.0.0.or.w(i).le.0.0) go to 7
      wh=w(i)*h
      k=x(i,jp)+.1
      mm(k,2)=mm(k,2)+1
      sp(mkp2,k)=sp(mkp2,k)+wh
      sp(mkp1,k)=sp(mkp1,k)+wh*(y(i)-yb)
      sp(m,k)=sp(m,k)+wh*h
      j=1
      go to 6
    5 j=j+1
    6 if((j).gt.(kr)) go to 7
      sp(j,k)=sp(j,k)+wh*db(i,j)
      go to 5
    7 continue
      do 8 j=1,nc
      mm(j,1)=j
    8 continue
      bof0=big
      ns=0
      jj=nc
      nrt=0
      k1=1
    9 bof1=big
      js=0
      do 18 j=1,jj
      k=mm(j,1)
      if(mm(k,2).eq.0) go to 18
      nr=nrt+mm(k,2)
      if(nr.le.me.or.ntt-nr.le.me) go to 18
      dy=sp(mkp1,n1)+sp(mkp1,k)
      a=sp(mkp2,n1)+sp(mkp2,k)
      dv=sp(m,n1)+sp(m,k)-a**2/sw
      if(dv .le. 0.d0) go to 17
      i=1
      go to 11
   10 i=i+1
   11 if((i).gt.(kr)) go to 12
      d(i,2)=sp(i,n1)+sp(i,k)
      go to 10
   12 a=0.d0
      b=a
      i=1
      go to 14
   13 i=i+1
   14 if((i).gt.(kr)) go to 15
      s=d(i,2)
      a=a+s*d(i,1)
      b=b+s**2
      go to 13
   15 b=dv-b
c      print *, 'b, eps, dv=', b, eps, dv
      if(b .le. eps*dv) go to 17
      nop=nop+1
      b=-(dy-a)**2/b
c      print *, 'b, bof0, bof1=',b, bof0, bof1
c     include a small amount for error
      if(b + 1.e-10 .ge. bof1) go to 16
      bof1=b
      js=j
   16 if(b + 1.e-10 .ge. bof0) go to 17
      bof0=b
      ns=jj
   17 continue
c      print *, b, bof0, bof1, kcp
      if(nc.eq.2) go to 19
   18 continue
   19 continue
c      print *, 'js=',js
      if(js.eq.0) go to 23
      k=mm(js,1)
      mm(js,1)=mm(jj,1)
      mm(jj,1)=k
      sp(mkp1,n1)=sp(mkp1,n1)+sp(mkp1,k)
      sp(mkp2,n1)=sp(mkp2,n1)+sp(mkp2,k)
      nrt=nrt+mm(k,2)
      sp(m,n1)=sp(m,n1)+sp(m,k)
      i=1
      go to 21
   20 i=i+1
   21 if((i).gt.(kr)) go to 22
      sp(i,n1)=sp(i,n1)+sp(i,k)
      go to 20
   22 jj=jj-1
      if(jj.le.2) go to 23
      go to 9
   23 tb(1,m)=bof0
      tb(3,m)=kcp
      do 24 j=1,nc
      cm(j+kcp)=0.0
   24 continue
      if(ns.eq.0) return
      do 25 j=ns,nc
      cm(mm(j,1)+kcp)=1.0
   25 continue
c      print *, 'in csp, kcp, ns, nc,mm=',kcp, ns, nc, (mm(j,1), j=ns,nc)
      print *, 'in csp cm=',(cm(j), j=kcp+1,kcp+nc)
      return
      end
      subroutine rspnpr (it,il,n,y,w,m)
      implicit none
      integer it,il,n,m(n)
      double precision y(n),w(n),wm(2)
      integer i,k
      double precision wt
      if(it.le.0) return
      if(il .ne. 1) go to 2
      wm(1)=0.0
      wm(2)=wm(1)
      do 1 i=1,n
      k=y(i)+1.1
      wm(k)=wm(k)+w(i)
    1 continue
      wt=wm(1)+wm(2)
      wm(1)=wm(1)/wt
      wm(2)=wm(2)/wt
      write(it,'(/,'' binary (0/1) response:  mass(0) ='',g12.4,
     1                           ''   mass(1) ='',g12.4)') wm(1),wm(2)
      return
    2 continue
      write(it,'(/,'' ordinal response:'')')
      write(it,'(''      min         n/4         n/2        3n/4
     1 max'')')
      do 3 i=1,n
      m(i)=i
    3 continue
      call psort(y,m,1,n)
      write(it,'('' '',5g12.4)') y(m(1)),y(m(n/4)),y(m(n/2)),y(m(n-n/4))
     1,y(m(n))
      return
      end
      subroutine ordpr (it,n,p,x,lx,m)
      implicit none
      integer it,n,p,lx(p),m(n,p)
      double precision x(n,p)
      integer no,j,n1,n2,n3
      if(it.le.0) return
      no=0
      do 1 j=1,p
      if(lx(j).gt.0) no=no+1
    1 continue
      if(no.eq.0) return
      write(it,'(/,'' there are'',i3,'' ordinal predictor variables.'',/
     1)') no
      write(it,'(''  var     min         n/4         n/2        3n/4
     1     max'')')
      n1=n/4
      n2=n/2
      n3=n-n1
      do 2 j=1,p
      if(lx(j).le.0) go to 2
      write(it,'('' '',i3,'' '',5g12.4)') j,  x(m(1,j),j),x(m(n1,j),j),x
     1(m(n2,j),j),x(m(n3,j),j),x(m(n,j),j)
    2 continue
      return
      end
      subroutine catpr(it,n,p,x,cm,mm)
      implicit none
      integer it,n,p,mm(*)
      double precision x(n,p),cm(*)
      integer nct,n2,np,i,j1,j2,j,nv,ic,k
      if(it.le.0) return
      nct=cm(1)+.1
      if(nct.eq.0) return
      n2=2*p+1
      np=0
      write(it,'(/,'' there are'',i3,'' categorical predictor variables.
     1'')') nct
      i=2
      go to 2
    1 i=i+(2)
    2 if((2)*((i)-(n2)).gt.0) go to 6
      np=np+1
      j1=cm(i)+.1
      if(j1.eq.0) go to 1
      j2=cm(i+1)+.1
      nv=j2-j1+1
      do 3 j=1,nv
      mm(j)=0
    3 continue
      do 4 j=1,n
      ic=x(j,np)+.1
      mm(ic)=mm(ic)+1
    4 continue
      write(it,'(/,'' categorical variable'',i3,'' has'',i3,'' values.''
     1)') np,nv
      write(it,'(''  value     internal code     counts'')')
      k=0
      do 5 j=j1,j2
      k=k+1
      write(it,'(f6.0,i13,i15)') cm(j),k,mm(k)
    5 continue
      go to 1
    6 return
      end
      subroutine holl (jp,cm,t,h)
      implicit none
      integer jp
      double precision t,cm(*)
      character*28 h
      integer j1,j2,k,j
      j1=cm(2*jp)+.1
      j2=cm(2*jp+1)+.1
      j2=j2-j1+1
      if(j2 .le. 28) go to 1
      h='   cat. factor > 28 values  '
      return
    1 h='                            '
      j1=(28-j2)/2
      j2=j1+j2-1
      k=t+.1
c      print *, 't, k, j1, j2 =', t, k, j1, j2
      do 3 j=j1,j2
c      print *, 'j=',j,cm(k+j-j1+1)
      if(cm(k+j-j1+1) .le. 0.0) go to 2
      h(j:j)='1'
      go to 3
    2 h(j:j)='0'
    3 continue
c      print *, h
      return
      end
      subroutine slova (nk,it,tb,ni,lp,lv)
      implicit none
      integer nk,it,ni,lp(3,*),lv(*)
      double precision tb(5,nk)
      integer m,na,k2,i2
c     write(it,4) ni
      call coll(nk,tb,lp,lv,lp(1,nk+1))
      m=1
    1 if(lp(1,m).eq.0) go to 2
      m=m+1
      go to 1
    2 na=m-1
      do 3 m=1,na
      k2=lp(2,m)
      i2=lp(1,m)+k2-1
c     write(it,5) m,lp(3,m),(lv(i),i=k2,i2)
    3 continue
      return
    4 format(/,' sliced anova decomposition on',i3,' basis functions:',/
     1,  '   fun      #bsfns      variable(s)')
    5 format('  ',i3,'         ',i2,'       ',20i4)
      end
      subroutine reducq (flg,x,nk,tb,cm,tc,kp,kv,lp,lv,r,td,sc,fc)
      implicit none
      integer nk,kp(5,*),kv(2,*),lp(3,*),lv(*)
      double precision flg,x(*),tb(5,nk),cm(*),tc(*),r(*),td(2,nk),sc(2,
     1*),fc(*)
      integer ll,la,laa,m,nv,jl,il,k,nt,jp,kp3,l1,l
      integer match
      ll=1
      la=ll
      l1=la
      laa=0
      do 1 m=1,nk
      td(1,m)=0.0
    1 continue
    2 if(kp(1,ll).lt.0) go to 9
      nv=0
      if(kp(1,ll) .le. 0) go to 4
      jl=kp(1,ll)
      do 3 il=1,jl
      k=kp(2,ll)+il-1
      nv=nv+1
      sc(1,nv)=kv(1,k)
      sc(2,nv)=kv(2,k)
    3 continue
      go to 5
    4 if(kp(3,ll) .gt. 0) go to 5
      ll=ll+1
      go to 2
    5 if(kp(3,ll) .gt. 0) go to 6
      m=match(nv,sc,nk,tb,cm,r,0)
      td(1,m)=tc(-kp(3,ll))
      ll=ll+1
      go to 2
    6 kp3=kp(3,ll)
      do 8 k=1,kp3
      l=lp(1,l1)
      nt=lp(3,l1)
      laa=laa+5*l*nt
      do 7 jp=1,nt
      call gtrm(1,jp,l,nt,lv(lp(2,l1)),flg,x,nk,tb,tc(la),sc(1,nv+1),fc)
      m=match(nv+l,sc,nk,tb,cm,r,0)
      td(1,m)=tc(jp+laa)
      call std(m,flg,x,l,sc(1,nv+1),fc,nk,tb,r,td)
    7 continue
      laa=laa+nt
      l1=l1+1
      la=la+nt*(5*l+1)
    8 continue
      ll=ll+1
      go to 2
    9 return
      end
      subroutine gtrm (il,jp,l,nt,jv,flg,x,nk,tb,tc,te,fc)
      implicit none
      integer il,jp,l,nt,nk,jv(l)
      double precision flg,x(*),tb(5,nk),tc(nt,*),te(2,*),fc(*)
      integer l2,nf,l3,l4,k,j,jj
      double precision cue
      l2=l+l
      nf=0
      l3=l2+l
      l4=l3+l
      do 1 k=1,l
      j=jv(k)
      jj=j
      if(tc(jp,k+l).gt.tc(jp,k+l2)) jj=-jj
      te(1,k)=jj
      te(2,k)=tc(jp,k)
      if(il.eq.2) go to 1
      if(x(j).eq.flg) go to 1
      nf=nf+1
      fc(nf)=cue(x(j),tc(jp,k+l),tc(jp,k),tc(jp,k+l2),tc(jp,k+l3),tc(jp,
     1k+l4))
    1 continue
      return
      end
      function match (nv,te,nk,tb,cm,r,iz)
      implicit none
      integer match
      integer nv,nk,iz
      double precision te(2,nv),tb(5,nk),cm(*),r(*)
      integer m,jg,j,jp,jp2,jp21,jpp,ig,ip,jqq,i1,
     1i2,kg,j1,j2,jq,nc,i
      double precision t,u
      integer nord,ieq
      match=0
      do 15 m=1,nk
      if(tb(1,m).eq.0.0) go to 15
      if(nord(m,tb).ne.nv) go to 15
      jg=0
      do 13 j=1,nv
      t=te(1,j)
      u=te(2,j)
      jp=abs(t)+.1
      jp2=2*jp
      jp21=jp2+1
      jpp=jp
      if(t.lt.0.0) jpp=-jpp
      ig=0
      ip=m
    1 if(ip.le.0) go to 12
      t=tb(2,ip)
      jq=abs(t)+.1
      jqq=jq
      if(t.lt.0.0) jqq=-jqq
      if(jp .eq. jq) go to 2
      ip=tb(4,ip)+.1
      go to 1
    2 if(cm(jp2) .ne. 0.0) go to 4
      if(jpp .ne. jqq .or. ieq(tb(3,ip),u,r(jp)) .ne. 1) go to 3
      ig=1
      go to 12
    3 ip=tb(4,ip)+.1
      go to 1
    4 nc=int(cm(jp21)+.1)-int(cm(jp2)+.1)+1
      i1=u+.1
      i2=tb(3,ip)+.1
      kg=0
      do 9 i=1,nc
      j1=cm(i1+i)
      j2=cm(i2+i)
      if(jpp .ge. 0) go to 6
      if(j1 .ne. 0) go to 5
      j1=1
      go to 6
    5 j1=0
    6 if(jqq .ge. 0) go to 8
      if(j2 .ne. 0) go to 7
      j2=1
      go to 8
    7 j2=0
    8 if(j1 .eq. j2) go to 9
      kg=1
      go to 10
    9 continue
   10 if(kg .ne. 0) go to 11
      ig=1
      go to 12
   11 ip=tb(4,ip)+.1
      go to 1
   12 if(ig .ne. 0) go to 13
      jg=1
      go to 14
   13 continue
   14 if(jg .ne. 0) go to 15
      match=m
      go to 16
   15 continue
   16 if(match.gt.0.or.iz.ne.0) return
c     write(6,17)
   17 format (' bug in match - term not found.')
      do 19 j=1,nv
c     write(6,18)j,te(1,j),te(2,j)
   18 format (' te(',i2,')=',2g12.4)
   19 continue
      do 21 j=1,nk
c     write(6,20)j,(tb(i,j),i=1,4)
   20 format (' tb(',i2,')=',4g12.4)
   21 continue
      stop
      end
      subroutine std (m,flg,x,l,te,fc,nk,tb,r,td)
      implicit none
      integer m,l,nk
      double precision flg,x(*),te(2,*),fc(*),tb(5,nk),r(*),td(2,nk)
      integer ip,j,jj,k,i,ig,jp
      double precision t,u
      integer ieq
      ip=m
    1 if(ip.le.0) go to 4
      t=tb(2,ip)
      j=abs(t)+.1
      jj=j
      if(t.lt.0.0) jj=-jj
      u=tb(3,ip)
      k=0
      ig=0
      do 2 i=1,l
      t=te(1,i)
      jp=abs(t)+.1
      if(x(jp).ne.flg) k=k+1
      if(t.lt.0.0) jp=-jp
      if(jj .ne. jp .or. ieq(te(2,i),u,r(j)) .ne. 1) go to 2
      ig=1
      go to 3
    2 continue
    3 if(ig.eq.1.and.x(j).ne.flg) td(2,ip)=fc(k)
      ip=tb(4,ip)+.1
      go to 1
    4 return
      end
      subroutine reducl (flg,x,nk,az,tb,cm,bz,td,r,azn,tbn,bzn,sc)
      implicit none
      integer nk
      double precision flg,az,bz,azn,bzn,x(*),tb(5,nk),cm(*),td(2,nk),r(
     1*),tbn(5,nk),sc(*)
      integer m,i,j,k,ip,iq,no,nv
      double precision t,u
      integer icat,match
      azn=az
      do 2 m=1,nk
      do 1 i=1,5
      tbn(i,m)=tb(i,m)
    1 continue
    2 continue
      bzn=bz
      do 9 m=1,nk
      t=tb(2,m)
      j=abs(t)+.1
      if(x(j).eq.flg) go to 9
      if(cm(2*j) .le. 0.0) go to 7
      k=icat(x(j),j,cm)
      if(k .ne. 0) go to 3
      u=0.0
      go to 4
    3 u=cm(k+int(tb(3,m)+.1))
    4 if(t .ge. 0.0) go to 6
      if(u .ne. 0.0) go to 5
      u=1.0
      go to 6
    5 u=0.0
    6 td(2,m)=u
      go to 8
    7 u=dmax1(0.0D0,sign(1.0D0,t)*(x(j)-tb(3,m)))
    8 sc(m)=u
    9 continue
      m=nk
      go to 11
   10 m=m+(-1)
   11 if((-1)*((m)-(1)).gt.0) go to 21
      ip=tbn(4,m)+.1
      t=tbn(2,m)
      j=abs(t)+.1
      if(x(j) .ne. flg) go to 15
      if(tbn(1,m) .eq. 0.0) go to 10
      iq=ip
   12 if(iq.le.0) go to 10
      t=tbn(2,iq)
      j=abs(t)+.1
      if(x(j) .eq. flg) go to 13
      tbn(1,m)=tbn(1,m)*sc(iq)
      td(1,m)=td(1,m)*td(2,iq)
   13 iq=tbn(4,iq)+.1
      go to 12
   15 k=m+1
      go to 17
   16 k=k+1
   17 if((k).gt.(nk)) go to 18
      if(int(tbn(4,k)+.1).eq.m) tbn(4,k)=tbn(4,m)
      go to 16
   18 if(tbn(1,m).eq.0.0) go to 10
      if(ip .ne. 0) go to 19
      azn=azn+tbn(1,m)*sc(m)
      bzn=bzn+td(1,m)*td(2,m)
      go to 10
   19 tbn(1,ip)=tbn(1,ip)+tbn(1,m)*sc(m)
      td(1,ip)=td(1,ip)+td(1,m)*td(2,m)
      go to 10
   21 no=nk
      m=nk
      go to 23
   22 m=m+(-1)
   23 if((-1)*((m)-(1)).gt.0) go to 31
      t=tb(2,m)
      j=abs(t)+.1
      if(x(j).eq.flg) go to 22
      k=m+1
      go to 25
   24 k=k+1
   25 if((k).gt.(no)) go to 30
      td(1,k-1)=td(1,k)
      do 26 i=1,5
      tbn(i,k-1)=tbn(i,k)
   26 continue
      i=k+1
      go to 28
   27 i=i+1
   28 if((i).gt.(no)) go to 24
      if(int(tbn(4,i)+.1).eq.k) tbn(4,i)=k-1
      go to 27
   30 no=no-1
      go to 22
   31 m=no+1
      go to 33
   32 m=m+1
   33 if((m).gt.(nk)) go to 34
      tbn(1,m)=0.0
      go to 32
   34 m=no
      go to 36
   35 m=m+(-1)
   36 if((-1)*((m)-(2)).gt.0) go to 39
      if(tbn(1,m).eq.0.0) go to 35
      nv=0
      ip=m
   37 if(ip.le.0) go to 38
      nv=nv+1
      sc(2*nv-1)=tbn(2,ip)
      sc(2*nv)=tbn(3,ip)
      ip=tbn(4,ip)+.1
      go to 37
   38 k=match(nv,sc,m-1,tbn,cm,r,1)
      if(k.eq.0) go to 35
      tbn(1,k)=tbn(1,k)+tbn(1,m)
      td(1,k)=td(1,k)+td(1,m)
      tbn(1,m)=0.0
      go to 35
   39 return
      end
      subroutine qslice (p,nk,tb,cm,td,kp,kv,lp,lv,tc,r,sc,js)
      implicit none
      integer p,nk,kp(5,*),kv(2,*),lp(3,*),lv(*),js(*)
      double precision tb(5,nk),cm(*),td(2,*),tc(*),r(p,2),sc(2,p)
      integer j,ll,kp3,m,l,nt,l1,la,le,laa,nv,jl,il,k,
     1jp
      double precision dum
      integer match
      do 1 j=1,p
      sc(1,j)=r(j,2)
      sc(2,j)=sc(1,j)+r(j,1)
    1 continue
      ll=1
      la=ll
      l1=la
    2 if(kp(1,ll).lt.0) go to 5
      if(kp(3,ll) .gt. 0) go to 3
      kp(5,ll)=0
      ll=ll+1
      go to 2
    3 kp3=kp(3,ll)
      kp(5,ll)=la
      do 4 m=1,kp3
      l=lp(1,l1)
      nt=lp(3,l1)
      call knts(l,nt,lv(lp(2,l1)),kp(1,ll),kv(1,kp(2,ll)),nk,tb,cm,tc(la
     1),js)
      call side(l,nt,lv(lp(2,l1)),sc,tc(la))
      l1=l1+1
      la=la+nt*(5*l+1)
    4 continue
      ll=ll+1
      go to 2
    5 le=la-1
      ll=1
      la=ll
      l1=la
      laa=0
    6 if(kp(1,ll).lt.0) go to 13
      nv=0
      if(kp(1,ll) .le. 0) go to 8
      jl=kp(1,ll)
      do 7 il=1,jl
      k=kp(2,ll)+il-1
      nv=nv+1
      sc(1,nv)=kv(1,k)
      sc(2,nv)=kv(2,k)
    7 continue
      go to 9
    8 if(kp(3,ll) .gt. 0) go to 9
      ll=ll+1
      go to 6
    9 if(kp(3,ll) .gt. 0) go to 10
      m=match(nv,sc,nk,tb,cm,r,0)
      le=le+1
      kp(3,ll)=-le
      tc(le)=td(1,m)
      ll=ll+1
      go to 6
   10 kp3=kp(3,ll)
      do 12 k=1,kp3
      l=lp(1,l1)
      nt=lp(3,l1)
      laa=laa+5*l*nt
      do 11 jp=1,nt
      call gtrm(2,jp,l,nt,lv(lp(2,l1)),dum,dum,nk,tb,tc(la),sc(1,nv+1),d
     1um)
      m=match(nv+l,sc,nk,tb,cm,r,0)
      tc(jp+laa)=td(1,m)
   11 continue
      laa=laa+nt
      l1=l1+1
      la=la+nt*(5*l+1)
   12 continue
      ll=ll+1
      go to 6
   13 return
      end
      function ieq(a,b,r)
      implicit none
      integer ieq
      double precision a,b,r
      ieq=0
      if(abs((a-b)/r).lt.1.e-5) ieq=1
      return
      end
      function lcm (p,nk,tb,cm)
      implicit none
      integer lcm
      integer p,nk
      double precision tb(5,nk),cm(*)
      integer m,ix,j,jj
      ix=0
      jj=0
      do 1 m=1,nk
      j=abs(tb(2,m))+.1
      if(cm(2*j).eq.0.0) go to 1
      if(int(tb(3,m)+.1) .le. ix) go to 1
      ix=tb(3,m)+.1
      jj=j
    1 continue
      if(ix .le. 0) go to 2
      lcm=ix+int(cm(2*jj+1)+.1)-int(cm(2*jj)+.1)+1
      return
    2 lcm=2*p+1
      do 3 j=1,p
      if(cm(2*j).eq.0.0) go to 3
      lcm=lcm+int(cm(2*j+1)+.1)-int(cm(2*j)+.1)+1
    3 continue
      return
      end
      function newb (m,tb)
      implicit none
      integer newb
      integer m
      double precision tb(5,m)
      integer mm1,k
      integer ieq
      newb=0
      mm1=m-1
      do 1 k=1,mm1
      if(ieq(tb(2,k),tb(2,m),1.0D0).eq.0) go to 1
      if(ieq(tb(3,k),tb(3,m),1.0D0).eq.0) go to 1
      if(ieq(tb(4,k),tb(4,m),1.0D0).eq.0) go to 1
      newb=1
      go to 2
    1 continue
    2 return
      end
      subroutine sscp (n,m,sc,y,w,yb,yv,sw,d,da)
      implicit none
      integer n,m
      double precision sc(n,*),y(n),w(n)
      double precision d(m,m),da(*),yb,yv,sw,s
      integer k,mm1,i,j, ii
      mm1=m-1
c      print *, 'sc=', mm1
c      do 999 k=1,mm1
c      print *, ' k=', k
c 999  print *,(sc(ii,k), ii=1,n)
      do 6 k=1,mm1
      s=0.d0
      do 1 i=1,n
      s=s+w(i)*sc(i,k)
    1 continue
      s=s/sw
c      print *, '   s=',s
      da(k)=s
      do 2 i=1,n
      sc(i,k)=sc(i,k)-s
    2 continue
      do 4 j=1,k
      s=0.d0
      do 3 i=1,n
      s=s+w(i)*sc(i,j)*sc(i,k)
    3 continue
      d(j,k)=s
    4 continue
      s=0.d0
      do 5 i=1,n
      s=s+w(i)*sc(i,k)*(y(i)-yb)
    5 continue
      d(k,m)=s
    6 continue
      d(m,m)=sw*yv
      return
      end
      subroutine lsf1 (d,m,xb,yb,al,rss,a,a0,dp)
      implicit none
      integer m
      double precision d(m,m),xb(*),yb,al,rss,a(*),a0,dp(*),eps,s
      integer i,mm1,im1,j, ii, jj
      data eps /1.d-4/
      mm1=m-1
      do 1 i=1,mm1
      dp(i)=d(i,i)
      d(i,i)=d(i,i)*(1.d0+al)
    1 continue
      do 5 i=1,mm1
      if(dp(i).le.0.d0) go to 5
      im1=i-1
      s=dp(i)
      j=1
      go to 3
    2 j=j+1
    3 if((j).gt.(im1)) go to 4
      if(d(j,j).lt.0.d0) s=s+dp(j)*d(j,i)**2
      go to 2
    4 if((d(i,i)-al*s)/dp(i).lt.eps) go to 5
      call sweep(d,m,i,-1.d0,dp(m))
    5 continue

      rss=0.d0
      a0=yb
      do 6 i=1,mm1
      a(i)=0.d0
      if(d(i,i).ge.0.d0) go to 6
      a(i)=d(i,m)
      a0=a0-a(i)*xb(i)
      rss=rss+dp(i)*a(i)**2
    6 continue
      rss=d(m,m)-al*rss
      return
      end
      subroutine bkstp (d,m,xb,yb,al,rss,a,a0,k,dp)
      implicit none
      integer m,k
      double precision d(m,m),xb(*),yb,al
      double precision a(*),a0,rss,dp(*),s
      integer mm1,i,j
      double precision big
      data big /9.9e30/
      mm1=m-1
      rss=big
      k=0
      do 4 i=1,mm1
      if(d(i,i).ge.0.d0) go to 4
      s=0.d0
      do 3 j=1,mm1
      if(d(j,j).ge.0.d0) go to 3
      if(j.eq.i) go to 3
      if(j .ge. i) go to 1
      a0=d(j,i)
      go to 2
    1 a0=d(i,j)
    2 s=s+dp(j)*(d(j,m)-a0*d(i,m)/d(i,i))**2
    3 continue
      s=d(m,m)-d(i,m)**2/d(i,i)-al*s
      if(s .gt. rss) go to 4
      rss=s
      k=i
    4 continue
      if(k.gt.0) call sweep(d,m,k,1.d0,dp(m))
      a0=yb
      rss=0.d0
      do 5 i=1,mm1
      a(i)=0.d0
      if(d(i,i).ge.0.d0) go to 5
      a(i)=d(i,m)
      a0=a0-a(i)*xb(i)
      rss=rss+dp(i)*a(i)**2
    5 continue
      rss=d(m,m)-al*rss
      return
      end
      subroutine sweep (a,m,k,fl,u)
      implicit none
      integer m,k
      double precision a(m,m),u(m),fl,c
      integer i,j
      c=a(k,k)
      do 1 i=1,k
      u(i)=a(i,k)
      a(i,k)=0.d0
    1 continue
      do 2 i=k,m
      u(i)=a(k,i)
      a(k,i)=0.d0
    2 continue
      u(k)=fl
      do 4 i=1,m
      do 3 j=i,m
      a(i,j)=a(i,j)-u(i)*u(j)/c
    3 continue
    4 continue
      return
      end
      subroutine array(p,n,i,j)
      implicit none
      integer p,n,i,j
      i=mod(p,n)
      if(i.eq.0) i=n
      j=(p-i)/n+1
      return
      end
      subroutine cvmars  (ix,n,p,x,y,w,nk,ms,df,fv,mi,lx,it,xm,xs,tb,cm,
     1sc,db,d,mm,wt,cv)
      implicit none
      integer ix,n,p,nk,ms,it,mm(n,*),lx(p)
      double precision x(n,p),y(n),w(n),xm(p),xs(p),tb(5,nk),cm(*),sc(*)
     1,wt(n,2),cv(nk,4)
      double precision db(n,*),d(nk,*)
      integer i,im,nr,nd,k,ir,m,mi,mk,ia1, ii, jj
      double precision eps,big,dfs,cvm,t,cv0,sw,wn,yv,r,fc,df,fv,yv1,wn1
     1,am,am1,az,dmx,cvl,dfu,gcv,a1,a2
      data eps,big,dfs,cvm,im /1.e-6,9.9e30,2*0.0,0/
      if(it.gt.0) write(it,'(/,'' sample reuse to estimate df:'')')
      if(ix .le. 0) go to 1
      nr=ix
      nd=nr
      if(it.gt.0) write(it,'('' '',i3,'' - fold cross-validation.'',/)')
     1 ix
      go to 2
    1 nr=1
      nd=-ix
      if(it.gt.0) write(it,'('' independent test set - every'',i4,'' obs
     1ervations.'',/)') nd
    2 do 3 i=1,n
      wt(i,1)=w(i)
      wt(i,2)=i
    3 continue
      do 4 i=1,n
      call rnms(r,1)
      k=(n-i+1)*r+i
      t=wt(i,2)
      wt(i,2)=wt(k,2)
      wt(k,2)=t
    4 continue
      do 6 i=1,3
      do 5 m=1,nk
      cv(m,i)=0.0
    5 continue
    6 continue
      cv0=0.0
      sw=cv0
      wn=sw
      yv=wn
      fc=yv
      do 14 ir=1,nr
      i=ir
    7 if(i.gt.n) go to 8
      wt(int(wt(i,2)+.1),1)=0.0
      i=i+nd
      go to 7
    8 continue
c      print *, 'before marsgo'
c      do 71 ii=1,5
c   71 print *,(tb(ii,jj),jj=1,nk)
c      do 711 ii=1,n
c  711 print *,(mm(ii,jj),jj=1,p)
      call marsgo (n,p,x,y,wt,nk,ms,df,fv,mi,lx,99,xm,xs,az,tb,cm,sc,db,
     1d,mm)
      yv1=sc(3)
      yv=yv+yv1
      wn1=sc(2)
      wn=wn+wn1
      fc=fc+sc(1)
      mk=sc((nk+1)**2+4)+.1
      i=ir
    9 if(i.gt.n) go to 10
      k=wt(i,2)+.1
      wt(k,1)=w(k)
      sw=sw+w(k)
      call cvmod(k,n,x,y,w,nk,mk,tb,cm,sc,cv0,cv(1,3))
      i=i+nd
      go to 9
   10 continue
      do 13 m=1,nk
      am=sc(m+4)
      cv(m,2)=cv(m,2)+am
      am1=yv1
      if(m.gt.1) am1=sc(m+3)
      if(am1/yv1 .le. eps) go to 11
      r=dsqrt(am/am1)
      go to 12
   11 r=1.0
   12 continue
      cv(m,1)=cv(m,1)+((wn1-1.0)*(1.0-r)/(m-r*(m-1))-1.0)/sc(1)
   13 continue
   14 continue
      do 15 m=1,nk  
      cv(m,1)=cv(m,1)/nr
      cv(m,2)=cv(m,2)/nr
      cv(m,3)=cv(m,3)/sw
   15 continue
      fc=fc/nr
      yv=yv/nr
      wn=wn/nr
      cv0=cv0/sw
      if(it.gt.0) write(it,21)
      im=0
      cvm=cv0
      dmx=-big
      cvl=cv(nk,1)
      m=nk
      go to 17
   16 m=m+(-1)
   17 if((-1)*((m)-(1)).gt.0) go to 19
      if(cv(m,1).le.dmx) go to 16
      dmx=cv(m,1)
      dfu=0.5*(cvl+cv(m,1))
      cvl=cv(m,1)
      if(cv(m,3) .gt. cvm) go to 18
      cvm=cv(m,3)
      df=dfu
      im=m
   18 continue
      gcv=cv(m,2)/(1.0-((dfu*fc+1.0)*m+1.0)/wn)**2
      if(it.gt.0) write(it,22) m,dfu,cv(m,2),gcv,cv(m,3)
      go to 16
   19 if(cv0 .gt. cvm) go to 20
      cvm=cv0
      df=dmx
      im=0
   20 dfs=df
      gcv=yv/(1.0-1.0/wn)**2
      if(it.gt.0) write(it,22) 0,dmx,yv,gcv,cv0
      if(it.gt.0) write(it,'(/,'' estimated optimal df('',i3,'') ='',f7.
     12,               '' with (estimated) pse ='',g12.4)') im,df,cvm
      return
      entry cvinfo(a1,a2,ia1)
      a1=dfs
      a2=cvm
      ia1=im
      return
   21 format('  #bsfns     df        asr           gcv           cv')
   22 format(' ',i5,d10.2,3g14.4)
      end
      subroutine cvmod (i,n,x,y,w,nk,mk,tb,cm,sc,cv0,cv)
      implicit none
      integer i,n,nk,mk
      double precision x(n,*),y(n),w(n),tb(5,nk),cm(*),sc(*),cv(nk,2),
     1cv0
      integer m,j,k,l,kp
      double precision t,s,u
      do 8 m=1,mk
      t=tb(2,m)
      j=abs(t)+.1
      if(cm(2*j) .le. 0.0) go to 5
      k=x(i,j)+.1
      if(k .ne. 0) go to 1
      u=0.0
      go to 2
    1 u=cm(k+int(tb(3,m)+.1))
    2 if(t .ge. 0.0) go to 6
      if(u .ne. 0.0) go to 3
      u=1.0
      go to 6
    3 u=0.0
      go to 6
    5 u=dmax1(0.0D0,sign(1.0D0,t)*(x(i,j)-tb(3,m)))
    6 l=tb(4,m)+.1
      if(l .le. 0) go to 7
      cv(m,2)=u*cv(l,2)
      go to 8
    7 cv(m,2)=u
    8 continue
      kp=nk+4
      cv0=cv0+w(i)*(y(i)-sc(4))**2
      do 10 m=1,nk
      kp=kp+1
      s=sc(kp)
      do 9 l=1,nk
      kp=kp+1
      if(l.le.mk) s=s+sc(kp)*cv(l,2)
    9 continue
      cv(m,1)=cv(m,1)+w(i)*(y(i)-s)**2
   10 continue
      return
      end
      subroutine nest (n,i,j,nv,vals)
      implicit none
      integer n,i,j,nv
      integer mlist,nlist
      parameter (mlist=200, nlist=2000)
      double precision vals(*),vm(nlist),tb(5,*),cm(*),x(n,*),bl(*)
      integer m(4,mlist),p,lx(*)
      save m,vm
      integer il,jl,ig,k,j1,jn,jv,jp,lm,nc,ll,l,it,k1,k2,
     1kk,ja,jb,ip,kp,lp,lk,mk,kg,lon,norm,kx,jg
      double precision t,ex,t1
      integer nord,ieqbf,icat
      data il,jl /2*0/
      if((i .ne. 0) .and. (j .ne. 0)) go to 1
      il=0
      jl=il
      return
    1 if(i.eq.j) return
      ig=0
      if(nv .le. 0) go to 8
      k=1
      go to 3
    2 k=k+1
    3 if((k).gt.(il)) go to 4
      if(m(1,k).eq.i.or.m(1,k).eq.j) return
      go to 2
    4 il=il+1
      if(il .le. mlist) go to 5
c     write(6,  '('' increase parameter mlist in subroutine nest to grea
c    1ter than'',               i5,'' and recompile.'')') il
      stop
    5 m(1,il)=i
      m(2,il)=j
      m(3,il)=nv
      m(4,il)=jl
      if(jl+nv .le. nlist) go to 6
c     write(6,  '('' increase parameter nlist in subroutine nest to grea
c    1ter than'',               i5,'' and recompile.'')') jl+nv
      stop
    6 do 7 k=1,nv
      jl=jl+1
      vm(jl)=vals(k)
    7 continue
      return
    8 k=1
      go to 10
    9 k=k+1
   10 if((k).gt.(il)) go to 12
      if(m(1,k) .ne. i .or. m(2,k) .ne. j) go to 9
      ig=1
   12 if(ig.eq.0) return
      il=il-1
      ll=k
      go to 14
   13 ll=ll+1
   14 if((ll).gt.(il)) go to 16
      do 15 l=1,4
      m(l,ll)=m(l,ll+1)
   15 continue
      go to 13
   16 return
      entry nstlst(it)
      if(it.le.0) return
      if(il.eq.0) return
c     write(it,'(/,'' variable nesting:'',/)')
      do 18 k=1,il
      if(m(3,k) .le. 5) go to 17
c     write(it,'('' '',i3,'': var('',i3,'') exists for var('',i3,'') =''
c    1)')  k,m(1,k),m(2,k)
c     write(it,'(100('' '',10f7.1))') (vm(l),l=m(4,k)+1,m(4,k)+m(3,k))
      go to 18
   17 continue
c     write(it,'('' '',i3,'': var('',i3,'') exists for var('',i3,'') =''
c    1,5f7.1)')  k,m(1,k),m(2,k),(vm(l),l=m(4,k)+1,m(4,k)+m(3,k))
   18 continue
      return
      entry oknest(it,p,lx,cm)
      if(it.le.0) return
      l=1
      go to 20
   19 l=l+1
   20 if((l).gt.(il)) go to 24
      j1=m(1,l)
      jn=m(2,l)
      jv=m(3,l)
      jp=m(4,l)
c     if(j1.lt.1.or.j1.gt.p) write(it,25) l,j1
c     if(jn.lt.1.or.jn.gt.p) write(it,26) l,jn
c     if(lx(jn).ge.0) write(it,27) l,jn,lx(jn)
      k1=cm(2*jn)+.1
      k2=cm(2*jn+1)+.1
      do 23 k=jp+1,jp+jv
      ig=0
      do 21 kk=k1,k2
      if(vm(k) .ne. cm(kk)) go to 21
      ig=1
      go to 22
   21 continue
   22 continue
c     if(ig.eq.0) write(it,28) l,vm(k),jn
   23 continue
      go to 19
   24 return
   25 format(' nesting entry',i3,', invalid variable',i3,' to be nested.
     1')
   26 format(' nesting entry',i3,', invalid nesting variable',i3,'.')
   27 format(' nesting entry',i3,', lx(',i3,') =',i2,'. must be < 0.')
   28 format(' nesting entry',i3,', categorical value ',g12.4,/,  ' not
     1among the data values for variable',i3,'.')
      entry isnstr(j,jb)
      jb=0
      k=1
      go to 30
   29 k=k+1
   30 if((k).gt.(il)) go to 32
      if(m(2,k) .ne. j) go to 29
      jb=m(1,k)
   32 return
      entry isfac (lm,j,mk,tb,cm,ja)
      ja=0
      ig=ja
      l=1
      go to 34
   33 l=l+1
   34 if((l).gt.(il)) go to 36
      if(j .ne. m(1,l)) go to 33
      ig=1
   36 if(ig.eq.0) return
      jn=m(2,l)
      if(cm(2*jn).eq.0.0) return
      jv=m(3,l)
      jp=m(4,l)
      ig=0
      ip=lm
   37 if(ip.le.0) go to 39
      j1=abs(tb(2,ip))+.1
      if(j1 .ne. jn) go to 38
      ig=1
      go to 39
   38 ip=tb(4,ip)+.1
      go to 37
   39 if(ig .eq. 0) go to 45
      nc=cm(2*jn+1)-cm(2*jn)+1.1
      t=tb(2,ip)
      kp=tb(3,ip)+.1
      do 44 l=1,nc
      lp=l+kp
      if(t .le. 0) go to 40
      if(cm(lp).eq.0.0) go to 44
      go to 41
   40 if(cm(lp).ne.0.0) go to 44
   41 ex=cm(int(cm(2*jn)+.1)+l-1)
      ig=0
      do 42 k=jp+1,jp+jv
      if(ex .ne. vm(k)) go to 42
      ig=1
      go to 43
   42 continue
   43 if(ig .ne. 0) go to 44
      ja=-1
      return
   44 continue
      return
   45 ja=l
      norm=nord(lm,tb)+1
      nc=cm(2*jn+1)-cm(2*jn)+1.1
      do 56 lk=1,mk
      if(nord(lk,tb).ne.norm) go to 56
      jg=0
      ip=lk
   46 if(ip.le.0) go to 55
      t1=tb(2,ip)
      j1=abs(t1)+.1
      if(j1 .ne. jn) go to 54
      kp=tb(3,ip)+.1
      kg=0
      do 52 l=1,nc
      lp=l+kp
      lon=cm(lp)+.1
      if(t1 .ge. 0.0) go to 48
      if(lon .ne. 0) go to 47
      lon=1
      go to 48
   47 lon=0
   48 ex=cm(int(cm(2*jn)+.1)+l-1)
      ig=0
      do 49 k=jp+1,jp+jv
      if(ex .ne. vm(k)) go to 49
      ig=1
      go to 50
   49 continue
   50 if(lon .ne. 1 .or. ig .ne. 0) go to 51
      kg=1
      go to 53
   51 if(lon .ne. 0 .or. ig .ne. 1) go to 52
      kg=1
      go to 53
   52 continue
   53 if(kg .ne. 0) go to 54
      jg=1
      go to 55
   54 ip=tb(4,ip)+.1
      go to 46
   55 if(jg.eq.0) go to 56
      if(ieqbf(lk,lm,tb,cm) .ne. 1) go to 56
      ja=-1
      return
   56 continue
      return
      entry cmpnst(ja,n,x,cm,bl)
      jn=m(2,ja)
      jv=m(3,ja)
      jp=m(4,ja)
      do 59 l=1,n
      kx=x(l,jn)+.1
      ex=cm(int(cm(2*jn)+.1)+kx-1)
      ig=0
      do 57 k=jp+1,jp+jv
      if(ex .ne. vm(k)) go to 57
      ig=1
      go to 58
   57 continue
   58 if(ig.eq.1) go to 59
      bl(l)=0.0
   59 continue
      return
      entry getnst(ja,cm,j,nv,vals)
      j=m(2,ja)
      jv=m(3,ja)
      jp=m(4,ja)
      nv=cm(2*j+1)-cm(2*j)+1.1
      do 60 k=1,nv
      vals(k)=0.0
   60 continue
      do 61 l=jp+1,jp+jv
      k=icat(vm(l),j,cm)
      if(k.gt.0) vals(k)=1.0
   61 continue
      return
      end
      subroutine blf0(l,ja,n,x,w,cm,sc,nnt,bl)
      implicit none
      integer l,ja,n,nnt
      double precision x(n,*),w(n),cm(*),sc(n,*),bl(n)
      integer i
      nnt=0
      call blf(l,n,sc,bl)
      if(ja.gt.0) call cmpnst(ja,n,x,cm,bl)
      do 1 i=1,n
      if(bl(i).gt.0.0.and.w(i).gt.0.0) nnt=nnt+1
    1 continue
      return
      end
      subroutine blf(l,n,sc,bl)
      implicit none
      integer l,n
      double precision sc(n,*),bl(n)
      integer i
      if(l .gt. 0) go to 2
      do 1 i=1,n
      bl(i)=1.0
    1 continue
      go to 4
    2 do 3 i=1,n
      bl(i)=sc(i,l)
    3 continue
    4 return
      end
      subroutine mnspan(ms,alf,nep,nnt,mn,me,mel)
      implicit none
      integer ms,nep,nnt,mn,me,mel
      double precision alf
      double precision al2,al25
      parameter(al2=0.693147,al25=1.732868)
      integer nst,nnr,nnl
      double precision allf,fmn,fme
      allf=-dlog(1.0-alf)
      fmn=-dlog(allf/(nep*nnt))/al25
      fme=-dlog(alf*0.125/nep)/al2
      
      if(ms .le. 0) go to 1
      me=ms*fme/fmn+0.5
      mn=ms
      go to 2
    1 me=fme+0.5
      mn=fmn+0.5
    2 me=max0(me,mn,2)
      nst=nnt-2*me-1
      nnr=nst/mn
      nnl=nst-nnr*mn
      nnr=(nnr+1)*mn-nst
      nst=min0(nnl,nnr)
      if(nnl .gt. nnr) go to 3
      nnl=1
      go to 4
    3 nnl=-1
    4 continue
      nnr=nst/2
      mel=me
      me=me+nnl*nnr
      mel=mel+nnl*nnr
      if(mod(nst,2).ne.0) mel=mel+nnl
      return
      end
      function ieqbf(lk,lm,tb,cm)
      implicit none
      integer ieqbf
      integer lk,lm
      double precision tb(5,*),cm(*)
      integer ipo,lg,jg,ko,nc,kp,lo,lon,lop,jo,ic,
     1ip,j1,kg,l,lp
      double precision to,t,t1
      integer ieq
      ipo=lm
      lg=0
      ko=0
      nc=0
    1 if(ipo.le.0) go to 16
      to=tb(2,ipo)
      jo=abs(to)+.1
      jg=0
      if(cm(2*jo) .ne. 0.0) go to 2
      t=tb(3,ipo)
      ic=0
      go to 3
    2 ko=tb(3,ipo)+.1
      nc=cm(2*jo+1)-cm(2*jo)+1.1
      ic=1
    3 ip=lk
    4 if(ip.le.0) go to 14
      t1=tb(2,ip)
      j1=abs(t1)+.1
      if(j1 .ne. jo) go to 13
      if(ic .ne. 0) go to 6
      if(to*t1 .le. 0.0) go to 13
      if(ieq(t,tb(3,ip),1.0D0) .ne. 1) go to 13
      jg=1
      go to 14
    6 kp=tb(3,ip)+.1
      kg=0
      do 11 l=1,nc
      lo=l+ko
      lp=l+kp
      lon=cm(lo)+.1
      lop=cm(lp)+.1
      if(to .ge. 0.0) go to 8
      if(lon .ne. 0) go to 7
      lon=1
      go to 8
    7 lon=0
    8 if(t1 .ge. 0.0) go to 10
      if(lop .ne. 0) go to 9
      lop=1
      go to 10
    9 lop=0
   10 if(lon .eq. lop) go to 11
      kg=1
      go to 12
   11 continue
   12 if(kg .ne. 0) go to 13
      jg=1
      go to 14
   13 ip=tb(4,ip)+.1
      go to 4
   14 if(jg .ne. 0) go to 15
      lg=1
      go to 16
   15 ipo=tb(4,ipo)+.1
      go to 1
   16 if(lg .ne. 0) go to 17
      ieqbf=1
      go to 18
   17 ieqbf=0
   18 return
      end
      function ibfext(m,tb,cm)
      implicit none
      integer ibfext
      integer m
      double precision tb(5,*),cm(*)
      integer mm1,norm,l
      integer nord,ieqbf
      mm1=m-1
      ibfext=0
      norm=nord(m,tb)
      do 1 l=1,mm1
      if(nord(l,tb).ne.norm) go to 1
      if(ieqbf(l,m,tb,cm) .eq. 0) go to 1
      ibfext=1
      return
    1 continue
      return
      end
      subroutine miss (n,p,x,lx,xm,flg,pn,xn,lxn,xs,xp)
      implicit none
      integer n,p,pn,lx(*),lxn(*)
      double precision flg,x(n,*),xm(*),xn(n,*),xs(*),xp(*)
      double precision s
      integer i,j,mf
      double precision ss
      pn=p
      xp(1)=p
      do 1 j=2,2*p+1
      xp(j)=0.0
    1 continue
      ss=0.0
      do 7 j=1,p
      lxn(j)=lx(j)
      xs(j)=flg
      if(lx(j).eq.0) go to 7
      s=0.d0
      mf=0
      do 4 i=1,n
      xn(i,j)=x(i,j)
      if(x(i,j) .ne. xm(j)) go to 2
      mf=mf+1
      go to 4
    2 if(lx(j) .ge. 0) go to 3
      ss=x(i,j)
      go to 4
    3 s=s+x(i,j)
    4 continue
      if(mf.eq.0) go to 7
      if(mf .ne. n) go to 5
      lxn(j)=0
      go to 7
    5 s=s/(n-mf)
      pn=pn+1
      lxn(pn)=-1
      xs(pn)=1.0
      xp(j+1)=pn
      call nest(n,j,pn,1,1.0D0)
      if(lx(j).gt.0) ss=s
      xp(j+p+1)=ss
      do 6 i=1,n
      xn(i,pn)=1.0
      if(x(i,j) .ne. xm(j)) go to 6
      xn(i,j)=ss
      xn(i,pn)=0.0
    6 continue
    7 continue
      return
      end
      subroutine mkmiss (n,p,x,y,w,xm,pm,nnx,nn,xnn,yn,wn,sc)
      implicit none
      integer nlist
      parameter(nlist=500)
      integer n,p,nnx,nn,m(nlist)
      double precision pm(p),xm(p),x(n,p),y(n),w(n),xnn(*),yn(*),wn(*),s
     1c(p,*)
      integer i,j,jp,km,in,k,nnk
      double precision tol,fin,cvx,r,val
      integer min0
      data tol /0.001/
      if(p .le. nlist) go to 1
c     write(6,'('' increase parameter nlist in subroutine mkmiss to '',i
c    15,                      '' and recompile.'')') p
      stop
    1 do 3 j=1,p
      m(j)=0
      do 2 i=1,n
      if(x(i,j).eq.xm(j)) m(j)=m(j)+1
      sc(j,i)=x(i,j)
      yn(i)=y(i)
      wn(i)=w(i)
    2 continue
    3 continue
      nn=n
      jp=0
      km=jp
    4 if(nn.ge.nnx.or.km.gt.p) go to 13
      jp=jp+1
      if(jp.gt.p) jp=1
      fin=nn*pm(jp)-m(jp)
      if(fin .le. 0.0) go to 5
      in=fin+0.5
      go to 6
    5 in=0
    6 in=min0(in,nnx-nn)
      if(in .le. 0) go to 7
      km=0
      go to 8
    7 km=km+1
      go to 4
    8 do 11 k=1,in
      call rnms(r,1)
      i=nn*r+1.0
      nnk=nn+k
      do 9 j=1,p
      sc(j,nnk)=sc(j,i)
    9 continue
      sc(jp,nnk)=xm(jp)
      yn(nnk)=yn(i)
      wn(nnk)=wn(i)
      do 10 j=1,p
      if(sc(j,nnk).eq.xm(j)) m(j)=m(j)+1
   10 continue
   11 continue
      nn=nn+in
      cvx=-9.9e30
      do 12 j=1,p
      cvx=dmax1(cvx,(nn*pm(j)-m(j))/float(nn))
   12 continue
      if(cvx.lt.tol) go to 13
      go to 4
   13 k=0
      do 15 j=1,p
      do 14 i=1,nn
      k=k+1
      xnn(k)=sc(j,i)
   14 continue
   15 continue
      return
      entry smktol(val)
      tol=val
      return
      end
      subroutine xmiss (n,x,xm,xp,xn)
      implicit none
      integer n
      double precision x(n,*),xm(*),xp(*),xn(n,*)
      integer p,i,j,k
      p=xp(1)+.1
      do 3 j=1,p
      k=xp(j+1)+.1
      do 2 i=1,n
      if(x(i,j) .eq. xm(j)) go to 1
      xn(i,j)=x(i,j)
      if(k.gt.0) xn(i,k)=1.0
      go to 2
    1 xn(i,j)=xp(j+p+1)
      if(k.gt.0) xn(i,k)=0.0
    2 continue
    3 continue
      return
      end
      function nnord(m,tb)
      implicit none
      integer nnord
      integer m
      double precision tb(5,*)
      integer ip,jb
      integer int
      ip=m
      nnord=0
    1 if(ip.le.0) go to 2
      call isnstr(int(abs(tb(2,ip))+.1),jb)
      if(jb.eq.0) nnord=nnord+1
      ip=tb(4,ip)+.1
      go to 1
    2 return
      end
      subroutine spofa(a,m,n,info)

      implicit none

      integer m,n,info

      double precision a(m,*),s,t,u

         integer i,j,j1,jm1,k,km1

         double precision dsqrt

         j1 = info

         do 30 j = j1, n

            info=j

            s = 0.0d0

            jm1 = j - 1

            if (jm1 .lt. 1) go to 20

            do 10 k = 1, jm1

               u=0.0

               km1=k-1

               if(km1.le.0) go to 40

               do 50 i=1,km1

                  u=u+a(i,k)*a(i,j)

   50          continue

   40          continue

               t = a(k,j) - u

               t = t/a(k,k)

               a(k,j) = t

               s = s + t*t

   10       continue

   20       continue

            s = a(j,j) - s

            if (s .le. 0.0d0)  return

            a(j,j) = dsqrt(s)

   30    continue

      info=0

      return

      end

      subroutine sposl(a,m,n,b)

      implicit none

      integer m,n

      double precision a(m,*),b(*),t

      integer i,k,kb,km1

      do 10 k = 1, n

         t = 0.0

         km1=k-1

         if(km1.le.0) go to 30

         do 40 i=1,km1

            t=t+a(i,k)*b(i)

   40    continue

   30    continue

         b(k) = (b(k) - t)/a(k,k)

   10 continue

      do 20 kb=1,n

      k=n+1-kb

      b(k)=b(k)/a(k,k)

      t=-b(k)

      km1=k-1

      if(km1.le.0) go to 50

      if(t.eq.0.0) go to 50

      do 60 i=1,km1

         b(i)=b(i)+t*a(i,k)

   60 continue

   50 continue

   20 continue

      return

      end

      subroutine psort (v,a,ii,jj)

c

c     puts into a the permutation vector which sorts v into

c     increasing order. the array v is not modified.

c     only elements from ii to jj are considered.

c     arrays iu(k) and il(k) permit sorting up to 2**(k+1)-1 elements

c

c     this is a modification of cacm algorithm #347 by r. c. singleton,

c     which is a modified hoare quicksort.

c

      implicit none

      integer ii,jj,a(jj)

      double precision v(jj)

      integer t,tt,l,m,i,j,k,ij,iu(20),il(20)

      double precision vt,vtt

      m=1

      i=ii

      j=jj

 10   if (i.ge.j) go to 80

 20   k=i

      ij=(j+i)/2

      t=a(ij)

      vt=v(t)

      if (v(a(i)).le.vt) go to 30

      a(ij)=a(i)

      a(i)=t

      t=a(ij)

      vt=v(t)

 30   l=j

      if (v(a(j)).ge.vt) go to 50

      a(ij)=a(j)

      a(j)=t

      t=a(ij)

      vt=v(t)

      if (v(a(i)).le.vt) go to 50

      a(ij)=a(i)

      a(i)=t

      t=a(ij)

      vt=v(t)

      go to 50

 40   a(l)=a(k)

      a(k)=tt

 50   l=l-1

      if (v(a(l)).gt.vt) go to 50

      tt=a(l)

      vtt=v(tt)

 60   k=k+1

      if (v(a(k)).lt.vt) go to 60

      if (k.le.l) go to 40

      if (l-i.le.j-k) go to 70

      il(m)=i

      iu(m)=l

      i=k

      m=m+1

      go to 90

 70   il(m)=k

      iu(m)=j

      j=l

      m=m+1

      go to 90

 80   m=m-1

      if (m.eq.0) return

      i=il(m)

      j=iu(m)

 90   if (j-i.gt.10) go to 20

      if (i.eq.ii) go to 10

      i=i-1

 100  i=i+1

      if (i.eq.j) go to 80

      t=a(i+1)

      vt=v(t)

      if (v(a(i)).le.vt) go to 100

      k=i

 110  a(k+1)=a(k)

      k=k-1

      if (vt.lt.v(a(k))) go to 110

      a(k+1)=t

      go to 100

      end

      subroutine stseed (iseed)

      implicit none

      integer iseed

      double precision x(1)

      double precision u

      integer i,j,n

      double precision dmod

      data i /987654321/

      i=iseed

      return

      entry rnms (x,n)

      do 1 j=1,n

      i=dmod(i*16807.d0,2147483647.d0)

      u=i

      u=u*.465661287d-9

      x(j)=u

    1 continue

      return

      end

