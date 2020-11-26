

pr2eir=function(x, b=0.55, r=1/200, k=4.2){
  ((1-x)^-k-1)*(r/k/b)
}

pr2k1 = function(x, mx=1.8, s1=4, pw=1.6){
  mx + s1*(1-x)^pw
}

pr2k2 = function(x, sl = 1.5, a=6.2, b=0.2){
  sl +  (a-sl)*(1-x)/(b+(1-x))
}

pr2k = function(x, fac=2.6){
  (pr2k1(x) + pr2k2(x))/fac
}

pr2kS = function(x, s1=6, s2 = 2, n=1, short=TRUE){
  pr = rep(x, each=n)
  mn = pr2k(pr)
  vr = .03 + s2*(1-pr)^1.3
  kk = pmax(1, rnorm(length(pr), mn, vr))
  if(short==TRUE) return(kk) else return(cbind(pr, kk))
}

x2kappa = function(x, k1=.3, k2=0, pw=1){
  k1*x^pw *exp(-k2*x)
}

pr2rc = function(x, b=0.55, r=1/200, s=1, c = 0.175, k=NULL, k2=0, pw=1){
  D = c/r
  if(is.null(k)) k = pr2k(x)
  eir = pr2eir(x, b, r, k)
  kappa = x2kappa(x, c, k2, pw)
  vc = eir*(1+s*kappa)/kappa
  Rc = b*D*vc
  list(pr=x, k=k, aeir=eir*365, kappa=kappa, vc=vc, Rc=Rc)
}

pr2ar = function(x, y){x}

ar2pr = function(x, y){x}

pr2rcS = function(x, rho=0, b=0.55, r=1/200, s=1, c = 0.175, k=NULL, k2=0, pw=1, n=1, Short=TRUE){
  pr = rep(x, each=n)
  alpha = pr2ar(pr, rho)
  pr1 = ar2pr(alpha, 0)
  L = length(pr)
  rr = 1/rnorm(L,1/r,5)
  cc = rbeta(L, c*50, (1-c)*50)
  D = cc/rr
  k = pr2kS(pr)
  bb = rbeta(L, b*100, (1-b)*100)
  eir = pr2eir(pr, bb, rr, k)
  kappa = x2kappa(pr, k1=cc)
  ss=s
  vc = eir*(1+ss*kappa)/kappa
  Rc = bb*D*vc
  if(Short==TRUE){return(Rc)} else{
    list(pr=pr, r=rr, c=cc, D=D, k=k, bb=bb, aeir=eir*365, kappa=kappa, s=ss, Rc=Rc)
  }
}
