#' Computes probability density function.
#'
#' @param t A numeric. The failure time at which the density should be computed.
#' @inheritParams ESPA_supp
#'
#' @return A numeric. The probability density estiamtion at t.
#' @export

ESPA_pdf=function(t,time, status){

    Y=log(time+1)


    delta=status[sort.list(Y)]
    z=sort(Y)
    n=length(z)

    pw=rep(0,n)
    pw[1]=delta[1]/n
    prod.term=1
    for (i in 2:n){
        prod.term = prod.term * (1 - delta[i-1]/(n-i+2))
        pw[i] = prod.term * delta[i]/(n-i+1)
    }


    p.str=sum(pw)

    phi=0
    if (delta[n]==0) {
      phi = -log(1-p.str)/z[n]
      }




    M.mgf=function(s){
        Mds=sum(pw*exp(s*z))
        Mcs=0
        if (delta[n]==0) {
          Mcs=phi*exp(z[n]*(s-phi))/(phi-s)
          }
        Ms=Mds+Mcs
        return(Ms)
    }



    K.cgf=function(s){
        Mds=sum(pw*exp(s*z))
        Mcs=0
        if (delta[n]==0) {
          Mcs=phi*exp(z[n]*(s-phi))/(phi-s)
          }
        log(Mds+Mcs)
    }



    K2.exp=function(s){
        Mds=sum(pw*exp(s*z))
        Mcs=0
        if (delta[n]==0) {
          Mcs=phi*exp(z[n]*(s-phi))/(phi-s)
          }
        Ms=Mds+Mcs

        Mds.prime=sum(pw*z*exp(s*z))
        Mcs.prime=0
        if (delta[n]==0) {
          Mcs.prime=Mcs*(z[n]+1/(phi-s))
          }
        Ms.prime=Mds.prime+Mcs.prime

        Mds.prime2=sum(pw*z^2*exp(s*z))
        Mcs.prime2=0
        if (delta[n]==0) {
          Mcs.prime2=Mcs.prime*(z[n]+1/(phi-s))+Mcs/(phi-s)^2
          }
        Ms.prime2=Mds.prime2+Mcs.prime2

        K.prime2=(Ms*Ms.prime2-Ms.prime^2)/Ms^2
        K.prime2
    }

    fh.com=function(t,s){
        K.exp=K.cgf(s)
        pdf.da=exp(K.exp-t*s)/sqrt(2*pi*K2.exp(s))
        return(pdf.da)
    }


    K.adj.exp=function(s){
        Mds=sum(pw*exp(s*z))
        Mcs=0

        if (delta[n]==0){
            Mcs=phi*exp(z[n]*(s-phi))/(phi-s)
        }
        Ms=Mds+Mcs

        if(round(Ms,3)==0|round(K2.exp(s),3)<=0){
            return(10000)
        }
        else{
            return(K.cgf(s)-log(t+1)*s)
        }
    }


    saddle.com<-function(){
        suppressWarnings(optim(0, K.adj.exp, method="Nelder-Mead",hessian=FALSE)$par)
    }


    sh=saddle.com()

    return((fh.com(log(t+1),sh))/(t+1))

}
