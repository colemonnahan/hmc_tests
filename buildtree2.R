buildtree2 <- function(theta, r, u, v, j, eps, theta0, r0, fn, gr,
                       delta.max=1000, info = environment() ){
  if(j==0){
    ## base case, take one step in direction v
    eps <- v*eps
    r <- r+(eps/2)*gr(theta)
    theta <- theta+eps*r
    r <- r+(eps/2)*gr(theta)
    ## verify valid trajectory
    H <- .calculate.H(theta=theta, r=r, fn=fn)
    s <- H-log(u) + delta.max > 0
    if(is.na(s) | is.nan(s)) s <- 0
    n <- log(u) <= H
    ## Useful code for debugging. Returns entire path to global env. as
    ## long as s=1
    if(s==1){ if(!exists('theta.trajectory'))
      theta.trajectory <<- theta
    else
      theta.trajectory <<- rbind(theta.trajectory, theta)
      }
    temp <- .calculate.H(theta=theta, r=r, fn=fn)-
      .calculate.H(theta=theta0, r=r0, fn=fn)
    alpha <- min(exp(temp),1)
    info$n.calls <- info$n.calls + 5
    return(list(theta.minus=theta, theta.plus=theta, theta.prime=theta, r.minus=r,
                r.plus=r, s=s, n=n, alpha=alpha, nalpha=1))
  } else {
    ## recursion - build left and right subtrees
    xx <- buildtree2(theta=theta, r=r, u=u, v=v, j=j-1, eps=eps,
                     theta0=theta0, r0=r0, fn=fn, gr=gr, info=info)
    theta.minus <- xx$theta.minus
    theta.plus <- xx$theta.plus
    theta.prime <- xx$theta.prime
    r.minus <- xx$r.minus
    r.plus <- xx$r.plus
    alpha <- xx$alpha
    nalpha <- xx$nalpha
    s <- xx$s
    if(is.na(s) | is.nan(s)) s <- 0
    nprime <- xx$n
    ## If it didn't fail, update the above quantities
    if(s==1){
      if(v== -1){
        yy <- buildtree2(theta=theta.minus, r=r.minus, u=u, v=v,
                         j=j-1, eps=eps, theta0=theta0, r0=r0,
                         fn=fn, gr=gr, info=info)
        theta.minus <- yy$theta.minus
        r.minus <- yy$r.minus
      } else {
        yy <- buildtree2(theta=theta.plus, r=r.plus, u=u, v=v,
                         j=j-1, eps=eps, theta0=theta0, r0=r0,
                         fn=fn, gr=gr, info=info)
        theta.plus <- yy$theta.plus
        r.plus <- yy$r.plus
      }
      ## This isn't in the paper but if both slice variables failed,
      ## then you get 0/0. So I skip this test. Likewise if model
      ## throwing errors, don't keep that theta.
      nprime <- yy$n+ xx$n
      if(!is.finite(nprime)) nprime <- 0
      if(nprime!=0){
        ## choose whether to keep this theta
        if(runif(n=1, min=0, max=1) <= yy$n/nprime)
          theta.prime <- yy$theta.prime
      }
      ## check for valid proposal
      test <- .test.nuts(theta.plus=theta.plus,
                         theta.minus=theta.minus, r.plus=r.plus,
                         r.minus=r.minus)
      ## if(!test) warning(paste("U turn at j=", j))
      ## check if any of the stopping conditions were met
      s <- xx$s*yy$s*test
    }
    return(list(theta.minus=theta.minus, theta.plus=theta.plus,
                theta.prime=theta.prime,
                r.minus=r.minus, r.plus=r.plus, s=s, n=nprime,
                alpha=alpha, nalpha=1))
  }
}
