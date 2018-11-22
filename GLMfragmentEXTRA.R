load("~/Downloads/rstuff/Normal_GLM.Rdata")
xsq = x^2
fit = lm(y~x+xsq)
sigsq = sum(fit$residuals^2)/47
b_1 = fit$coefficients/sigsq
b_2 = c(1/sigsq,0,0)
b=c(b_1,b_2)
lmlik = lik(b_1,b_2,y,model.matrix(fit),model.matrix(fit))

D2l=function(b)  {D2lik(b[1:3],b[4:6],y,model.matrix(fit),model.matrix(fit))}
Dl=function(b)  {Dlik(b[1:3],b[4:6],y,model.matrix(fit),model.matrix(fit))}
l=function(b)  {lik(b[1:3],b[4:6],y,model.matrix(fit),model.matrix(fit))}
h_m=function(b){-solve(D2l(b))%*%Dl(b)}


#opt = optim(b,l,gr=Dl,control= list(fnscale = -1))


#Newton-Raphson! This now works. It produces a better minimum than 'optim.' Probably because this
#is a second-order method.
b=c(b_1,b_2)
h=h_m(b)
del=1
while(del>10^-35){ # 10^-35 is around 100 times the size of single-precision 'machine zero.'
  #In reality this will be satisfied when l(b) and l(b+h) agree to ~16 digits,
  #which will happen much much sooner.
while(is.nan(l(b+h))||l(b+h)<l(b)){
  h = h/2
}
del = l(b+h)-l(b)
b = b+h
h = h_m(b)
}
bH2 = b

#Load 'b' with the \beta coefficients!! This computes each of your 'eta's
eta1 = model.matrix(fit)%*%b[1:3] #\eta_{1,i} = X_i*\beta_1
eta2 = model.matrix(fit)%*%b[4:6] #\eta_{2,i} = X_i*\beta_2

#plots plots plots!!

#Plot of 'y' values, fitted values from the GLM, and fitted values from the linear model
#Black is the 'y' values, red is the linear model and green is the GLM
plot(y)
points(fit$fitted.values,col="red")
points(eta1/eta2, col="green")

#Residual plots
#Residuals from the GLM
plot(y - eta1/eta2)
lines(sqrt(1/eta2))
lines(-sqrt(1/eta2))
#Residuals from the linear model
plot(fit$residuals)
points(sqrt(sigsq)*replicate(50,1),pch="-")
points(sqrt(sigsq)*replicate(50,-1),pch="-")
#It looks good, since there are lots of points within the +/- 1 stdev bars...
#But actually there are TOO MANY points within the bars! We should only have ~70% between the bars!
#~30% of the points SHOULD be outside the bars, but only 7 of 50 (14%) points are outside the bars.



#Fitting the GLM with x^2 covariant for eta_1 equal to zero.
#
#WARNING: The 'h' step we're using in Newton-Raphson is NOT the gradient of the function.
#This means we can't just delete the associated entry in the gradient vector to hold one
#variable constant! We need to use the Hessian without that variable, AND the gradient without
#that variable, so that's what's going on in the new h_m function.
#Since the covariant we want is the third in the unrolled list, we can just take D2l(b)
#without the third column and row, and Dl(b) without its third entry
#We can't just set those columns and rows equal to zero since we need to invert the matrix!
h_mp1=function(b){
  H = D2l(b)
  H2=solve(H[-3,-3])
  G = Dl(b)[-3]
  h_m = -H2%*%G
  c(h_m[1:2],0,h_m[3:5])
  }
fit2 = lm(y~x)
b=c(fit2$coefficients,0,b_2)

h=h_mp1(b)
del=1
while(del>10^-35){ #
  while(is.nan(l(b+h))||l(b+h)<l(b)){
    h = h/2
  }
  del = l(b+h)-l(b)
  b = b+h
  h = h_mp1(b)
  print(del)
}
bH0 = b
#Boom. bH0 is now our vector of betas for the hypothesis above.


#We fit the model for the hypothesis, the beta associated to x^2 in eta_2 is zero
#Like before, we have to delete the variable out of the Hessian AND the gradient
#to compute 'h_m.' To do this we take the Hessian and delete its sixth column and row,
#and delete the sixth entry from the gradient. 
h_mp2=function(b){
  H = D2l(b)
  H2=solve(H[-6,-6])
  G = Dl(b)[-6]
  h_m = -H2%*%G
  c(h_m,0)
}
b=c(b_1,b_2)
h=h_mp2(b)
del=1
while(del>10^-35){
  while(is.nan(l(b+h))||l(b+h)<l(b)){
    h = h/2
  }
  del = l(b+h)-l(b)
  b = b+h
  h = h_mp2(b)
  
}
bH1 = b

#The big one; part e!!
yt = (c(1.0,0.75,0.75^2)%*%bH2[1:3])/(c(1.0,0.75,0.75^2)%*%bH2[4:6])
N=100
ysims1 = vector("numeric",N)
sigsims = vector("numeric",N)
ysims2 = vector("numeric",N)
for(i in 1:N){
  #First, simulate data:
  #N(m,s) = m + \sqrt(s)N(0,1)
  n  = rnorm(50)
  sim = eta1/eta2 + sqrt(1/eta2)*n
  #Set up the gradient functions so we can fit the GLM
  D2ls=function(b)  {D2lik(b[1:3],b[4:6],sim,model.matrix(fit),model.matrix(fit))}
  Dls=function(b)  {Dlik(b[1:3],b[4:6],sim,model.matrix(fit),model.matrix(fit))}
  ls=function(b)  {lik(b[1:3],b[4:6],sim,model.matrix(fit),model.matrix(fit))}
  h_ms=function(b){-solve(D2ls(b))%*%Dls(b)}
  #Intialize values for fitting
  b=c(b_1,b_2)
  h=h_ms(b)
  del=1
  #Fit with Newton-Raphson
  while(del>10^-35){
    while(is.nan(ls(b+h))||ls(b+h)<ls(b)){
      h = h/2
    }
    del = ls(b+h)-ls(b)
    b = b+h
    h = h_ms(b)
  }
  bH0=b
  #Collect the just-fitted model's prediction at x=0.75
  ysims1[i] = (c(1.0,0.75,0.75^2)%*%bH0[1:3])/(c(1.0,0.75,0.75^2)%*%bH0[4:6])
  sigsims[i]=(c(1.0,0.75,0.75^2))%*%(1/bH0[4:6])
  #and the standard linear model's prediction, too
  ysims2[i] = c(1.0,0.75,0.75^2)%*%lm(sim~x+xsq)$coefficients
  print(i)
}
mean(ysims1-yt)
sd(ysims1-yt)
#~0.016 and ~0.122
mean(ysims2-yt)
sd(ysims2-yt)
#~0.695 and ~0.127
#So the linear model is more biased, and has basically the same standard deviation. Terrible.

#Simulation to find an empirical p-value;
#This function runs ONE simulation from the hypothesis H

llsim=function(H){
  eta1H = model.matrix(fit)%*%H[1:3]
  eta2H = model.matrix(fit)%*%H[4:6]
  sims = eta1H/eta2H + sqrt(1/eta2H)*rnorm(50)
  D2ls=function(b)  {D2lik(b[1:3],b[4:6],sims,model.matrix(fit),model.matrix(fit))}
  Dls=function(b)  {Dlik(b[1:3],b[4:6],sims,model.matrix(fit),model.matrix(fit))}
  ls=function(b)  {lik(b[1:3],b[4:6],sims,model.matrix(fit),model.matrix(fit))}
  h_ms=function(b){-solve(D2ls(b))%*%Dls(b)}
  b=H
  h = h_ms(b)
  del=1
  while(del>10^-35){
    while(is.nan(ls(b+h))||ls(b+h)<ls(b)){
      h = h/2
    }
    del = ls(b+h)-ls(b)
    b = b+h
    h = h_ms(b)
  }
  2*abs(ls(b)-ls(H))
}
# This takes a LONG time to run! Uncomment at your own risk. Or try with a smaller number.
# B = replicate(1000,llsim(b))
# plot.ecdf(B)
# x = (0:10000)/100
# lines(x, pchisq(x,5), col='red')
# BBH0 = replicate(10000,llsim(bH0))
# BBH1 = replicate(10000,llsim(bH1))
# BBlin = replicate(10000, llsim(b))

