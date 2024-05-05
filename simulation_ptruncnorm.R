#Simulation pnormtrunc
library(truncnorm)
library(rlist)
          
ptruncnorm(q=1.58, a=1.424, b=6.29, mean=0, sd=0.17)
ptruncnorm(q=1.58, a=1.424, b=6.29, mean=0, sd=0.2)


test_seq=seq(from=0.17, to=0.2, by=0.0001)
test_seq=rev(test_seq)

a<-ptruncnorm(q=1.58, a=1.424, b=6.29, mean=0, sd=test_seq)

plot(a)

diffs<-rep(0,length(a)-1)

for (i in 2:length(test_seq)){
  diffs[i]=a[i]-a[i-1]
}


diffs

#We can see that since some diffs are >0 and some are <0 that the  relationship (at least computationally) is not monotonours.
#I don't think that this matters for our purposes though.
