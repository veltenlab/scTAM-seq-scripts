data {
  //the actual data
  int n0; #number of zeros
  int n1; #number of ones
  real<lower=0,upper=1> p; #observed dropout rate
}

parameters {
  real<lower=0,upper=1> m;
}


model {
  n1 ~ binomial(n0 + n1, (1-p) * m);
  m ~ beta(1,1); # flat prior
}
