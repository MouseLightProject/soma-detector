sx =  1 
sy =  2 
rho = 0
Sigma = [ sx*sx rho*sx*sy ; ...
          rho*sx*sy sy*sy ] 
k = 3
Sigma_b = [ sx*sx rho*sx*k*sy ; ...
            rho*sx*k*sy k*k*sy*sy ] 
        
[S, Lambda] = eig(Sigma)      

[S_b, Lambda_b] = eig(Sigma_b)
