clear;clc;
#number of grid points along x1
p = 9;
#number of grid points alon x2
q = 9;
#lower/upper bounds of uncertainty with uniform distribution
lambda0 = 1;
lambda1 = 16;
h=0.125;
number_of_samples = 6750;
I = eye(p*q,p*q);
test_percentage = 0.15;
test_indexes = randperm(number_of_samples,ceil(test_percentage*number_of_samples));
k = 1;
while true
   p = q = 9;
    # initializing the uncertainty over the domain.
    a = (lambda1 - lambda0).*rand(p, q) + lambda0;
    #generating the dense matrix 
    OA = generate_dense_matrix(p, q, h, a);

    %getting the smallest eigen values with its corresponding eigen vector
    [V, E] = eigs(OA,1,'sm');
    #V = V*h;
    #fprintf('Begining with E0: %d\n',E);

    for s=0:0.4:2
        A =  update_A(OA, V, s,h);
        %[E,V,iter,relres] = invpower(CA, V, E, 1e-5, 1e4);
          z0=V;
          mu=E;
          tol  = 1e-5;
          nmax = 1e4;
          
          M=A-mu*eye(size(A)); [L,U,P]=lu(M);
          q=z0/norm(z0); q2=q'; sigma=[ ];
          relres=tol+1; iter=0;
          while relres(end)>=tol && iter<=nmax
              iter=iter+1;
              b=P*q;
              y=L\b; z=U\y;
              q=z/norm(z); z=A*q; sigma=q'*z;
              b=q2'; y=U'\b; w=L'\y;
              q2=w'*P; q2=q2/norm(q2); costheta=abs(q2*q);
            if costheta>=5e-2
              temp=norm(z-sigma*q)/costheta; relres=[relres,temp];
            else
              fprintf('Multiple eigenvalue');break;
            end
            x=q;
          end
          E = min(sigma);
          V = x;
        #V = V*h;
        #fprintf('sigma = %d,  E0= %d\n',s, E);
    end
    E = E/h^2;
    fprintf('Iteration (%d) ended with  E0: %d\n',k, E);
    #S.(sprintf('%d',k)) = [V',E]
    #save -6 data.mat -append -struct S
    if E>9.7 && E<10.99
    if ismember(k, test_indexes)
      save(sprintf('data/test/data%d.mat',k),'-6', '-append', "[V',E, a]");
    else
      save(sprintf('data/train/data%d.mat',k),'-6', '-append', "[V',E, a]");
    end
    if k == number_of_samples
      break;
    else
      k = k+1;
    end
    end
    
end
