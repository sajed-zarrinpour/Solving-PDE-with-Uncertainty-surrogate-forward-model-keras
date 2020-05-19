#https://math.unice.fr/~frapetti/CorsoF/cours4part2.pdf

function [sigma,x,iter,relres]=invpower(A,z0,mu,tol,nmax)
  %INVPOWER Inverse power method
  % [SIGMA,X,ITER,RELRES]=INVPOWER(A,Z0,MU,TOL,NMAX) computes the
  % eigenvalue LAMBDA of smallest module of the matrix A and the
  % corresponding eigenvector X of unit norm. TOL specifies the tolerance of the
  % method. NMAX specifies the maximum number of iterations. X0 specifies
  % the initial guess. MU is the shift. ITER is the iteration number at which
  % X is computed.
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
      fprintf('Multiple eigenvalue'); break;
    end
    x=q;
  end

return