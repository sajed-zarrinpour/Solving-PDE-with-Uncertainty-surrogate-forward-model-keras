## Copyright (C) 2020 samim
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} powerrm (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: samim <samim@samim-pc>
## Created: 2020-05-06
function [lambda,x,iter,relres]=powerm(A,z0,tol,nmax)
  %POWERM Power method
  % [LAMBDA,X,ITER,RELRES]=POWERM(A,Z0,TOL,NMAX) computes the
  % eigenvalue LAMBDA of largest module of the matrix A and the corresponding
  % eigenvector X of unit norm. TOL specifies the tolerance of the method.
  % NMAX specifies the maximum number of iterations. Z0 specifies the initial
  % guess. ITER is the iteration number at which X is computed.
  q=z0/norm(z0); q2=q;
  relres=tol+1; iter=0; z=A*q;
  while relres(end)>=tol & iter<=nmax
      h=0.125;
      q=z/norm(z); z=A*q;
      lambda=q'*z; x=q;
      z2=q2'*A; q2=z2/norm(z2); q2=q2';
      y1=q2; costheta=abs(y1'*x);
    if costheta >= 5e-2
      iter=iter+1;
      temp=norm(z-lambda*q)/costheta;
      relres=[relres; temp];
    else
      fprintf('Multiple eigenvalue'); break;
    end
  end
return
