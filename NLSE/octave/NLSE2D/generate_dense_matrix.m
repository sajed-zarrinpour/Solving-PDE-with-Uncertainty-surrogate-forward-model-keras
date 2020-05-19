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
## @deftypefn {} {@var{A} =} generate_dense_matrix (@var{p}, @var{q}, @var{h}, @var{a})
##
## @seealso{}
## @end deftypefn

## Author: samim <sa.zarrinpour@iasbs.ac.ir>
## Created: 2020-05-05

function A = generate_dense_matrix (p, q, h, a)
    # this line is neccessary since I dont want matlab to confuse with the type pf this two variable.  comment it if you want to know what will happen. 
    i=j=1; 
    %Declaring the problem
    A = zeros(p^2,q^2);
    for r = 1:p*q
      #for each row
      if mod(r,q)==0
        j = q;
      else
        j = mod(r,q);
      end
      i = floor((r-1)/p) + 1;
      if i>p
        i=mod(i,p);
      end
      if j>q
        j = mod(j,q);
       end
      # i+1, j
     if i == p 
      A(r, j) = -1;
     else
      A(r,i*p + j) = -1;
     end
     # i-1, j
     if i==1
      A(r, (p-1)*p + j) = -1;
     else
      A(r, (i-2)*p +j) = -1;
     end
     # i, j+1
     if j==q
      A(r, (i-1)*p + 1) = -1;
     else
      A(r, (i-1)*p + (j+1)) = -1;
     end
     # i, j-1
     if j==1
      A(r, (i-1)*p + q) = -1;
     else
      A(r, (i-1)*p + (j-1)) = -1; 
     end
     A(r, r) = h^2 * a(i,j) + 4;
    end
endfunction
