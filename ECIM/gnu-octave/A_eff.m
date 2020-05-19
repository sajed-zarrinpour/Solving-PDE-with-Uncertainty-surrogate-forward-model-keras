function A = A_eff(u,u_T,a,L,b,h,d)
 A = h^d * (u_T * L * u - 2* u_T * b + a.' * ones(length(a),1));
end 