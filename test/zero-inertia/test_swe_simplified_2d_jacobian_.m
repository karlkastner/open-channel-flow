% 2024-09-27 11:20:43.942523355 +0200
K  = 1;
dt = 1;
L  = 4*[1,1];
n  = 4*[1,1];
h = ones(n(1)*n(2),1);

J = swe_simplified_jacobian(h,K,dt,L,n);

[Dx,Dy,D2x,Dxy,D2y] = derivative_matrix_2d(n,L-1,2,{'circular','circular'});
J_ = speye(prod(n)) - K*dt*(D2x+D2y);

full(J)
full(J_)
rms(J-J_)
