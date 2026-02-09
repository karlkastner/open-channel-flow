% 2020-09-03 09:14:27.683951002 +0200
syms h w g A t S cd positive;
 syms Q(t);
 s = dsolve(diff(Q,t) + g*(h*w)*S - cd*w*Q^2/(h*w)^2,Q(0)==0);
 simplify(s);
 f=matlabFunction(s);
 t=linspace(0,86400/10);
 plot(t,f(5e-5,2.5e-3,9.81,10,t,1))

