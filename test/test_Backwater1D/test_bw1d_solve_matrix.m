% Sun 10 Nov 17:48:39 +08 2019

bw = Backwater1D();

Xi = [0,5e5];
x = linspace(Xi(1),Xi(end))';
Q0 = 1e4;
w = 500;
zb = -15*(1 - x/3e5);
zs = [];

x0 = x(1);
z0 = 0;

bw.solve_matrix(x,zs,Q0,Qt,Qmid,Qhr,chezy,width,zb,x0,z0)

