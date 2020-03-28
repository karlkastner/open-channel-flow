% 2016-05-17 23:00:06.927153039 +0200
% Karl Kastner, Berlin

function test_inverse_backwater_curve()

Q   = 6.5e3;
%W   = 500;
W   = @(x) 300 + 300*x/3e5;
C   = 50;
%S0 = @(x) 1e-5 + (1e-4-1e-5)*x/3e5;
%S0  = @(x) 1e-5*(x < 1.5e5) + 1e-4*(x >= 1.5e5);
S0 = @(x) 1e-4;
X   = [0 3e5];

h0  = normal_flow_depth(Q,W(0),C,S0(0));

% forward backwater curve
[x h] = backwater_curve(Q,C,W,S0,h0,X);

% surface slope
yb = + cumsum(S0(x).*cdiff(x));
ys = h + yb;


% test
% surface level and surface slope interpolator (second order)
[ys1 dys_] = interp1_slope(x,ys,x);
norm(ys-cvec(ys1))
ys2 = cumsum(cvec(dys_).*cdiff(x)) + ys(1);
norm(ys-cvec(ys2))
figure(1)
clf
plot(x,[yb ys ys1' ys2])

% inverse backwater curve
yb0 = yb(1) + 5;
[xi ybi] = inverse_backwater_curve(Q,C,W,@(x0) interp1_slope(x,ys,x0),yb0,X);

figure(2)
clf();
plot(x,yb);
hold on
plot(xi,ybi,'.');
plot(x,ys,'k-');

ybi_ = interp1(xi,ybi,x,'linear');
norm(yb-ybi_)./norm(yb)

end

