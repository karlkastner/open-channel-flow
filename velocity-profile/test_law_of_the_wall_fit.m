% 2017-01-04 00:48:07.865083573 +0100

function test

p = 6;

kappa=0.4;
us = 0.1;
z0 = 0.1;
H = 10;
n = 100;
z = linspace(z0,H,n)';

% cubic function, top value and derivative, bottom value and derivative (mean)
u = us/kappa*log(z/z0);
% TODO analytic
U = mean(u);
U
for idx=0:p
	A = vander_1d(z,idx);
	c = A \ u;
	u_(:,idx+1) = A*c;
end
clf
subplot(2,2,1)
plot(z,[u u_])
%view([90 -90])

% asymptotically at zt (not mean preserving, diverges)
% only valid within zt<z<2zt
 zt = 0.5*H;
 %zt = H;
% zt = exp(-1)*H;
Ar  = [];
rhs = [];
u_  = [];
C   = [];
A   = vander_1d((z-zt),p);
D = diag(1./factorial(0:p))
for idx=0:p
	Ar(:,idx+1) = vanderd_1d(z0,idx,3)
	% regress
	if (0==idx)
		rhs = us/kappa*log(zt/z0);
	else
		rhs(idx+1,1) = -(-1)^idx*factorial(idx-1)*us/kappa*1/zt^idx;
		%rhs(idx+1,1) = -(-1)^idx*factorial(idx-1)*us/kappa*1/z0^idx;
	end
	c = D(1:idx+1,1:idx+1)*rhs;
	%c = (Ar(1:idx+1,:)) \ rhs;
	%c = (D*Ar(1:idx+1,:)) \ rhs;
	C(1:idx+1,idx+1) = c;
end
u_ = A*C;
%u_(:,idx+1) = A(:,1:idx+1)*c;
C
subplot(2,2,2)
plot(z,[u u_])
%view([90 -90])
ylim([0 2])

% double asymptote: match n derivatives at z0 and at surface

