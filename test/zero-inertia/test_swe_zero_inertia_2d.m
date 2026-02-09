% Mon  7 Oct 21:37:40 CEST 2024
% TODO linear extrap for z at boundary instead of circular to have circular in dz/dx not z
% Timing bicgstabl-ilu, euler_fw (midpoint might already improve)
% n    dt   precision  runtime	iter-gn	iter-bicg
% 2^10 0.25 double     13.8	4	2.5
% 
if (0)
	n = 2^4*[1,1]+1;
	L = n;
	h0 = 1;
	C  = 10;
else
	L  = 2^8*[1,1]+1;
	n  = L;
	s  = 10*2*L./n;
	% 10 mm water depth
	h0max  = 1e-2;
	h0min  = 1e-2*h0max;
	zrange = 0.*[0;1];
	zrange = [1;0];
	%zrange = [1;1];
	zrange = zrange/(sqrt(eps)+norm(zrange));
	% K = 10/d
	C      = 1.33e3;
	dz_dx  = [-1e-2,0,1e-2];
end
	%method = 'iterative';
	%method = 'mg';
	
	mode_C = {'central'}
	mode_C = {'upwind'};
%	mode_C = {'optimal'};
%	method = 'optimal';
%	method = 'central';
	nvar = 1;
maxiter = 1e2;
dt = 1
dx = L./n;

x = ((0:n(1)-1)'+1/2)*L(1)/n(1);
y = x';
%zx = zrange/2*(1 + cos(pi*x/L(1)));
zx = zrange(1)/2*(1 + sin(2*pi*x/L(1)));
%zx = zrange(1)*(1-x/L(1));
zy = zrange(2)*(1-x/L(1));
%zy = ones(1,n(2));
z = cvec(zx)+rvec(zy);
%z = 1/(2*pi)*(1 + sin(2*pi*x/L(1)))*ones(1,n(2));
%z = z';
z  = z+0*1e-2*randn(n);
%z  = [];
zC = {z};
%z = z+z';
%z = -0*x/L(1)*ones(1,n(1));

mu = L/2+sqrt(eps);
% slow convergence when offset is near zero
hx = normpdf(x,mu(1),s(1))/normpdf(0,0,s(1));
hy = normpdf(y,mu(2),s(2))/normpdf(0,0,s(2));%/normpdf(0,0,s)^2;
%hy = ones(1,n(2));
%hx = cvec(tukeywin(n(1)))>0.99;
%hx = 0.5*(1+sin(2*pi*x/L(1)+0.75));
%hy = hx';
%hy = rvec(tukeywin(n(2)))>0.99;
%hx = 2*rand(n(1),1).^0;
%hy = 2*rand(1,n(2));
if (0)
hx(1) = 1;
hx(end)=1;
hy(1) = 1;
hy(end)=1;
for idx=1:0
hx = meanfilt1(hx,3);
hy = meanfilt1(hy,3);
end
end
%hy = ones(1,n(2));
h0 = cvec(hx)*rvec(hy);
%h0 = h0.^2;
h0 = h0/max(h0,[],'all');
h0 = h0min+h0max*h0;
h0 = h0';
%h0 = 0+10*(cvec(normpdf(x,L(1)/2,s))*ones(1,n(1)))';
h=h0;
%dz_dt1 = reshape(swe_zero_inertia_dzdt(flat(h),flat(z),K,L,n,true),n);
%dz_dt2 = reshape(swe_zero_inertia_dzdt(flat(h),flat(z'),K,L,n,true),n);
%z
%h
%dz_dt1 
%dz_dt2t = dz_dt2'
%rms(dz_dt1-dz_dt2','all')
%rms(dz_dt2 - fliplr(dz_dt2),'all')

% 1 - dt*d*4/dx^2 - dt*d*2/dy^2
% dtmax = dx^2/(4*d)

dt_max = min((dx(1)^2+dx(2).^2)./(2*C*h0max))

gnsolver = 'ilu-ks';
gnsolver = 'direct';
for mdx=1:length(mode_C)
method = mode_C{mdx};

for jdx=2
	rmse = [];

tic
switch (jdx)
case{1}
	% explicit solver
	nt = 40*ceil(dt/dt_max);
	dt_ = dt/nt;
	z = flat(z);
	h0 = flat(h0);
	h = h0;
	for idx=1:nt
		% explict midpoint
		dh_dt = swe_zero_inertia_2d_dhdt(h,z,C,L,n,method,true);
		[idx, min(dt_*(dh_dt))]
		h = h + 1.0*dt_*dh_dt;
		%dz_dt = swe_zero_inertia_2d_dhdt(h_,z,C,L,n,method,true);
		%h = h + dt_*dz_dt;
		%imagesc(log(reshape(h,n)))
		%pause(1)
		if (min(h,[],'all')<-1e-3)
			idx
			error('dt too large')
		end
	end
case {2} % krylov-subspace
	z  = flat(z);
	h0 = flat(h0);
	h  = h0;
	returnmat         = true;
	tic
%	J = jfun(h0);
	resfun            = @(h) euler_implicit_residual(@(h) swe_zero_inertia_2d_dhdt(h,z,C,L,n,method,returnmat), dt, h, hold);
	%res = resfun(h0);
	%rms(res)
	[dhdt,bl,bu] = swe_zero_inertia_2d_dhdt(h0,z,C,L,n,method,returnmat);
	[J,out] = swe_zero_inertia_2d_jacobian(h0,z,C,L,n,method,returnmat);
	rms(out.vu,'all')
	rms(out.vl,'all')
	rms(J(:))
if (0)
	
	sum(~isfinite(J(:)))
	sum(abs(J(:))>0)
pause
end
	%nt_ = round(nt/5);
	nt_ = 1;
	hold = h0;
	for idx=1:nt_
	idx
	resfun            = @(h) euler_implicit_residual(@(h) swe_zero_inertia_2d_dhdt(h,z,C,L,n,method,returnmat), dt/nt_, h, hold);
	jfun              = @(h) euler_implicit_jacobian(@(h) swe_zero_inertia_2d_jacobian(h,z,C,L,n,method,returnmat), dt/nt_, h, returnmat);
	[h,rmse,res,iter] = gauss_newton(resfun,h,[],maxiter,gnsolver,jfun);
	h = real(h);
	h = max(h,0);
	hold = h;
	end
	figure(10+jdx)
	semilogy(rmse,'.')
	toc
if (0)
	rms(flat(out.J-out.J'))
	D = abs(diag(out.J));
	max((sum(abs(out.J),2)-D)./D)
	eigs(out.J,1,'largestreal')
	eigs(@(x) bicgstabl(out.J,x,100),prod(n),2,'largestabs')
	figure(10)
	clf
	subplot(2,2,1)
		abssumcol  = abs(sum(J,1)-1);
		sumabscol = sum(abs(J),1)-1;
		maxrelcolsum = max(abssumcol./sumabscol);
	s = abs(sum(out.J,1)-1); %./(sum(abs(out.J),1)-1);
	plot(s)
	s= reshape(s,n);
	subplot(2,2,2)
	imagesc(log10(s))
	colorbar
	J = out.J;
	%whos J
	subplot(2,2,3)
	A = out.J - speye(size(out.J));
	fdx = abs(A) > sqrt(eps);
	hist(abs(A(fdx)),n(1))
	%max(flat(abs(J-J')))
	%spy(J-J')
	%m = 24;
	%full(J(1:m,1:m)-eye(m))
	%Jt = J';
	%full(J(1:m,1:m)-Jt(1:m,1:m))
%pause
end
case {3} % multigrid
	%fun = @(h) trapezeudal_residual(@(h) swe_zero_inertia_dzdt(h,z,C,L,n), dt, h, h0);
	%fun = @(h) trapezeudal_residual(@(h) swe_zero_inertia_dzdt(h,z,C,L,n), dt, h, h0);
	returnmat = false;
	z = reshape(z,n);
	h0 = reshape(h0,n);
	h      = h0;
	hz     = h;
	hz(:,:,2) = z;
	method = 'mg';
	resfun = @(h) euler_implicit_residual(@(h) swe_zero_inertia_dzdt(h,z,C,L,n,returnmat), dt, h, h0);
	jfun   = @(L,n,nvar,hz) euler_implicit_jacobian(@(h) swe_zero_inertia_jacobian(h,hz(:,:,2),C,L,n,returnmat), dt, hz(:,:,1), returnmat);
	n 
	[h,rmse,res,iter,out] = gauss_newton(resfun,h,[],maxiter,method,jfun,L,n,nvar,z);
	out.mg.check();
case {4}
	h0 = flat(h0);
	h  = h0;
	z  = flat(z);
	method = 'mg-ks';
	% mg as preconditioner
	returnmat = true;
	resfun         = @(h) euler_implicit_residual(@(h) swe_zero_inertia_dzdt(h,z,C,L,n,returnmat), dt, h, h0);
	jfun           = @(h) euler_implicit_jacobian(@(h) swe_zero_inertia_jacobian(h,z,C,L,n,returnmat), dt, h, returnmat);
	jfun2          = @(L,n,nvar,hz) euler_implicit_jacobian(@(h) swe_zero_inertia_jacobian(h,hz(:,:,2),C,L,n,false), dt, hz(:,:,1), false);
	%[h,rmse,res,iter] = gauss_newton(resfun,h,[],maxiter,method,jfun);
	[h,rmse,res,iter] = gauss_newton(resfun,h,[],maxiter,method,jfun,L,n,nvar,reshape(z,n),jfun2);
end
t(jdx) = toc;
if (0)

tic
h = h0;
rmse = [];
m = 2e1;
for idx=1:100
	[res] = swe_zero_inertia_resfun(h,h0,z,C,dt,L,n);
	h = h-0.01*res;
	rmse(idx) = rms(res);
	[idx,rmse(idx)]
	h = max(h,0);
end % for idx
toc

end % if 0

h0 = reshape(h0,n);
h = reshape(h,n);
if (~isempty(z))
z = reshape(z,n);
end

figure(jdx);
clf
subplot(2,5,1);
imagesc(h0)
axis equal
axis tight
colorbar

subplot(2,5,2);
imagesc(h)
axis equal
axis tight
colorbar

subplot(2,5,3)
%imagesc(z)
dx = L./n;
imagesc(z)
%cdiff(z)/dx(1))
colorbar
subplot(2,5,4)
%imagesc(z)
dx = L./n;
imagesc(cdiff(z)/dx(1))
colorbar

subplot(2,5,5)
%imagesc(z)
dx = L./n;
imagesc(cdiff(cdiff(z))/dx(1)^2)
colorbar



subplot(2,4,5)
plot(x,h0(round(end/2),:))
hold on
plot(x,h(round(end/2),:))
if (0)
plot(x,max(h0,[],1))
hold on
plot(x,max(h,[],1))
plot(x,max(h,[],1)*max(h0(:))/max(h(:))) %*dx(2))
vline(L(1)/2)
yyaxis right
if (~isempty(z))
plot(x,mean(z,1))
end
end

subplot(2,4,6)
plot(x,h0(:,round(end/2),:))
hold on
plot(x,h(:,round(end/2),:))
if (0)
plot(x,max(h0,[],2)) %*dx(2))
hold on
plot(x,max(h,[],2)) %*dx(2))
plot(x,max(h,[],2)*max(h0(:))/max(h(:))) %*dx(2))
vline(L(1)/2)
yyaxis right
if (~isempty(z))
plot(x,mean(z,2))
end
end
h0_ = imrotate(h0,45,'crop','bilinear');
h_  = imrotate(h,45,'crop','bilinear');


subplot(2,4,7)
plot(x,max(h0_,[],1))
hold on
plot(x,max(h_,[],1))
plot(x,max(h_,[],1)*max(h0_(:))/max(h_(:))) %*dx(2))
%vline(L(1)/2)
%yyaxis right
%if (~isempty(z))
%plot(x,mean(z,1))
%end

subplot(2,4,8)
plot(x,max(h0_,[],2)) %*dx(2))
hold on
plot(x,max(h_,[],2)) %*dx(2))
plot(x,max(h_,[],2)*max(h0_(:))/max(h_(:))) %*dx(2))
%vline(L(1)/2)
%yyaxis right
%if (~isempty(z))
%plot(x,mean(z,2))
%end
%subplot(2,3,6)
%semilogy(rmse)
drawnow

xc(jdx,:) = [(sum(h)*x)/sum(h,'all'), (sum(h')*x)/sum(h,'all')]

end % 
end % for mode_C

%function x = precondition(x)
%end

