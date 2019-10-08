% 2015-02-23 11:42:21.016246040 +0100 test_PowerRatingCurve.m


n = 100;
h = 10*(1:n)'/(n+1);
if (0)
 c = [1 1 1];
 q = c(1) + c(2)*h + c(3).^h.^2;
 rc  = PolyRatingCurve(2);
else
 c = [1 0 1.5];
 q = c(1)*(h - c(2)).^c(3);
 rc  = PowerRatingCurve();
end
q0 = q + 0.01*mean(q)*randn(size(q));
rc.fit(h,q0);
rc.param.val
rc.param.val0

rc.lin.C
g_num = grad(@(c) rc.rcfunc(c,[],h)-q0 ,rc.param.val0)
[f g]= rc.rcfunc_(rc.param.val0, [], h) 
%g = grad(@(c) norm(rc.rcfunc(c,[],h)-q0) ,rc.param.val0)
% note f = g*c = A*c is only true for linear functions !!!

powerrc=rc;

clf
N = nchoosek(1:3,2);
C0 = powerrc.lin.C;
C1 = powerrc.param.C;
C1 = squeeze(C1);
%C1 = C1*C1';
p0 = powerrc.param.val0; 
p1 = squeeze(powerrc.param.val1);
cc = colormap('lines');
for idx=1:3
	subplot(1,3,idx);
	id = N(idx,1);
	jd = N(idx,2);
	plot(p0(id), p0(jd),'.b');
	% plot estimated parameter by linearisation
	c = [C0(id,id), C0(id,jd), C0(jd,jd), p0(id), p0(jd)];
	hold on
	% plot parameter error ellipse for linearisation
	plot_ellipse(c);
	c = [C1(id,id), C1(id,jd), C1(jd,jd), p0(id), p0(jd)]
	plot_ellipse(c,[],'color',cc(2,:));
	% plot jackknife parameter
	plot(p1(id,:),p1(jd,:),'.','color',cc(2,:));
	xlabel(id)
	ylabel(jd)
end

