% 2024-06-18 15:06:22.851137405 +0200
function [Q,Q_] = unit_hydrograph(t,k,n)
	%Q = 1./(k.*gamma(n)).*(t./k).^(n-1).*exp(-t./k);
	%Q = 1./(k.*gamma(n)).*exp(log(t./k).*(n-1)-t./k);
	Q = 1./k.*exp(-gammaln(n)+log(t./k).*(n-1)-t./k);
end

