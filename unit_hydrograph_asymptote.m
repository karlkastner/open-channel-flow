% 2024-06-18 15:06:21.871118350 +0200

	% approximation for large n
	% gammaln ~ n*ln(n) - n - 1/2*ln(n) + 1/*log(2*pi)
	%Q_ = 1./k.*exp(-n.*log(n) + n + 1/2*log(n) - 1/2*log(2*pi) + log(t./k).*(n-1) - t./k);
	%k = K/n;
	K = k.*n;
	Q_ = n./K.*exp(-n.*log(n) + n + 1/2*log(n) - 1/2*log(2*pi) + log(t.*n./K).*(n-1) - t.*n./K);
	Q_ = n./K.*exp(-n.*log(n) + n + 1/2*log(n) - 1/2*log(2*pi) + (log(t./K) + log(n)).*(n-1) - t.*n./K);
	Q_ = n./K.*exp(            +n + 1/2*log(n) - 1/2*log(2*pi) +  log(t./K).*(n-1) - log(n)  - t.*n./K);
	Q_ = n./K.*exp(            +n + 1/2*log(n) - 1/2*log(2*pi) +  log(t./K).*(n-1) - log(n)  - t.*n./K);
	Q_ = n./K.*exp(            +n - 1/2*log(n) - 1/2*log(2*pi) +  log(t./K).*(n-1)           - t.*n./K);
	Q_ = 1./K.*exp( log(n)     +n - 1/2*log(n) - 1/2*log(2*pi) +  log(t./K).*(n-1)           - t.*n./K);
	Q_ = 1./K.*exp( 1/2*log(n) +n - 1/2*log(2*pi) + (log(t./K)).*(n-1) - t.*n./K);
	Q_ = 1./K.*exp( 1/2*log(n) +n - 1/2*log(2*pi) + (log(t./K)-t./K).*(n-1) - t./K);
	% series expansion of log(t/K) near t=K
	Q_ = 1./K.*exp( 1/2*log(n) +n - 1/2*log(2*pi) + (log(t./K)-t./K).*(n-1) - t./K);
	Q_ = 1./K.*exp( 1/2*log(n) +n - 1/2*log(2*pi) + (log(1+(t-K)./K)-t./K).*(n-1) - t./K);
	Q_ = 1./K.*exp( 1/2*log(n) +n - 1/2*log(2*pi) + ((t-K)./K - 1/2*((t-K)./K).^2 -t./K).*(n-1) - t./K);
	Q_ = 1./K.*exp( 1/2*log(n) +n - 1/2*log(2*pi) + (-1 - 1/2*((t-K)./K).^2).*(n-1) - t./K);
	Q_ = 1./K.*exp( 1/2*log(n) +n - 1/2*log(2*pi) + -(n-1) - (1/2*((t-K)./K).^2).*(n-1) - t./K);
	Q_ = 1./K.*exp( 1/2*log(n) +1 - 1/2*log(2*pi) + - (1/2*((t-K)./K).^2).*(n-1) - t./K);
	Q_ = 1./(sqrt(2*pi).*K).*exp( 1/2*log(n) + 1 - (1/2*((t-K)./K).^2).*(n-1) - t./K);
	% near t=K
	Q_ = 1./(sqrt(2*pi).*K).*exp( 1/2*log(n) - (1/2*((t-K)./K).^2).*(n-1));
	Q_ = sqrt(n)./(sqrt(2*pi).*K).*exp(-1/2*((t-K)./K).^2.*(n-1));
	Q_ = sqrt(n)./(sqrt(2*pi).*k.*n).*exp(-1/2*((t-k.*n)./(k.*n)).^2.*(n-1));
	% as n-1 ~ n
	% so, the unique hydrograph takes the shape of a gaussian,
	% with zero at tc = k*n, and standard deviation k*sqrt(n)
	% so if reservoirs are split into smaller reservoirs with k_n = k1/n
	% we have tc = k_n*n = k_1*n/n = k_1, and sd = k_n*sqrt(n) = k_1/sqrt(n)
	Q_ = 1./(sqrt(2*pi*n).*k).*exp(-1/2*((t-k.*n)./(k.*sqrt(n))).^2);
	% so k_1 is a proxy for celerity
	% and k_1/sqrt(n) is a proxy for diffusion

	%Q = n./K.*exp(log(t/K).*n + n + 1/2*log(n) - t.*n./K);
	%Q_ = n./K.*exp(log(t./K).*n + n + log(sqrt(n)) + 1/2*log(2*pi) - t.*n./K);

