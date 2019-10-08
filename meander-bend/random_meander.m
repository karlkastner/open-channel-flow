% Fri Apr 29 18:20:08 CEST 2016
% Karl Kastner, Berlin
%
%% generate a pseudo random meander
%
% note: Fr and Fi are uncorrelated
function [Fr, Fi] = random_meander(mu,n,m)
	% the first and last 10*mu bins of the fft/ifft are corrupted, so padd
	o = 10*mu;
%	o = 0;
	n_ = n+2*o;
	
	% mean frequency shifts accordingly
	mu_ = mu*n_/n;
	
	% TODO, padd with zeros until next power of 2
	
	%[a b] =logn_mode2param(20,5);
	% r=lognrnd(a,b,n,1);
	Rr = raylrnd(mu_,n_,m);
	Ri = raylrnd(mu_,n_,m);
	Hr = hist(Rr,0:n_-1);
	Hi = hist(Ri,0:n_-1);
	Hr = Hr.*sign(randn(n_,m));
	Hi = Hi.*sign(randn(n_,m));
	Hr = [Hr;  flipud(Hr)];
	Hi = [Hi; flipud(Hi)]; % do not negate second img part
	H  = Hr + 1i*Hi;
	F  = fft(H); %f = sqrt(n)*ifft(h);
	% normalise
	F  = 4/n*F;
	% extract first parts / remove redundant part
	Fr = real(F(1:n_,:));
	Fi = imag(F(1:n_,:));
	% undo padding
	Fr = Fr(o+1:o+n,:);
	Fi = Fi(o+1:o+n,:);
end

