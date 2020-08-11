% Mon  7 Jan 19:00:30 CET 2019
function [ps1,r1,ps2,r2] = sediment_division_wang(Q1,Q2,w1,w2,k)
	% S2/S1 = f = (Q1/Q2)^k(w1/w2)^(1-k)
	% (S0-S1) = S1 f
	% S0 = S1*(1+f)
	% S1 = S0/(1+f)
%	k   = 5;
	% c2/c1 = f = (Q1/Q2)^(k+1)(w1/w2)^(1-k)
	f   = 1./bsxfun(@times,(Q1./Q2).^k,(w1./w2).^(1-k));
	% ps2 = S2/S0
	ps1 = 1./(1+f);
	ps2 = 1./(1+1./f);
	pq1 = Q1./(Q1+Q2);
	pq2 = Q2./(Q1+Q2);
	% c1  = r1 = ps1/pq1;
	r1  = ps1./pq1;
	r2  = ps2./pq2;
end 

