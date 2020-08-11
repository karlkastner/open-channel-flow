% Wed  8 Aug 20:52:28 CEST 2018
% Wed  5 Sep 08:54:14 CEST 2018
% function r = sediment_division_meijer_ksiazek(Qin,Qs,w0,ws)
function swdr = sediment_division_meijer_ksiazek(Qin,Qs,w0,ws)
	if (issym(Q0))
		syms m n
	else
		m = 1.13;
		n = 0.39;
	end

	% constant main channel width
	wd = w0;
	% residual discharge in downstream main branch
	Qd         = Qin - Qs;
	% ratio of sediment in downstream and side branch
	Ss_div_Sd  = (Qs./Qd).^m.*(wd./ws).^n;
	% ratio of sediment in approaching and side branch
	Ss_div_Sin = Ss_div_S1./(1 + Ss_div_Sd);
	% sediment in approaching channel (dummy)
	c  = 1;
	% approaching sediment (dummy)
	Sin = c.*Qin;
	% diverted sediment
	Ss  = Ss_div_Sin.*Sin;
	% sediment to water division ratio
	swdr = (Qin./Sin).*(Ss./Qs);
end

