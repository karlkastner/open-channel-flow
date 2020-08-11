% Wed  5 Sep 08:26:06 CEST 2018
%
% for alpha > 0.25
% r : 0.49/alpha + 1.4
% function swdr = sediment_division_raudkivi(alpha)
function swdr = sediment_division_raudkivi(alpha)
	if (~issym(alpha))
	ab = 0.4
	bb = 1.15
	as = 0.81
	else
		syms ab bb as
	end
	
	% dummy width
	ws = 1;

	% width at surface (eq 3.9)
	Bs = ws*(as*alpha);

	% width at bottom (3.10)
	Bb = ws*(ab + bb*alpha);

	% width of diverted flow in the approaching channel, by definition
	Bq = ws*alpha;

	% ratio of near bed to depth averaged streamline
	% i.e. sediment to water division ratio
	swdr = Bb./Bq;
end

