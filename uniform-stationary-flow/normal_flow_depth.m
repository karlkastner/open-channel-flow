% 2016-05-17 13:37:38.441944219 +0200
% 2017-09-05 00:41:54.131459023 +0800
% Karl Kastner, Berlin
%
%% normal flow depth for uniform stationary flow
%% function H = normal_flow_depth(Q,W,C,S)
function H = normal_flow_depth(Q,W,cf,S,type)
	if (isnumeric(type))
		ismanning = type;
		if (ismanning)
			type = 'manning';
		else
			type = 'chezy';
		end
	end
	if (nargin()<5 || isempty(type))
		type = 'chezy';
	end
	switch (lower(type))
	case {'drag','cd'}
		Cz = drag2chezy(cf);
		H = (Q.^2./(Cz.^2.*W.^2.*abs(S))).^(1/3);
	case {'chezy','cz'}
		Cz = cf;
		H = (Q.^2./(Cz.^2.*W.^2.*abs(S))).^(1/3);
	case {'manning','n'}
		n = cf;
		H = ( (Q.*n)./(sqrt(S).*W) ).^(3/5);
	otherwise
		error('here');
	end
end

