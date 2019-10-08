% Thu 22 Mar 11:48:57 CET 2018
% Karl Kastner, Berlin
%
%% calibrate bend geometry to given profile
function obj = calibrate(obj,Q,Qs)
	% match given river discharge and sediment transport by
	% adjusting centreline depth and mean velocity	

	val = [obj.Hc,obj.u_bar];
	opt.RelTol = 1e-2;
	val = lsqnonlin(@fobj,val,[],[],opt);

	obj.Hc    = val(1);
	obj.u_bar = val(2);
	%[r,H,D]   = obj.bed_profile();

	function res = fobj(val)
		obj.Hc    = val(1);
		obj.u_bar = val(2);

		[r,h,d] = obj.bed_profile();
	
		% integrate discharge
		u   = obj.u(r);
		Q_ = sum(mid(u.*h).*diff(r));
		%Q_ = centre(u.*h)'*diff(r);
	
		% integrate sediment transport
		qs  = total_transport_engelund_hansen(obj.C,1e3*d,u,h,1);
		Qs_ = sum(mid(qs).*diff(r));
		res = [Q_  - Q;
		       Qs_ - Qs];
	end	
end % calibrate

