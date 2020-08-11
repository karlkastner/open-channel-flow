% Wed 26 Feb 21:01:47 +08 2020
%
% predict discharge through a cross-section with polynomial shape
%
% predict bedforms and roughness from flow conditions and bed material
%
% zb : bed level at points spaced equally distributed across section between [-1/2, 1/2]*w
% 
% symmetric cross section, only half the cross section is discretized and then mutiplied with 2
function Q = compute_discharge(obj,setflag,varargin)
	if (nargin()<2)
		setflag = true;
	end

	obj = parse_arguments(obj, varargin{:});

%	zs = obj.zs;
%	nt = length(zs);
	if (obj.force_ordered ...
		&& any(diff(obj.zb.c.param)<0))
		%&& any(diff(obj.n.zb(1,:))<0))
		Q = NaN(size(obj.zs));
		return;
	end
%	zb = sort(zb);

	% reassign grain size
	obj.init();

	[Q, qn, Cn, d50n, d90n, zbn, weight, area, width_,xn] = discharge_(obj);

	if (obj.verbose)
		disp([num2str(rvec(obj.zb.c.param)),' ',num2str(mean(Q))])
		%disp([num2str(rvec(obj.n.zb(1,:))),' ',num2str(mean(Q))])
	end

	if (setflag)
		obj.discharge = Q;
		obj.area     = area;	
		obj.width_   = width_;

		obj.n.x      = xn;
		obj.n.zb     = zbn;
		obj.n.weight = weight;
		obj.n.q      = qn;
		obj.n.C      = Cn;
		obj.n.zb     = zbn;
		obj.n.d50    = d50n;
		obj.n.d90    = d90n;
	end
end % function

% Wed 26 Feb 21:01:47 +08 2020
function [Q, qn, Cn, d50n, d90n, zbn, weightn, area, width_,xn] = discharge_(obj)

	nn = obj.nn;
	zs = obj.zs;
	nt = length(zs);
	S  = obj.slope;

	% discretize cross section
	[zbn, d50n, d90n, weightn, width_, xn] = obj.wetted_cross_section();


	% depth
	h     = cvec(obj.zs)-zbn;
	h     = max(0,h);

	if (~isscalar(S))
	      S = repmat(S,1,nn);
	end % if

	if (strcmp(obj.method,'fixed-profile'))
		%zs_ = flat(repmat(zs,1,obj.nn));
		% fit roughness
		[n]            = Compound_Cross_Section.roughness(obj.discharge,zs,zbn,width_,S,true);
		%cs(idx).n.C
		Cn             = manning2chezy(n,h);
		%Cn = obj.n.C;
		[Q, qn] = Compound_Cross_Section.discharge(zs,zbn,1,Cn,S,false);
	else
		[qn,Cn] = stage_discharge(h,1,S,d50n,d90n,obj.T_C,'log',obj.method);
	end

	% reshape discharge to time x cross-section
%	Cn      = reshape(Cn,nt,nn);
%	hn      = reshape(h,nt,nn);
%	qn      = reshape(qn,nt,nn);

	% integrate discharge across section
	Q       = width_.*sum(qn.*weightn,2);
	area    = width_.*sum(h.*weightn,2);
end

