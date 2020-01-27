% Sat  1 Sep 11:13:38 CEST 2018
% Karl Kästner, Berlin

% t  : time
% xy : time
% pf : potential flow object
% f  : potential flow functions
function [t,xy,pf,f,h] = streamline(obj,shape,alpha,beta)
	if (nargin() < 4)
		key    = [shape,' ',num2str(alpha)];
		ufield = 'u';
		vfield = 'v';
		beta = [];
	else
		% when beta is given, determine streamline along bed
		key    = [shape,' ',num2str([alpha,beta])];
		ufield = 'ubed';
		vfield = 'vbed';
	end
	if (~obj.map.isKey(key))
		obj.change = true;
		%pf  = Potential_Flow_Analytic();
		pf  = Lateral_Diversion_Wide_Channel();
		f   = pf.lateral_outflow(shape,alpha,beta);
		pf.alpha = alpha;
		pf.beta  = beta;

		e = obj.erel;
		x0  = pf.fun.y0([],[],alpha,[]);
		y0  = e*alpha;
	
		% integration limit
		e = obj.erel;
		l = alpha*tan(pi/2 - pi*e);
		u0 = 1;
		T = 2*l/u0;
	
		% TODO choose left limit so that int(Q,0,Q/u0)- Qs < 1e-3 Qs
		xlim_ = [-l,1]; %ws];
		ylim_ = [-inf,inf]; %-sqrt(eps)*ws];
		opt = odeset(...
		        'event', @(t,xy) boxevent(t,xy,xlim_,ylim_));

		[t,xy]    = pf.streamline([0,-T],[x0,y0],ufield,vfield,opt);
		%[t,xy]    = pf.streamline([0,-T],[(1-e)*x0,y0],ufield,vfield,opt);
		
		% relative error
		% TODO this is not a good error estimate
		erel = abs((xy(end,2)-xy(end-1,2))/xy(end,2));

		if (0)
			[t_,xy_]    = pf.streamline([0,-T],[(1+e)*x0,y0],ufield,vfield,opt);
			%erel(2) = (xy_(end,2)-xy(end,2))/(xy_(end,2)+xy(end,2))
			erel(2) = (abs(xy_(end,2))-abs(xy(end,2)))/(abs(xy_(end,2))+abs(xy(end,2)));
		end			
	
		sol.pf = pf;
		sol.f  = f;
		sol.t  =  t;
		sol.xy = xy;
		sol.erel  = erel;
		obj.map(key) = sol;
		%sol.xy_ = xy_;
		%sol.t_ = t_;
	else
		% load precomputed streamline
		sol = obj.map(key);
		t   = sol.t;
		f   = sol.f;
		xy  = sol.xy; % ws
		%xy  = ws*sol.xy;
		pf   = sol.pf;
		erel = sol.erel;
	end
end % function map

