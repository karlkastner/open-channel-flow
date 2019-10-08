% Sat  1 Sep 11:13:38 CEST 2018
% Karl Kästner, Berlin
%%
%% wrapper to store precomputed streamlines of potential flows
%%
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
classdef Potential_Flow_Map < handle
	properties
		map
		filename = '';
		change = false;
		%erel = 2.5e-4;
		erel = 1e-3;
		sfilename = 'streamlines.mat';
	end
	methods
		function obj=Potential_Flow_Map()
			obj.map = containers.Map();
		end

		function streamline_dimensional()
			% TODO when ws is scaled / set, xy0 must be scaled as well
			ws = 1;
			%h  = 1;
			u0 = 1;
			if (nargin() < 4)
				beta = 1;
				key    = [shape,' ',num2str(alpha)];
				ufield = 'u';
				vfield = 'v';
			else
				% when beta is given, determine streamline along bed
				key    = [shape,' ',num2str([alpha,beta])];
				ufield = 'ubed';
				vfield = 'vbed';
			end
			fs = 2/0.41^2;
			h  = beta*ws/fs
			Qs = alpha*u0*ws*h;
				%ws = h*fs/beta;
				%h = (0.41^2/2)*Qs/u0*beta;
				%h_div_ws = (0.41^2/2)*alpha*beta;

		end

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
				pf  = Potential_Flow_Analytic();
				f   = pf.lateral_outflow(shape,alpha,beta); %Qs,ws,shape);
	
				switch (1)
				case {1}
					e = obj.erel;
					%x0 = 0.999*f.(shape).fun.y0([],[],u0,Qs,ws);
					%x0  = (1-e)*f.(shape).fun.y0([],[],alpha,[]); %u0,Qs,h,ws);
					%x0_ = (1+e)*f.(shape).fun.y0([],[],alpha,[]); %u0,Qs,h,ws);
					x0  = (1-e)*pf.fun.y0([],[],alpha,[]); %u0,Qs,h,ws);
					x0_ = (1+e)*pf.fun.y0([],[],alpha,[]); %u0,Qs,h,ws);
					y0  = e*alpha; % war -e
				case {2}
					y0=-sqrt(eps)*alpha;
					[x0,fval] = fzero(@(x) pf.ubed(x,y0),1) %1/2)
					x0 = (1+sqrt(eps))*x0 + sqrt(eps);
					pf.vbed(x0,y0)
				case {3}
					x0  = f.(shape).fun.y0([],[],u0,Qs,h,ws)
					% 1e-4 is too low
					%[y0,fval] = fzero(@(y) pf.vbed(x0,y),-10*alpha);
					if (~isempty(h))
					x0   = x0+sqrt(eps);
					y0   = alpha;
					vbed = pf.vbed(x0,y0)
					while (pf.vbed(x0,y0/2)>0)
						y0 = y0/2;
					end
					if (0)
					y0/alpha
					vbed = pf.vbed(x0,y0)
					vbed = pf.vbed(x0,y0/2)
	
					figure(1)
					clf
					y = linspace(-1,0,1e3);
					%plot(y,pf.vbed(x0,y)); hold on
					plot(y,pf.vbed(0,y));
					%plot(y,pf.vbed(1/2*0.99,y));
					%vline(y0)
					pause
					end
					% = 2*y0;
					else
						y0 = sqrt(eps)*alpha;
						y0 = 1e-3*alpha;
					end
					%$alpha/y0
					%pf.vbed(x0,y0)
					%pf.vbed(xy0(1),xy0(2))
					%pf.vbed(x0,1.1*y0)
					%pf.vbed(x0,-1e-3*alpha)
					%pf.vbed(x0,-1e-2*alpha)
					%pf.vbed(x0,-1e-1*alpha)
					%xy0 = [x0,-1e-3*alpha];
					%xy0 = [x0,1.001*y0];
				end % switch
				xy0 = [x0,y0];
				%-1e-4*ws];

				% integration limit
				e = obj.erel;
				l = alpha*tan(pi/2 - pi*e);
				%l = (alpha*ws)/(sqrt(e*pi))
				T = 2*l/1;% u0;
				%[T,l]
				%l  = 1e3*ws;
				%T    	 = 20000;

				% TODO choose left limit so that int(Q,0,Q/u0)- Qs < 1e-3 Qs
				xlim_ = [-l,1]; %ws];
				ylim_ = [-inf,inf]; %-sqrt(eps)*ws];
				opt = odeset(...
					... % 'RelTol',sqrt(eps), ... fails sometimes
				        'event', @(t,xy) boxevent(t,xy,xlim_,ylim_));
				%ufield = 'u'; vfield = 'v';
				[t,xy]    = pf.streamline([0,-T],xy0,ufield,vfield,opt);
				if (0)
					figure(1);
					clf
					[t_,xy_]    = pf.streamline([0,-T],xy0,'u','v',opt);
					plot(xy(:,1),xy(:,2));
					hold on
					plot(xy_(:,1),xy_(:,2));
					plot(xy0(1),xy0(2),'o');
					hline(alpha)
					pause(1)
				end
				
				% relative error
				erel = abs((xy(end,2)-xy(end-1,2))/xy(end,2));
	
				if (1)
					[t_,xy_]    = pf.streamline([0,-T],[x0_,y0],ufield,vfield,opt);
					%erel(2) = (xy_(end,2)-xy(end,2))/(xy_(end,2)+xy(end,2))
					erel(2) = (abs(xy_(end,2))-abs(xy(end,2)))/(abs(xy_(end,2))+abs(xy(end,2)));
				end			

				if (0)
					[t_,xy_] = f.streamline([0,+T],xy0,'ubed','vbed',opt);
					t        = [flipud(t_);t];
					xy       = [flipud(xy_bed); xy_bed_];
				end
				sol.pf = pf;
				sol.f  = f;
				sol.t  =  t;
				sol.xy = xy/1; %ws;
				sol.erel  = erel;
				obj.map(key) = sol;
				sol.xy_ = xy_;
				sol.t_ = t_;
			else
				sol = obj.map(key);
				t   = sol.t;
				f   = sol.f;
				xy  = sol.xy; % ws
				%xy  = ws*sol.xy;
				pf   = sol.pf;
				erel = sol.erel;
			end
		end % function map

		function obj = load(obj)
			what_ = what(class(obj));
			path_ = [what_.path,filesep(),obj.sfilename];

			%if (nargin()>0)
			%	obj.filename = filename;
			%end
			%obj.slb = containers.Map();
			% obj.map  = containers.Map();

			if (exist(path_,'file'))
				load_ = load(path_,'map');
				obj.map = load_.map;
			else
				obj.map = containers.Map();
			end
			obj.change = false;
		end % load()
		function obj = save(obj)
			if (obj.change)
				what_ = what(class(obj))
				path_ = [what_.path,filesep(),obj.sfilename]
				map   = obj.map;
				save(path_,'map');
				obj.change = false;
			end
		end % save()
	end % methods()
end % Potential_Flow_Map

