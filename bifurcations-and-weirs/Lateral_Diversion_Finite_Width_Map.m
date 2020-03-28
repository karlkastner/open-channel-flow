% Mon 27 Jan 10:36:44 +08 2020

% TODO, interpolation
% TODO, order
classdef Lateral_Diversion_Finite_Width_Map < handle
	properties
		matfile
%		keymap
		matfilename
	end % properties
	methods
		function obj = Lateral_Diversion_Finite_Width_Map(matfilename)
			if (nargin()>0)
				obj.matfilename = matfilename;
			end
		end

		function key = key(obj,alpha,beta,gamma)
			key1 = sprintf('%0.3e%0.3e%0.3e',[cvec(alpha),cvec(beta),cvec(gamma)]');
			key = reshape(key1,27,[])';
		end % key

		function index = index(obj,key)
			index = zeros(size(key,1),1);
			map = obj.matfile.map;
			for idx=1:size(key,1)
				if (map.isKey(key(idx,:)))
					% load index of precomputed streamline
					index(idx) = map(key(idx,:));
				end
			end
		end % index

		function dsbi = streamline(obj,a,b,g)
			siz = size(a);
			a = flat(a);
			b = flat(b);
			g = flat(g);

			key_a   = obj.key(a,b,g);
			index   = obj.index(key_a);
			flag    = (index > 0);
			% fetch precomputed streamlines
			fdx     = find(flag);
			for idx=1:length(fdx)
				dsbi(fdx(idx)) = obj.matfile.dsb(1,index(fdx(idx)));
			end
			% compute non-existing streamlines, vectorized
			% TODO process blockwise
			if (obj.matfile.Properties.Writable)
			fdx = find(~flag);
			if (length(fdx)>0)
				l       = Lateral_Diversion_Finite_Width('n',obj.matfile.kmax);

				l.alpha = a(fdx);
				l.beta  = b(fdx);	
				l.gamma = g(fdx);

				x0 = l.stagnation_point();
				y0 = obj.matfile.y0*l.alpha;

				tic();
				% TODO, no magic numbers
				tol = 1e-3;
				opt      = odeset('RelTol',tol);
				[t,xy] = l.streamline(obj.matfile.T,[x0,y0],'bed',opt);
				toc();
				x = xy(:,1:end/2);
				y = xy(:,end/2+1:end);
				tic()
				map = obj.matfile.map;
				last = obj.matfile.last;
				for idx=1:length(fdx)
					dsbi(fdx(idx)).alpha = l.alpha(idx);
					dsbi(fdx(idx)).beta  = l.beta(idx);
					dsbi(fdx(idx)).gamma = l.gamma(idx);
					dsbi(fdx(idx)).t     = t;
					dsbi(fdx(idx)).x     = x(:,idx);
					dsbi(fdx(idx)).y     = y(:,idx);
					% these are redundant
					%dsbi(fdx(idx)).x0    = x0(idx);
					%dsbi(fdx(idx)).y0    = y0(idx);
					%dsbi(fdx(idx)).T     = T;
					% store
					last                    = last+1;
					map(key_a(fdx(idx),:))  = last;
					obj.matfile.dsb(1,last) = dsbi(fdx(idx));
				end % for idx
				obj.matfile.map  = map;
				obj.matfile.last = last;
				toc()
			end % if
			end

			% reshape
			dsbi = reshape(dsbi,siz);
		end % streamline

		% TODO, this does not remove existing dsb
		function obj = create(obj,kmax,T,y0)
			obj.matfile = matfile(obj.matfilename,'Writable',true);
			
			obj.matfile.map  = containers.Map();
			obj.matfile.last = 0;
			obj.matfile.kmax = kmax;
			obj.matfile.T    = T;
			obj.matfile.y0   = y0;
		end % create

		function [flag, obj]  = load(obj,wflag)
			if (nargin()<2)
				wflag = true;
			end
			if (   ~isempty(obj.matfilename) ...
			     && exist(obj.matfilename,'file'))
				obj.matfile = matfile(obj.matfilename,'Writable',wflag);
				flag = true;
			else
				flag = false;
			end
		end % load()
	end % methods
end % Lateral_Diversion_Finite_Width_Map

