% Thu 23 Apr 15:54:16 +08 2020
function obj = parse_arguments(obj,varargin)
	for idx=1:2:length(varargin)
		f   = varargin{idx};
		v   = varargin{idx+1};
		obj = setfield_deep(obj,f,v);
	end
end

