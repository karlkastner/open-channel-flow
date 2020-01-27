% Thu 30 Aug 18:12:46 CEST 2018
classdef Lateral_Diversion_Wide_Channel < Potential_Flow_Analytic
	properties
		% scaled velocity of diverted flow : v0 = alpha*uin
		alpha

		% scaled depth : beta ~ fs h/ws
		beta

		% secondary flow strength factor
		fs = 2/Constant.Karman^2;

		% scale of x-coordinate
		ws

		% scale of flow velocity
		uin
	end % properties
	
	methods
		function obj = Lateral_Diversion_Wide_Channel(varargin)
			obj = obj@Potential_Flow_Analytic();
			for idx=1:2:length(varargin)
				obj.(varargin{idx}) = varargin{idx+1};
			end
		end

		function Qs = Qs(obj)
			Qs = obj.uin*obj.alpha*obj.h*obj.ws;
		end

		function h = h(obj)
			 h = obj.beta*obj.ws/obj.fs;
		end

		function scales(uin,h,Qs,ws)
			
			
		end

	end % methods
end % Lateral_Diversion_Finite_Width

