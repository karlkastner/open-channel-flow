% Mon 20 Jan 20:08:31 +08 2020
function val = evalk(obj,fname,x,y)
	val = 0;
	for k=-obj.n:obj.n
		val = val + obj.fun.k.(fname)(x,y,obj.alpha,obj.beta,obj.gamma,k);
	end
%		function u = u(obj,name,x,y)
%			u = obj.evalk(obj.fun.(name)
%		end
end


