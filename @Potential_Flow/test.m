function test(obj)
	class(obj)
	if (isa(obj,'Potential_Flow_Analytic'))
		'yes'
	end
end
