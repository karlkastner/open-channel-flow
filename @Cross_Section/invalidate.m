% Fri 24 Apr 16:49:56 +08 2020
function obj = invalidate(obj,fdx)
	ff = {'n', 'Qs'};
	for idx=1:length(ff)
		f = fieldnames_deep(obj.(ff{idx}));
		for jdx=1:length(f)
			v = getfield_deep(obj.(ff{idx}),f{jdx});
			if (size(v,1) == length(fdx))
				v(fdx,:) = NaN;	
			end
			obj.(ff{idx}) = setfield_deep(obj.(ff{idx}),f{jdx},v);
		end
	end
	f = {'area','discharge'};
	for idx=1:length(f)
		if (length(obj.(f{idx})) == length(fdx))
			obj.(f{idx})(fdx) = NaN;
		end
	end
	
end

