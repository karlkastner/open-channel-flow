% Thu 15 Mar 10:02:36 CET 2018
%% solve for the flow potential
% TODO rename determine
function [c, obj] = solve_potential(obj)
	% matrix is not symmetric, minres cannot be used
	%phi = gmres(AA,bb,sum(n),[],prod(n)); %,[],prod(n)-1,[],[],phi);

	c = condest(obj.A);
	fprintf('Condition number is %g\n',c);
	obj.phi_ = obj.A\obj.b;

	% reset
	obj.ubed_ = [];
	obj.vbed_ = [];
end

