% Thu 15 Mar 10:02:36 CET 2018
%% solve for the potential
function [c] = solve_potential(obj)
	
	if (isempty(obj.A) || isempty(obj.b))
		obj.assemble_potential_matrix();
	end

	c        = condest(obj.A);
	fprintf('Condition number is %g\n',c);

	% solve for potential
	switch (obj.solver)
	case {'gmres'}	
		% minres cannot be used without making the dirichlet condition
		% symmetric by partiall gaussian elimination
		%phi = gmres(AA,bb,sum(n),[],prod(n)); %,[],prod(n)-1,[],[],phi);
		obj.phi_ = gmres(obj.A,obj.b);
		% TODO is a SPD to use even cg inestead of minres?
		phi = gmres(A,rhs);
	otherwise
		% TODO at least reorder equations
		obj.phi_ = obj.A\obj.b;
	end

	% reset
	obj.u_ = [];
	obj.v_ = [];
	obj.R_ = [];
	obj.ubed_ = [];
	obj.vbed_ = [];
end

