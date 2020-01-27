syms Q W S n; U = normal_flow_velocity(Q,W,n,S,1); H = normal_flow_depth(Q,W,n,S,1); isequal(Q,U*H*W)
