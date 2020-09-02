Qlim=[1e3;1e4]; q = []; for idx=1:6; [w,b]=int_1d_gauss(idx); Q=inv_hydrograph(b*[0;1],Qlim(1),Qlim(2)); q(idx)=w'*Q.^2, end, qtrue=quad(@(t) (1e3+(1e4-1e3)*1/2*(1+cos(2*pi*t))).^2,0,1 )

