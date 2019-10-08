% 2015-06-08 20:12:01.962718678 +0200
% Karl Kastner, Berlin

%% transverse profile of the streamwise velocity in a meander bend

% this is equation 2.69 in olesen (1987),
% who adapted it from Kalkwijk
% - his choice of q is to small for the wide bend in sanggau
% - fails to predict change of gradient with the Dean number (as computed by de Vriend)
%
% this does not change with H/R, does the author implicitely assume uniform flow? (infinitely long bend?)

% radius of curvature (from Sanggau bend)
Rc = 1500;	
% width (from Sanggau bend)
w = 600;
% depth (from Sanggau bend)
%h = 7;
H = linspace(5,12,3);

% acelleration by gravity
g = 9.81;
% Chezy (good guess)
C  = 50;
% spanwise coordinate
n  = linspace(-w/2,w/2,100);

% side wall adaptation factor, olesen chooses 2, but this is clearly too small,
% as he writes himself, q should depend on depth and width
%q = 2;
q = 6;
% some factor, explained more detailed in Kalkwijk
ksn = 0.4;
u = [];
for idx=1:length(H)
	h=H(idx);
	u(:,idx) = (1-n/Rc)./sqrt(1 - 2*q*ksn*C^2/g*(1/Rc)*(1/w)*sign(n).*h.^2.*(2*n/w).^(q-1));
end
plot(n,u);

