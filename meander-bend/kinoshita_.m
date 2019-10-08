% 2015-09-29 14:36:03.575085215 +0200
% Bart Vermeulen
%
% angle in planform teta(s) = teta0 cos phi − teta0^3 ( c_f cos(3 phi) + c_s sin(3 phi))
%
% curvature : C(s) = dteta(s)/ds
%
function [x,y,theta,c]=kinoshita(s,L,cs,cf,theta0,varargin)


validateattributes(s,{'numeric'},{'nonempty','finite','vector'},'kinoshita','s',1);
validateattributes(L,{'numeric'},{'nonempty','finite','vector','positive'},'kinoshita','L',2);
validateattributes(cs,{'numeric'},{'nonempty','finite','vector'},'kinoshita','cs',3);
validateattributes(cf,{'numeric'},{'nonempty','finite','vector'},'kinoshita','cf',4);
validateattributes(theta0,{'numeric'},{'nonempty','finite','vector','positive'},'kinoshita','theta0',6);

P=inputParser;
P.parse(varargin{:});
delete(P)

L=vectinput(L,s);
cs=vectinput(cs,s);
cf=vectinput(cf,s);
theta0=vectinput(theta0,s);

if isrow(s), vecdim=2; else vecdim=1; end

sth=(s(1:(end-1))+s(2:end))/2;
L=(L(1:(end-1))+L(2:end))/2;
cs=(cs(1:(end-1))+cs(2:end))/2;
cf=(cf(1:(end-1))+cf(2:end))/2;
theta0=atan2((sin(theta0(1:(end-1)))+sin(theta0(2:end)))/2,(cos(theta0(1:(end-1)))+cos(theta0(2:end)))/2);

phi=2*pi*sth./L;


theta=theta0.*cos(phi)-theta0.^3.*(cf.*cos(3*phi) + cs.*sin(3*phi));

c=diff(theta)./diff(sth);

dx=cos(theta).*diff(s);
dy=sin(theta).*diff(s);
 
x=cumsum(cat(vecdim,0,dx));
y=cumsum(cat(vecdim,0,dy));


function y=vectinput(x,s)

y=ones(size(s)).*x;

