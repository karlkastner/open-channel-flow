% 2015-09-29 14:36:03.575085215 +0200
% Bart Vermeulen

classdef Kinoshita < handle
    % KINOSHITA Kinoshita curves generator
    %
    %   Kinoshita Properties:
    %       s - along channel coordinate
    %       L - meander length
    %       cs - skewing parameter (typically between -0.1 and 0.1)
    %       cf - fattening paramter (typically between -0.07 and 0.1)
    %       theta0 - Planform angle (typically between pi/12 and pi/2)
    %       W - Width of the planform
    %       x - planform x coordinate (read only)
    %       y - planform y coordinate (read only)
    %       c - planform curvature (read only)
    %       theta - planform angle (read only)
    %
    %   Kinoshita Methods:
    %       plot - plot the Kinoshita curve
    %       polygon - returns a Polygon object given a width of the channel
    
    %%% Public properties
    properties(SetObservable)
        s % along channel coordinate
        L % meander length
        cs % skewing parameter (typically between -0.1 and 0.1). This parameter should either be scalar or match the size of s
        cf % fattening paramter (typically between -0.07 and 0.1). This parameter should either be scalar or match the size of s
        theta0 % Planform angle amplitude in radians (typically between pi/12 and 2/3 pi). This parameter should either be scalar or match the size of s
    end % properties
    
    %%% Public get properties
    properties(Dependent)
        x % x coordinate of kinoshita curve
        y % y coordinate of kinoshita curve
        c % curvature of kinoshita curve
        theta % angle of kinoshita curve with horizontal
    end
    
    %%% Private properties
    properties(Access=protected,Dependent)
        cs_vec % vectorized cs (same size as s)
        cf_vec % vectorized cf (same size as s)
        theta0_vec % vectorized theta0 (same size as s)
        L_vec % vectorized L (same size as s)
        s_c % centered s (average of consecutive values)
        L_c % centered L (average of consecutive values)
        cs_c % centered cs (average of consecutive values)
        cf_c % centered cf (average of consecutive values)
        theta0_c  % centered theta0 (average of consecutive values)
    end
    
    methods
        %%% Constructor
        function obj=Kinoshita(varargin)
            if nargin==5
                if all(cellfun(@iscell,varargin))
                    sizin=cellfun(@size,varargin,'UniformOutput',false);
                    assert(isequal(sizin{:}),'Kinoshita:bad_input','All input cells must have equal size')
                    sizin=sizin{1};
                    obj(numel(varargin{1}))=Kinoshita;
                    obj=reshape(obj,sizin);
                    [obj(:).s]=varargin{1}{:};
                    [obj(:).L]=varargin{2}{:};
                    [obj(:).cs]=varargin{3}{:};
                    [obj(:).cf]=varargin{4}{:};
                    [obj(:).theta0]=varargin{5}{:};
                else
                    obj.s=varargin{1};
                    obj.L=varargin{2};
                    obj.cs=varargin{3};
                    obj.cf=varargin{4};
                    obj.theta0=varargin{5};
                end % if all input is cell
            else % if nargin == 5
                [obj.s, obj.L, obj.cs, obj.cf, obj.theta0]=deal(double.empty(0,1));
            end % if nargin=5
        end
    end
    
    methods
        %%% Setters and getters
        % s
        function val = get.s(obj), val = obj.s; end % get method for s
        function set.s(obj,val) % set method for s
            validateattributes(val,{'numeric'},{'finite','vector'},'Kinoshita','s'); % check variable
            val=val(:); % vectorize
            if ~issorted(val), val=sort(val); end % sort if necessary
            obj.s=val; % assign value to property
        end % set.s
        
        % L
        function val = get.L(obj), val = obj.L; end % get method for s
        function set.L(obj,val) % set method for s
            validateattributes(val,{'numeric'},{'finite','vector','positive'},'Kinoshita','L'); % check variable
            obj.L=val(:); % assign value to property
        end % set.s
        
        % cs
        function val = get.cs(obj), val = obj.cs; end % get method for cs
        function set.cs(obj,val) % set method for cs
            validateattributes(val,{'numeric'},{'finite','vector'},'Kinoshita','cs'); % check input
            obj.cs=val(:); % assign value to property
        end % set.cs
        
        % cf
        function val = get.cf(obj), val = obj.cf; end %
        function set.cf(obj,val) % set method for cf
            validateattributes(val,{'numeric'},{'finite','vector'},'Kinoshita','cf'); % check input
            obj.cf=val(:); % assign value to property
        end % set.cs
        
        % theta0
        function val = get.theta0(obj), val = obj.theta0; end %
        function set.theta0(obj,val) % set method for theta0
            validateattributes(val,{'numeric'},{'finite','vector','positive'},'kinoshita','theta0');
            obj.theta0=val(:); % assign value to property
        end % set.theta0
        
        % s_c
        function val=get.s_c(obj), val=obj.consec_mean(obj.s); end
        
        % L_vec
        function val=get.L_vec(obj), val=obj.vectinput('L'); end
        
        % L_c
        function val=get.L_c(obj), val=obj.consec_mean(obj.L_vec); end
        
        % cf_vec
        function val=get.cf_vec(obj), val=obj.vectinput('cf'); end
        
        % cf_c
        function val=get.cf_c(obj), val=obj.consec_mean(obj.cf_vec); end
        
        % cs_vec
        function val=get.cs_vec(obj), val=obj.vectinput('cs'); end
        
        % cs_c
        function val=get.cs_c(obj), val=obj.consec_mean(obj.cs_vec); end
        
        %theta0_vec
        function val=get.theta0_vec(obj), val=obj.vectinput('theta0'); end
        
        % theta0_c
        function val=get.theta0_c(obj), val=atan2(obj.consec_mean(sin(obj.theta0_vec)),obj.consec_mean(cos(obj.theta0_vec))); end
        
        % theta
        function val=get.theta(obj)
            if isempty(obj.s_c) || isempty(obj.L_c), val=double.empty(0,1); return, end
            phi=2*pi*obj.s_c./obj.L_c;
            val=obj.theta0_c.*cos(phi)-obj.theta0_c.^3.*(obj.cf_c.*cos(3*phi) + obj.cs_c.*sin(3*phi));
        end
        
        % c
        function val=get.c(obj)
            if isempty(obj.theta) || isempty(obj.s_c), val=double.empty(0,1); return, end
            val=diff(obj.theta)./diff(obj.s_c);
        end
        
        % x
        function val=get.x(obj)
            if isempty(obj.theta) || isempty(obj.s), val=double.empty(0,1); return, end
            dx=cos(obj.theta).*diff(obj.s);
            val=cumsum([0;dx]);
        end
        
        % y
        function val=get.y(obj)
            if isempty(obj.theta) || isempty(obj.s), val=double.empty(0,1); return, end
            dy=sin(obj.theta).*diff(obj.s);
            val=cumsum([0;dy]);
        end
        
        %%% generic methods
        function varargout=plot(obj,varargin)
            % Plot Plots the Kinoshita curve
            %
            %   Plot(Kinoshita) plots
            allpl=reshape([reshape({obj(:).x},1,[]); reshape({obj(:).y},1,[])],[],1);
            hp=plot(allpl{:},varargin{:});
            if nargout==1
                varargout{1}=hp;
            end
            axis equal
        end
        
        function [polyout_ar shp] =polygon(objin,Win)
            % Polygon Returns a Polygon given a channel width
            %
            %   Polygon(K,W) returns a Polygon object with the
            %   banklines of a channel with the Kinoshita curve K as centreline
            %   and W as width. W can be a scalar or a vector with the same
            %   size as s.
            nout=numel(objin);
            polyout_ar(nout)=Polygon;
            polyout_ar=reshape(polyout_ar,size(objin));
            for cp=1:numel(objin)
                obj=objin(cp);
                polyout=polyout_ar(cp);
                W=Win;
                try
                    W=ones(size(obj.s)).*W;
                catch err
                    continue
                end
                ang=[obj.theta(1);obj.consec_mean(obj.theta);obj.theta(end)];
                Nx=-sin(ang);
                Ny=cos(ang);
                xlb=obj.x+W/2.*Nx;
                ylb=obj.y+W/2.*Ny;
                xrb=obj.x-W/2.*Nx;
                yrb=obj.y-W/2.*Ny;
                polyout.x=[xlb;flipud(xrb);xlb(1)];
                polyout.y=[ylb;flipud(yrb);ylb(1)];
		% individual boundaries
		shp = struct();
		shp(1).X = rvec(xlb);
		shp(1).Y = rvec(ylb); 
		shp(1).side = 'north';

		shp(2).X = [xrb(1), xlb(1)];
		shp(2).Y = [yrb(1), ylb(1)];
		shp(2).side = 'west';

		shp(3).X = fliplr(rvec(xrb));
		shp(3).Y = fliplr(rvec(yrb));
		shp(3).side = 'south';

		shp(4).X = [xlb(end), xrb(end)];
		shp(4).Y = [ylb(end), yrb(end)];
		shp(4).side = 'east';

            end
        end
        
    end % methods
    methods(Access=protected)
        function val=vectinput(obj,property)
            % make sure a property is the same size as s
            if isempty(obj.(property)) || isempty(obj.s), val = double.empty(0,1); return, end
            try
                val=bsxfun(@times,ones(size(obj.s)),obj.(property));
            catch err
                err2=MException('Kinoshita:InconsistenProperties','Property %s is inconsistent with size of s',property);
                addCause(err2,err);
                throw(err2);
            end
        end
    end
    methods(Static, Access=protected)
        function val=consec_mean(in)
            % Computes the mean of all consecutive elements in input
            if (numel(in)<2);
                val=double.empty(0,1); return
            end
            val=(in(1:end-1)+in(2:end))/2;
        end
    end
end
