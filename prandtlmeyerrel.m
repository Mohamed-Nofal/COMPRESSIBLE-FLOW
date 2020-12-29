function varargout=prandtlmeyerrel(varargin)
if nargin==0
    help prandtlmeyerrel
    varargout={};
    return
elseif nargin==2
    gamma=varargin{2};
    choice='M';
elseif nargin==1
    gamma=1.4; %assume air
    choice='M';
elseif nargin==3
    choice=varargin{3};
    if ~isempty(varargin{2})
        gamma=varargin{2};
    else
        gamma=1.4;
    end
else
    error('Inputs not accepted')
end
if ~isnumeric(varargin{1}) || ~isnumeric(gamma) || ~any(strcmpi(choice,{'M','P','A'}))
    error('Inputs not accepted')
end
%%%%%%%%%%%%%%%%%%%%SOLVE FOR MISSING DATA
if strcmpi(choice,'M') && all(varargin{1}>=1)
    M=reshape(varargin{1},numel(varargin{1}),1);
    mu=asin(1./M ).*180./pi;  %eq (4.1)
    up=(sqrt((gamma+1)/(gamma-1)).*atan(sqrt(((gamma-1)/(gamma+1)).*(M.^2-1)))-atan(sqrt(M.^2-1))).*180./pi;    %eq 4.44
elseif strcmpi(choice,'P') && all(varargin{1}>=0)
    [~,p,~]=prandtlmeyerrel(1e100,gamma,'M');
    if all(varargin{1}<p)
        up=reshape(varargin{1},numel(varargin{1}),1);
        for i=length(up):-1:1    %solve for corresponding M and then use that M to find the rest of the values
            M(i,1)=fzero(@(x) (sqrt((gamma+1)/(gamma-1)).*atan(sqrt(((gamma-1)/(gamma+1)).*(x.^2-1)))-atan(sqrt(x.^2-1))).*180./pi - up(i),1);
        end
        mu=asin(1./M ).*180./pi;  %eq (4.1)
    else
        error('Input Out of Range')
    end
elseif strcmpi(choice,'A') && all(varargin{1}>=0) && all(varargin{1}<=90)
    mu=reshape(varargin{1},numel(varargin{1}),1);
    M=1./sind(mu);
    up=(sqrt((gamma+1)/(gamma-1)).*atan(sqrt(((gamma-1)/(gamma+1)).*(M.^2-1)))-atan(sqrt(M.^2-1))).*180./pi;    %eq 4.44
else
    error('Input Out of Range')
end
%%%%%%%%%%%%%%%%%%%%FORMAT OUTPUTS
if nargout<=1 %work with it if they dont wana differentiate
    varargout{1}=[M,up,mu];
elseif nargout==3 %put it back how you found it if they give enough output info
    varargout{1}=reshape(M,size(varargin{1}));
    varargout{2}=reshape(up,size(varargin{1}));
    varargout{3}=reshape(mu,size(varargin{1}));
else %probably a mistake
    error('Innaproiate Number of Output Arguements')
end
end