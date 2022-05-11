function [mixed, abf] = getSynDataPatial(varargin)

%
% Generate synthetic data.
% The spectra of simulated data is obtained from the USGS library "signatures"
%
% Input
%   - A: matrix of reflectances (endmember signatures)
%   - dim: size of image (dim/8==0)
%   - win: 0 - size of smoothing filter
%   - pure: no pure pixels, 'no'/'yes'
%   - purity: intensity of abundance (maximum purity:1)
%   - qpiror: number of weak signal
%   - wk: intensity of abundance of weak signal
%   - disp: visually display results
%
% Output
%   - mixed: generated synthetic mixed data
%   - abf: actual abundance fractions
%
% The pure pixels can be removed by adding the following two lines
%        ----Index = ceil(find(abf>0.8)/c);
%        ----abf(:,Index) = 1/c*ones(c,1)*ones(1,length(Index));
%
% This matlab script is modified by Xiangfei Shen (xfshen95@outlook.com)
% for generating synthetic datasets containing weak signals
% Example:
% [mixed, abf] = getSynDataPatial('signatures',A,'dim', 64,'win',7,'ispure', ...
%    'no','purity',0.8,'priors',2,'weakness',0.15,'isdisplay','yes');

%%
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'SIGNATURES',    A=varargin{i+1};
            case 'DIM',    dim=varargin{i+1};
            case 'WIN',    win=varargin{i+1};
            case 'ISPURE',      pure=varargin{i+1};
            case 'PURITY',      purity=varargin{i+1};
            case 'PRIORS',    qprior=varargin{i+1};
            case 'WEAKNESS',    eta=varargin{i+1};
            case 'ISDISPLAY',    disp=varargin{i+1};
            otherwise
                error(['Unrecognized option: ',varargin{i}]);
        end
    end
end



[band, p] = size(A);
p=p-qprior;
av=1/p;

%dim = 64;
%dim = 160;
label = ones((dim/8)^2,1);
num = floor(length(label)/p);

for i=1:p-1
    label((i-1)*num+1:i*num) = (i+1);
end

ridx = randperm(length(label));
label = label(ridx)';
label = reshape(label,dim/8,dim/8);
abf = zeros(dim,dim,p);
img = zeros(dim,dim);
for i=1:dim
    for j=1:dim
        for cls = 1:p
            if label(floor((i-1)/8)+1,floor((j-1)/8)+1) == cls
                tmp = zeros(p,1);
                tmp(cls) = 1;
                abf(i,j,:) = tmp;
                img(i,j) = p;
            end
        end
    end
end

%low pass filter
H = ones(win,win)/(win*win);
img_fil = filter2(H,img);
for i=1:p
    abf(:,:,i) = filter2(H,abf(:,:,i));
end
abf = abf(ceil(win/2):end-floor(win/2),ceil(win/2):end-floor(win/2),:);


% generate mixtures
[M,N,p] = size(abf);
abf = reshape(abf,M*N,p)';

% remove pure pixels
if strcmp(pure,'no')
    if purity~=1
        Index = ceil(find(abf>purity)/p);
        abf(:,Index) = 1/p*ones(p,1)*ones(1,length(Index));
        
    end
end

wl=repmat(zeros(1,M*N),qprior,1);


% generate abundance maps for weak signals 
if qprior>p
    disp('Q should be less than P in this dataset!!')
elseif qprior>=1
    
    for ii=1:qprior
        
        wl(ii,abf(ii,:)==av)=eta;
        abf(:,abf(ii,:)==av)=(1-eta)/p;
        av=(1-eta)/p; 
    end
    
end

abf=[abf;wl];

if strcmp(disp,'yes')
    
    Y=(A*abf);
    
    figure
    plot(Y(50,:),Y(100,:),'ko')
    xlabel('Band 100');
    ylabel('Band 50');
    hold on
    h=plot(A(50,1:p),A(100,1:p),'ro');
    set(h,'MarkerFaceColor',get(h,'color'));
    hh=plot(A(50,p+1:end),A(100,p+1:end),'go');
    set(hh,'MarkerFaceColor',get(hh,'color'));
    
    figure;
    for j=1:p+qprior;
        subplot_tight(1, p+qprior, j,[.01 .01]);
        imagesc(reshape(abf(j,:)',M, N),[0,1]);axis image;axis off;
    end 
  
end
mixed = reshape((A*abf)',M,N,band);
end
function h=subplot_tight(m,n,p,margins,varargin)

if (nargin<4) || isempty(margins)
    margins=[0.01,0.01]; % default margins value- 1% of figure
end

if length(margins)==1
    margins(2)=margins;
end

%note n and m are switched as Matlab indexing is column-wise, while subplot indexing is row-wise :(
[subplot_col,subplot_row]=ind2sub([n,m],p);


height=(1-(m+1)*margins(1))/m; % single subplot height
width=(1-(n+1)*margins(2))/n;  % single subplot width

% note subplot suppors vector p inputs- so a merged subplot of higher dimentions will be created
subplot_cols=1+max(subplot_col)-min(subplot_col); % number of column elements in merged subplot
subplot_rows=1+max(subplot_row)-min(subplot_row); % number of row elements in merged subplot

merged_height=subplot_rows*( height+margins(1) )- margins(1);   % merged subplot height
merged_width= subplot_cols*( width +margins(2) )- margins(2);   % merged subplot width

merged_bottom=(m-max(subplot_row))*(height+margins(1)) +margins(1); % merged subplot bottom position
merged_left=min(subplot_col)*(width+margins(2))-width;              % merged subplot left position
pos_vec=[merged_left merged_bottom merged_width merged_height];

% h_subplot=subplot(m,n,p,varargin{:},'Position',pos_vec);
% Above line doesn't work as subplot tends to ignore 'position' when same mnp is utilized
h_subplot=subplot('Position',pos_vec,varargin{:});

if nargout~=0
    h=h_subplot;
end
end
