% Towards Weak Signal Analysis in Hyperspectral Data: A Semi-supervised Unmixing Perspective
%
% Usage: [W,H,results,his] = hyperweak(V,varargin);
%==========================================================================
% Input:  		   V  -  L*N hyperspectrtal data;
%             W_INIT  -  L*p initial endmember matrix;
%             H_INIT  -  p*N initial abundance matrix;
%             LAMBDA  -  sparsity regularization parameter;
%                TAU  -  regularization parameter for controling weak signal degradation;
%             PRIORS  -  number of priors
%                TOL  -  stop condition 1: tolerance for stoping optimization
%            MAXITER  -  stop condition 2: maximum number of iterations
%               SIZE  -  image size: [nr,nc,nb]
%            VERBOSE  -  display current information

%
% Output: W  -  L*p estimated endmember matrix
%         H  -  p*N estimated abundance matrix
%   results  -  record optimization information
%       his  -  record W in each iteration
%==========================================================================
%
% Copyright (C) 2021, Xiangfei Shen (xfshen95@outlook.com)
%                     Haijun Liu
%					  Jian Qin
%					  Fangyuan Ge
%					  Xichuan Zhou* (zxc@cqu.edu.cn)
%                     School of Microelectronics and Communication Engineering, Chongqing University£¬
%                     All rights reserved.
%
%
% References:
% Submitted to IEEE TGRS (Under Review)
%
%
% Last Modified:
% 10 January, 2022 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [W,H,results,his] = hyperweak(V,varargin)

% test 4 mur
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'W_INIT',    W=varargin{i+1};
            case 'H_INIT',    H=varargin{i+1};
            case 'LAMBDA',    lambda=varargin{i+1};
            case 'TAU',      tau=varargin{i+1};
            case 'PRIORS',      q=varargin{i+1};
            case 'TOL',      tol=varargin{i+1};
            case 'SIZE',      imgsize=varargin{i+1};
            case 'MAXITER',    maxiter=varargin{i+1};
            case 'VERBOSE',     verbose=varargin{i+1};
            otherwise
                error(['Unrecognized option: ',varargin{i}]);
        end
    end
end

his=cell(1);
%%
% calculate initial gradient
gradW = W*(H*H') - V*H'; gradH = (W'*W)*H - W'*V;
initgrad = norm([gradW; gradH'],'fro');
if strcmp(verbose,'on')
    fprintf('Init gradient norm %f\n', initgrad);
end

%tolerance for a relative stopping condition
tolW = max(0.001,tol)*initgrad;
tolH = tolW;
objhistory =  sum(sum((V-W*H).^2));
objhistory = [objhistory 0];
inc = 0;
inc0 = 0;

%initial obj
results=[];
[L,p]=size(W);
AA=[zeros(p-q,q);eye(q,q)];%auxiliary matrix
WWh=W(:,p-q+1:end);

for iter=1:maxiter
    % stopping condition
    Wup=W;
    Hup=H;
    objdiff=objhistory(end-1)-objhistory(end);
    if objdiff>0.0001
        inc = 0;
    else
        disp('inc');
        inc = inc+1;
        inc0 = inc0+1;
    end
    if iter < 5, inc = 0; end
    if inc >= 5 && inc0 >= 10, break; end
    
    %% --------------Updata Weights ---------------------------

    WW=RWSpasity(H,imgsize);
    
    %%  ------------- Update W -------------------------------
    
    [W,iterW,gradW] = subprob_W_OG(V',H',W',WWh,AA,tau,10,200,tolW); W = W'; gradW = gradW';
    if iterW==1, tolW = 0.1 * tolW; end
    
    
    %%  ------------- Update H -------------------------------
    tW = ASC(W,q); tV = ASC(V,q); % ASC
    
    [H,iterH,gradH]=subprob_H_OG(tV,tW,H,lambda,WW,10,200,tolH);
    if iterH == 1, tolH = 0.1 * tolH; end
 
    %%  ------------- Results Display ------------------------
    projnorm = norm([gradW(gradW<0 | W>0); gradH(gradH<0 | H>0)],'fro');
    if strcmp(verbose,'on')
        fprintf('\nIter = %d [obj:%.4f --objdiff:%.4f-- proj-grad norm %.4f]\n', iter, sum(sum((V-W*H).^2)),objdiff,projnorm);
    end
    objhistory = [objhistory sum(sum((V-W*H).^2))];
    % results recording
    results(iter,1) = sum(sum((V-W*H).^2));
    results(iter,2) = objdiff;
    results(iter,3) = sum(sum((W-Wup).^2));
    results(iter,4) = sum(sum((H-Hup).^2));
    results(iter,5) = norm([gradW(gradW<0 | W>0); gradH(gradH<0 | H>0)],'fro');
    his{iter}=W;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [MA] = ASC(M,thresh)
[L,P] = size(M);
MA = thresh*ones(L+1,P);
MA(1:L,:) = M;
end

%%

%%
function W = RWSpasity(S,imgsize)
nr=imgsize(1);nc=imgsize(2);

[p,N]=size(S);

S=reshape(S',nr,nc,p);
for i=1:p;
    
    S(:,:,i) = wiener2(S(:,:,i),[3 3]);
    
end

S=reshape(S,nr*nc,p)';
W=1./(abs(S)+eps);


end

%%




function [H,iter,Grad]=subprob_H_OG(V,W,Z,beta,WW,iterMin,iterMax,tol)

STOP_RULE = 1;
WtV = W'*V; WtW = W'*W;

if ~issparse(WtW),
    L=norm(WtW);	% Lipschitz constant
else
    L=norm(full(WtW));
end
H=Z;    % Initialization
Grad=WtW*Z-WtV+beta*WW;     % Gradient
alpha1=1;

for iter=1:iterMax,
    H0=H;
    H=max(Z-Grad/L,0);    % Calculate sequence 'Y'
    alpha2=0.5*(1+sqrt(1+4*alpha1^2));
    Z=H+((alpha1-1)/alpha2)*(H-H0);
    alpha1=alpha2;
    Grad=WtW*Z-WtV+beta*WW;
    
    % Stopping criteria
    if iter>=iterMin,
        pgn=GetStopCriterion(STOP_RULE,Z,Grad);
        if pgn<=tol,
            break;
        end
    end
end

Grad=WtW*H-WtV+beta*WW;



end

function [H,iter,Grad]=subprob_W_OG(V,W,Z,WWh,AA,tau,iterMin,iterMax,tol)

STOP_RULE = 1;
WtV = W'*V; WtW = W'*W;

if ~issparse(WtW),
    L=norm(WtW);	% Lipschitz constant
else
    L=norm(full(WtW));
end
H=Z;    % Initialization
Grad=WtW*Z-WtV-tau*(WWh*AA'-Z'*(AA*AA'))';     % Gradient
alpha1=1;

for iter=1:iterMax,
    H0=H;
    H=max(Z-Grad/L,0);    % Calculate sequence 'Y'
    alpha2=0.5*(1+sqrt(1+4*alpha1^2));
    Z=H+((alpha1-1)/alpha2)*(H-H0);
    alpha1=alpha2;
    Grad=WtW*Z-WtV-tau*(WWh*AA'-Z'*(AA*AA'))';
    
    % Stopping criteria
    if iter>=iterMin,
        pgn=GetStopCriterion(STOP_RULE,Z,Grad);
        if pgn<=tol,
            break;
        end
    end
end

Grad=WtW*H-WtV-tau*(WWh*AA'-Z'*(AA*AA'))';


end

