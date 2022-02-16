
% demo code for submitted paper titiled "Towards Weak Signal Analysis 
% in Hyperspectral Data: A Semi-supervised Unmixing Perspective"
%%
clear
close all

load jasperRidge2_R198_3ws

p=7;
q=3;
S=H;
A=W;

%% --------------------HyperWeak-----------------------------------
        
[What,Hhat,results] = hyperweak(X,...% L*N hyperspectral data
    'W_INIT',Minit,...% initilized endmembers
    'H_INIT',Rinit,...% initilized abundances
    'LAMBDA',1e-3,...% sparse regularization parameter
    'tau',500,...% parameter for prior knowledge degradation
    'priors',q,...% number of prior knowledge
    'size',[nr,nc],...% size of hyperspectral data in 2d perspective
    'TOL',0.05,...% tolerance
    'MAXITER',500,...% number of maximum iterations
    'VERBOSE','on');% display results
   
%% scaling        
        perm = permute_corr(A,What);
        What = What * perm;
        Hhat = (Hhat' * perm);
        What = What./repmat(max(What), size(What,1), 1);
        What = What.*repmat(max(A), size(What,1), 1);
        Hhat   =  Hhat./max(repmat(max(Hhat), size(Hhat,1), 1),eps);
        Hhat   = (Hhat.*repmat(max((S')), size(Hhat,1), 1))';
        
        Minitp = Minit * perm;
        Rinitp = (Rinit' * perm);
        Minitp = Minitp./repmat(max(Minitp), size(Minitp,1), 1);
        Minitp = Minitp.*repmat(max(A), size(Minitp,1), 1);
        Rinitp   =  Rinitp./max(repmat(max(Rinitp), size(Rinitp,1), 1),eps);
        Rinitp   = (Rinitp.*repmat(max((S')), size(Rinitp,1), 1))';
        
        sad=valSAD(A, What);
        [rmse,mrmse]=valRmse(S, Hhat);
        paraSam=sad(3,1:end-1);%sad values of total objects
        paraRmse=mrmse;%mrmse values 
        paraRmseq=rmse(p-q+1:end);%mrmse values of weak signals
        
        Data.What=What;Data.Hat=Hhat;Data.Sam=paraSam;Data.Rmse=paraRmse;Data.Rmseq=paraRmseq;
        
 
%% display
% signature comparisons among gt,roubst-osp, and hyperweak
% reflectance values are scaled
for ii=1:p
    
    figure;
    title('Method')
    plot(A(:,ii),'k-','LineWidth',2.25)
    hold on;
    plot(What(:,ii),'r-.','LineWidth',1.75)
    plot(Minitp(:,ii),'g-','LineWidth',1.25)
    set(gca,'FontSize',12);
    set(gcf,'unit','normalized','position',[0.4,0.4,0.25,0.3]);
    xlabel('Band')
    ylabel('Reflectance')
    legend('Library','HyperWeak','Robust-OSP','Location','best')
    %axis([wlen(1) wlen(end) 0 1]);
    grid off
    box on
    
end


% objective values
figure
plot(results(:,1),'ro-')
handle=legend('$\left \| Y-WH \right \|_{F}^{2}$');
set(handle,'Interpreter','latex', 'FontSize',12);
xlabel('Iteration');
ylabel('Values');

% residual values with regard to W
figure
plot(results(:,3),'ro-')
handle=legend('$\left \|W_{t}- W_{t-1} \right \|_{F}^{2}$');
set(handle,'Interpreter','latex', 'FontSize',12);
xlabel('Iteration');
ylabel('Values');

% residual values with regard to H
figure
plot(results(:,4),'gs-')
handle=legend('$\left \|H_{t}- H_{t-1} \right \|_{F}^{2}$');
set(handle,'Interpreter','latex', 'FontSize',12);
xlabel('Iteration');
ylabel('Values');



dispResults=[paraSam';mean(paraSam,2);mean(paraRmse);mean(paraRmseq)];
open dispResults

% Plot signatures (correspond to wavelength)
load jasperBands
load waveJasper
WW=zeros(224,7); 
WW(jbands,:)=W;
WW(jband_nochoose,:)=NaN;
W0=WW; 

Wh=zeros(224,7);
Wh(jbands,:)=What;
Wh(jband_nochoose,:)=NaN;
Ae=Wh;
wave=waveRef/1000;
for ii=1:p
    
    figure;
    plot(wave,W0(:,ii),'k-','LineWidth',2.25)
    hold on;
    plot(wave,Ae(:,ii),'r-.','LineWidth',1.75)
    set(gca,'FontSize',14);
    set(gcf,'unit','normalized','position',[0.4,0.4,0.25,0.3]);
    xlabel('Wavelength(\mum)')
    ylabel('Reflectance')
    legend('Library','HyperWeak','Location','best')
    axis([wave(1) wave(end) 0 1]);
    grid on
    box on
end


figure;
for j=1:p
    subplot_tight(2, p, j,[.01 .01]);
    imagesc(reshape(S(j,:)',nr, nc),[0,1]);axis image;axis off;
    subplot_tight(2, p, j+p,[.01 .01]);
    imagesc(reshape(Hhat(j,:)',nr, nc),[0,1]);axis image;axis off;
end
colormap(jet)


