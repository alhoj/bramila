function cfg=bramila_glm_np(cfg)
% simple permutation based univariate parametric GLM for first level analysis
% Usage:
%   cfg = bramila_glm_np(cfg);
%
%   Input:
%   cfg.infile= path to nifti file
%   cfg.regressor = a 2D design matrix with number of rows equal to number of time
%       points, and number of columns equal to number of regressors. No
%       orthogonalization happens, i.e. regressors are tested independently
%   cfg.cdtP = p value of cluster defining threshold (default 0.05, uses t-distribution )
%   cft.cdtR = r value of cluster defining threshold (overrides cdtP)
%   cfg.NPERM = amount of permutations (default is 5000)
%   cfg.seed = random generator seed
%   cfg.toi = a vector of indeces with times of interest (set to 1:size(regressor,1) if you don't need it) 
%
%   cfg.CR =[0 0.5]; range to estimate correlations
%%
%   Output:
%   cfg.vol = one volume for each regressor, cluster corrected
%   cfg.cth = cluster threshold

% input testing - still missing

regressor=cfg.regressor(cfg.toi);
T=size(regressor,1);
nii=load_nii(cfg.infile);
mask=sign(var(nii.img,0,4));
inmask=find(mask>0);

% estimates the average autocorrelation of regressors
DF=bramila_autocorr(regressor,regressor);

AC=T./DF;
mAC=mean(AC);

cfg.AC=AC;
cfg.mAC=mAC;
cfg.DF=DF;

block_size=round(mAC);

% estimate cfg.cdtR
rho=cfg.CR(1):0.001:cfg.CR(2);

for rhoID=1:length(rho)
    tempp(rhoID)=pvalPearson('r', rho(rhoID), DF); % it's not the mean(AC) but the DF
end
cfg.cdtR=rho(min(find(tempp<cfg.cdtP)))


% creating permuted regressors
surrog=zeros(T,cfg.NPERM+1);


rng(cfg.seed);

for perms=1:cfg.NPERM+1
    if(perms==1)
        surrog(:,perms)=regressor;
    else
        surrog(:,perms)=bootres(regressor,0,block_size);
    end
end

cfg.surrog=surrog;

csurr=zeros(cfg.NPERM+1,1);
data=reshape(nii.img(:,:,:,cfg.toi),[],T);
data=data'; % Time in 1st dimension

for perms=1:cfg.NPERM+1
    disp(['Perm # ' num2str(perms)]);
    rsurrtemp=corr(data(:,inmask),squeeze(surrog(:,perms)));
    rsurr=zeros(size(mask));
    rsurr(inmask)=rsurrtemp;
    if(perms==1)
        vol=rsurr;
    end

    [x y z]=ind2sub(size(mask),find(rsurr>cfg.cdtR));
    a=[x y z];
    clust=spm_clusters(a',18);

    if(perms==1)
        rvol=rsurr;
        volclust=clust;
        vola=a;
    end

    uu=unique(clust); 
    uu(find(uu==0))=[];
    cs=histc(clust,uu); % cluster sizes
    csum=[];
    for clusti=1:length(uu)
        csum(clusti)=sum(rsurrtemp(find(clust==uu(clusti)))); % summed correlation coefficients within clusters
    end

    switch cfg.clusterstatistic
        case 'maxsize'
            if(isempty(max(cs)))
                csurr(perms)=0;
            else
                csurr(perms)=max(cs);
            end
        case 'maxsum'
            if(isempty(max(csum)))
                csurr(perms)=0;
            else
                csurr(perms)=max(csum);
            end
        otherwise
            error('cfg.clusterstatistic has to be either ''maxsize'' or ''maxsum''');
    end

end
%%
cfg.cth=prctile(csurr(2:end),100*(1-cfg.alpha));

newmask=zeros(size(mask));

for c = 1:length(volclust)
    cids=find(c==volclust);
    
    switch cfg.clusterstatistic
        case 'maxsize'
            csize=length(cids);
            if(csize<cfg.cth)
                for cc=1:length(cids)
                    newmask(vola(cids(cc),1),vola(cids(cc),2),vola(cids(cc),3))=0;
                end
            else
                for cc=1:length(cids)
                    newmask(vola(cids(cc),1),vola(cids(cc),2),vola(cids(cc),3))=c;
                end
            end
            
        case 'maxsum'   
            csum=0;
            for cc=1:length(cids)
                csum=csum+rvol(vola(cids(cc),1),vola(cids(cc),2),vola(cids(cc),3));
            end
            if(csum<cfg.cth)
                for cc=1:length(cids)
                    newmask(vola(cids(cc),1),vola(cids(cc),2),vola(cids(cc),3))=0;
                end
            else
                for cc=1:length(cids)
                    newmask(vola(cids(cc),1),vola(cids(cc),2),vola(cids(cc),3))=c;
                end
            end
    end


    
end

cfg.vol(:,:,:)=vol;
cfg.cmask(:,:,:)=newmask;

cfg.csurr=csurr;


%% helpers

function tsout=bootres(ts,type,MBS)
L=length(ts);
if(type==0)
    notok=1;
    tempids=(1:L)';
    while(notok==1)
        B=MBS+round(100*rand);	% block size
        R=ceil(L/B);
        temp=zeros(R*B,1);
        temp(1:L)=1:L;
        temp=circshift(temp,round(MBS/2)+round((length(temp)-round(MBS/2))*rand));
        mat=reshape(temp,B,R);
        mat=mat(:,randperm(R));
        temp=mat(:);
        
        tempids=temp(find(temp>0));
        if(min(abs(tempids-(1:L)'))>round(MBS/2))
            notok=0;
        end
    end
    tsout=ts(tempids);
else
    tsout=circshift(ts,round(MBS/2)+round((L-round(MBS/2))*rand));
end
