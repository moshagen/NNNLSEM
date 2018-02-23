function [out] = gendata( phi, lambda, N, f, g, acc, xi)
%GENDATA Generates non-normally distributed random variables conforming to
%a given structural equation model
%   Detailed explanation goes here


y=cell(2,1);
y{1}='lat(:,jjj)';
y{2}='err(:,i)';


%Initialize
x=norminv(10^(-acc):10^(-acc):(1-10^(-acc)),0,1);
nLat = size(phi,2);
nMan = 0;
for iii=1:nLat
    nMan=nMan+size(lambda{iii},1);
end
decompLat=chol(phi);
vark=zeros(2,nMan);
cork=zeros(2,nMan);
hMan=zeros(length(x),nMan);
hErr=zeros(length(x),nMan);
corrMat=zeros(nMan,nMan);

%Check if error covariance is specified
if (~exist('xi', 'var'))
    xi = eye(nMan);
end


% h and hz store the functions that are needed to estimate the distortion (hz) 
% and to calculate the manifest variable (h)
hz=cell(2,1);
hz{1}=f; hz{2}=g;
h=cell(2,1);
h{1}=cell(nMan,1); h{2}=cell(nMan,1);


%Calculate true covariance matrix and latent correction matrix
trueCov=zeros(nMan);
i1=1; i2=1;
for iii=1:nLat
    for jjj=1:size(lambda{iii},1)
        for kkk=1:nLat
            for lll=1:size(lambda{kkk},1)
                trueCov(i1,i2)=...
                    lambda{iii}(jjj).*lambda{kkk}(lll).*phi(iii,kkk);
                corrMat(i1,i2)=phi(iii,kkk);
                i2=i2+1;
            end
        end
        i1=i1+1; i2=1;
    end
end
errsigma=diag(sqrt(1-diag(trueCov)));
trueCov=trueCov+errsigma*xi*errsigma;

%Calculate the variance and correlation correction factors
for jjj=1:2
    if (size(hz{jjj},1)==1)%change input in case function is the same for every variable
        k=hz{jjj};
        hz{jjj}=cell(nMan,1);
        for iii=1:nMan
            hz{jjj}{iii}=k;
        end       
    end
        for iii=1:size(hz{jjj},1)
            vark(jjj,iii)=1/sqrt(var(eval(hz{jjj}{iii}))); %correcting for variance
            cork(jjj,iii)=corr(x',eval(hz{jjj}{iii})'); %correcting for correlation
            h{jjj}{iii}=strrep(hz{jjj}{iii}, 'x', y{jjj}); 
        end
    
end

%Calculate the remaining deviation in the manifest variables and an error 
%covariance matrix correcting the deviation
i=1;
for jjj=1:nLat
    for iii=1:size(lambda{jjj},1)
        hMan(:,i)=eval(hz{1}{i})'.*(lambda{jjj}(iii)./cork(1,i)).*vark(1,i);
        hErr(:,i)=eval(hz{2}{i}).*vark(2,i).*sqrt(1-(lambda{jjj}(iii)./cork(1,i))^2);
        i=i+1;
    end
end
latDevi=cov(hMan).*corrMat;
errCovDiff=trueCov-latDevi;
errCorr=(errCovDiff./(cov(hErr)));
decompErr=chol(errCorr);

lat=randn(N,nLat)*decompLat; %Create latent variables
err=randn(N,nMan)*decompErr; %Create errors

%Calculate manifest variables
man=zeros(N,nMan);
i=1;
for jjj=1:nLat
    for iii=1:size(lambda{jjj},1)
        man(:,i)=...
            eval(h{1}{i}).*(lambda{jjj}(iii)./cork(1,i)).*vark(1,i)+...
            eval(h{2}{i}).*sqrt(1-(lambda{jjj}(iii)./cork(1,i))^2).*vark(2,i);
        i=i+1;
    end
end


out.trueCov=trueCov;
out.decomLat=decompLat;
out.decomErr=decompErr;
out.nLat=nLat;
out.nMani=nMan;
out.lat=lat;
out.err=err;
out.vark=vark;
out.cork=cork;
out.man=man;
end
