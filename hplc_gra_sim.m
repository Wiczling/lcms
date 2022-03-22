function  samples2 = hplc_gra_sim(samples,datastruct,hplcparam_design)

nObs = size(hplcparam_design,1);
nSamples= size(samples.logkwHat,1);
nAnalytes=datastruct.nAnalytes;

hplcparam=hplcparam_design;
mod=hplcparam(:,7);
steps=hplcparam(:,13); 

nDiss=datastruct.R;
chargesA=datastruct.chargesA;
chargesB=datastruct.chargesB;
 
samples2.trHatCond=zeros(nSamples, nAnalytes, nObs); 
samples2.trHatPred=zeros(nSamples, nAnalytes, nObs);
samples2.trObsCond=zeros(nSamples, nAnalytes, nObs); 
samples2.trObsPred=zeros(nSamples, nAnalytes, nObs);
 
for z = 1:nSamples
  z
dlogkT=samples.dlogkT(z,:);
logkwx=squeeze(samples.logkwx(z,:,:));
S1mx=squeeze(samples.S1mx(z,:,:));
S1ax=squeeze(samples.S1ax(z,:,:));
pKaw =squeeze(samples.pKaw(z,:,:));
alpham =squeeze(samples.alpham(z,:,:));
alphaa =squeeze(samples.alphaa(z,:,:));
S2mHat =samples.S2mHat(z,:);   
S2aHat =samples.S2aHat(z,:);  
nu = 3;
apH =squeeze(samples.apH(z,:));
sigma =samples.sigma(z,:);

dlogkTPred=samples.dlogkTPred(z,:);
logkwxPred=squeeze(samples.logkwxPred(z,:,:));
S1mxPred=squeeze(samples.S1mxPred(z,:,:));
S1axPred=squeeze(samples.S1axPred(z,:,:));
pKawPred=squeeze(samples.pKawPred(z,:,:));
alphamPred=squeeze(samples.alphamPred(z,:,:));
alphaaPred=squeeze(samples.alphaaPred(z,:,:));
sigmaPred =samples.sigmaPred(z,:);

trHatCond=zeros(nAnalytes,nObs);
trHatPred=zeros(nAnalytes,nObs);

 for j=1:nAnalytes
 for i=1:nObs

 if (mod(i)==1)
  y_hat_Cond = chromgratrapz(steps(i), ...
               logkwx(j,:) + dlogkT(j)*hplcparam(i,11),  ...
               S1mx(j,:),  ...
               pKaw(j,:),  ...
               alpham(j,:), ...
               S2mHat, ...
               apH,...
               chargesA(j,:), ...
               chargesB(j,:), ...
               nDiss(j), ...
               hplcparam(i,:));
 end

 if (mod(i)==2)
  y_hat_Cond = chromgratrapz(steps(i),  ...
               logkwx(j,:) + dlogkT(j)*hplcparam(i,11),  ...
               S1ax(j,:),  ...
               pKaw(j,:),  ...
               alphaa(j,:), ...
               S2aHat, ...
               apH,...
               chargesA(j,:), ...
               chargesB(j,:), ...
               nDiss(j), ...
               hplcparam(i,:));
 end

  trHatCond(j,i) = hplcparam(i,3) + hplcparam(i,4) + y_hat_Cond; 
  
 end

 for i=1:nObs
 if (mod(i)==1)
  y_hat_Pred = chromgratrapz(steps(i), ...
               logkwxPred(j,:) + dlogkTPred(j)*hplcparam(i,11), ...
               S1mxPred(j,:),  ... 
               pKawPred(j,:), ...
               alphamPred(j,:), ...
               S2mHat, ...
               apH,...
               chargesA(j,:), ...
               chargesB(j,:), ...
               nDiss(j), ...
               hplcparam(i,:));
 end

 if (mod(i)==2)
  y_hat_Pred = chromgratrapz(steps(i),  ...
               logkwxPred(j,:) + dlogkTPred(j)*hplcparam(i,11), ...
               S1axPred(j,:),  ...
               pKawPred(j,:), ...
               alphaaPred(j,:), ...
               S2aHat, ...
               apH,...
               chargesA(j,:), ...
               chargesB(j,:), ...
               nDiss(j), ...
               hplcparam(i,:));
 end
 
 trHatPred(j,i) = hplcparam(i,3) + hplcparam(i,4) + y_hat_Pred;
 end

samples2.trHatCond(z,:,:)=trHatCond;
samples2.trHatPred(z,:,:)=trHatPred;
samples2.trObsCond(z,:,:)=trHatCond+(sigma(:)*ones(1,nObs)).*trnd(nu,size(trHatCond));
samples2.trObsPred(z,:,:)=trHatPred+(sigmaPred(:)*ones(1,nObs)).*trnd(nu,size(trHatPred));
 
 end
end

function sol = gra_state(t, hplcparam) 

tg = hplcparam(1);
td = hplcparam(2);
fio = hplcparam(5);
fik = hplcparam(6);
pHo = hplcparam(8);
alpha1 = hplcparam(9);
alpha2 = hplcparam(10);

fi = fio+(fik-fio)/tg*(t-td);

if (t<td)
fi = fio;
end
if (t>tg+td)
fi = fik;
end

sol(1)=fi;
sol(2)=pHo+alpha1*fi+alpha2*fi.^2;

function logki=funlogki(logkwx, S1, pKaw, alpha, S2, apH, nDiss, chargesA, chargesB, fipH)

fi = fipH(1);
pH = fipH(2);

logkix   = logkwx-S1*fi/(1+S2*fi)+ chargesA*apH(1)*(pH-7) + chargesB*apH(2)*(pH-7);
pHmpKa    = pH-(pKaw+alpha*fi);

logki = NaN;

if (nDiss==0)
  logki = logkix(1); 
end
if (nDiss==1)
  logki = logkix(1) + ...
          log1p_exp(log(10)*(pHmpKa(1)+logkix(2)-logkix(1)))/log(10)- ...
          log1p_exp(log(10)*(pHmpKa(1)))/log(10);
end
if (nDiss==2)
  logki = logkix(1) + ...
          log1p_exp(log(10)*(pHmpKa(1)+logkix(2)-logkix(1)) +  ...
          log1p_exp(log(10)*(pHmpKa(2)+logkix(3)-logkix(2))))/log(10)- ...
          log1p_exp(log(10)*(pHmpKa(1)) +  ...
          log1p_exp(log(10)*(pHmpKa(2))))/log(10);
end

function cki_b = areaandslope(time1, time2, invki1, invki2)

if (invki2>1.001*invki1)
bo = (log(invki2)-log(invki1))/(time2-time1);
cki = (invki2-invki1)/bo;
else
bo  = 0.001/(time2-time1);
cki = (time2-time1)*(invki2+invki1)/2;
end

cki_b(1) = cki;
cki_b(2) = bo;

function tr = chromgratrapz(steps, logkwx, logkmx, pKaw, alpha, S2mHat, apH, chargesA, chargesB, nDiss, hplcparam)

tg = hplcparam(1);
td = hplcparam(2);
to = hplcparam(3);

dt = tg/steps;

time1 = 0;
time2 = td;

fipH1 = gra_state(time1,  hplcparam);
% fipH2 = fipH1;

logki1 = funlogki(logkwx, ...
                  logkmx, ... 
                  pKaw, ...
                  alpha,...
                  S2mHat,...
                  apH, ...
                  nDiss,...
                  chargesA, chargesB, ... 
                  fipH1);

% logki2 = logki1;

invki1 = 1/to/10^logki1;
invki2 = invki1;

cumki1 = 0;
cumki2 = td*invki1;
bo     = 0.001/td;  

for x=1:steps 
   if (cumki2<1) 
    time1 = time2;
    time2 = time2+dt;
%     fipH1 = fipH2;
    fipH2 = gra_state(time2,  hplcparam);
%     logki1 = logki2;
    logki2 = funlogki(logkwx, ...
                  logkmx,  ...
                  pKaw, ...
                  alpha,...
                  S2mHat,...
                  apH, ...
                  nDiss,...
                  chargesA, chargesB, ... 
                  fipH2);
    invki1 = invki2;
    invki2 = 1/to/10^logki2;
    cki_b = areaandslope(time1, time2, invki1, invki2);
    cumki1 = cumki2;
    cumki2 = cumki2 + cki_b(1); 
    bo      = cki_b(2); 
   end
end

if (cumki2>=1) 
    tr = time1+log1p((1-cumki1)*bo/invki1)/bo;
elseif (cumki2<1)
    tr = time2 + (1-cumki2)/invki2;
else
end

function y = log1p_exp(x)
y=log(1+exp(x));

function y = log1p(x)
y=log(1+x);