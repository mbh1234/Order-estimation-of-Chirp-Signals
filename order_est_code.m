
% MODEL ORDER ESTIMATION: AIC,AICC, BIC, Corrected BIC; MAP, WIC, PAL, Hannan-Quinn
% SINUSOIDAL MODEL SIMULATIONS
clear;
a(1)=4.5;a(2)=3.3;
b(1)=3.6;b(2)=2.4;
om(1)=.3; om(2)=0.35;  
freqd=om;
sz=size(a);
numpar=sz(1,2);

maxcomp=8; % maximum number of components

for numsamp =50;
 'Sample Size',numsamp
sigsq=2;
stdev=sqrt(sigsq); 
df=2; 

stdev1=sqrt(1); stdev2=sqrt(4);  % good result
% Number of simulations
nsim=10;

%Initialization of counters

for comp=1:maxcomp
   bic(comp)=0;
   bicc(comp)=0;
   pal(comp)=0;
   map(comp)=0;
   aic(comp)=0;
   aicc(comp)=0;
   wic(comp)=0;
   hq(comp)=0;
end


for isim=1:nsim,
    'simulation #'; isim

%Data generation

%independent normal error
y=data_nor(numsamp,numpar,a,b,om,stdev);      


%independent t error
%y=data_t(numsamp,a,b,om,df,numpar);      

% normal mixture error
%y=data_n_mix(numsamp,numpar,a,b,om,stdev1,stdev2); 


%Estimation of sigma from full model
yadj=y;

gridnum=1000; % **** number of grid points for periodogram maximizer 

for ks=1:maxcomp
   inits=per_maxm(yadj,numsamp,gridnum); %usual periodogram initialization
   seq_ao(ks)=inits(1);
   seq_bo(ks)=inits(2);
   seq_omo(ks)=inits(3);
for ns=1:numsamp
    yadj(ns)=(yadj(ns)-inits(1)*cos(inits(3)*ns)-inits(2)*sin(inits(3)*ns));
end

end

initial=[seq_ao',seq_bo',seq_omo'];
L2_full=fminsearch('obj_L2_fun',initial,[],y,numsamp,maxcomp);

sigsq_null=0;
sigsq_full=0;
for k=1:numsamp
    yp(k)=0.0;
    for kp=1:maxcomp
        yp(k)=yp(k)+(L2_full(kp,1)*cos(L2_full(kp,3)*k)+L2_full(kp,2)*sin(L2_full(kp,3)*k));
    end
sigsq_null=sigsq_null+y(k)^2;    
diff_full_sq(k)=(y(k)-yp(k))^2;
sigsq_full=sigsq_full+diff_full_sq(k);
end
sigsq_null=sigsq_null/numsamp;
sigsq_full=sigsq_full/numsamp;


   
% MODEL SELECTION USING USUAL NONROBUST METHODS
for comp=1:maxcomp

nummax=comp;
yadj=y;

%gridnum=1000; % number of grid points for periodogram maximizer 

for ks=1:nummax
%   inits=per_max(yadj,numsamp,gridnum);
   seq_a(ks)=seq_ao(ks);%inits(1);
   seq_b(ks)=seq_bo(ks);%inits(2);
   seq_om(ks)=seq_omo(ks);%inits(3);
for ns=1:numsamp
    yadj(ns)=(yadj(ns)-seq_om(ks)*cos(seq_omo(ks)*ns)-seq_bo(ks)*sin(seq_omo(ks)*ns));
end
end

%comp
initial=[seq_a',seq_b',seq_om'];
L2=fminsearch('obj_L2_fun',initial,[],y,numsamp,comp);
clear seq_a;
clear seq_b;
clear seq_om;
clear initial;


sum_sq=0.0;
for k=1:numsamp
    yp(k)=0.0;
    for kp=1:comp
        yp(k)=yp(k)+(L2(kp,1)*cos(L2(kp,3)*k)+L2(kp,2)*sin(L2(kp,3)*k));
    end
        diff(k)=(y(k)-yp(k));
        sum_sq=sum_sq+diff(k)^2;

end
est_var=sum_sq/numsamp;
pal_est_sigsq(comp)=est_var;
if comp==1
    pal_rn=numsamp*log(sigsq_null/sigsq_null);
    pal_rhon=numsamp*log(sigsq_null/sigsq_full);
else
    pal_rn=numsamp*log(sigsq_null/pal_est_sigsq(comp-1));
    pal_rhon=numsamp*log(pal_est_sigsq(comp-1)/sigsq_full);
end

bic_ic_obj(comp)=numsamp*log(est_var)+(log(numsamp))*(3*comp);
bic_cor_obj(comp)=numsamp*log(est_var)+(log(numsamp))*(5*comp);
pal_obj(comp)=numsamp*log(est_var)+(3*comp)*log(3*maxcomp)*(log((pal_rn)+1)/log((pal_rhon)+1));

map_obj(comp)=numsamp*log(est_var)+10*comp*log(numsamp);
aic_obj(comp)=numsamp*log(est_var)+2*(3*comp);
aicc_obj(comp)=numsamp*log(est_var)+2*(3*comp*numsamp)/(numsamp-3*comp-1);
ca=2*(3*comp*numsamp)/(numsamp-3*comp-1); cb=(log(numsamp))*(5*comp);  % check the expression
wic_obj(comp)=((ca)/(ca+cb))*aicc_obj(comp)+((cb)/(ca+cb))*bic_cor_obj(comp);
hq_obj(comp)=numsamp*log(est_var)+2*(3*comp)*log(log(numsamp));

end % loop for all possible number of components

% USUAL BIC
[CR,I]=min(bic_ic_obj);
estcomp=I;
'USUAL BIC', estcomp


    if estcomp == 1
        bic(1)=bic(1)+1;
    elseif estcomp == 2
        bic(2)=bic(2)+1;
    elseif estcomp == 3
        bic(3)=bic(3)+1;
    elseif estcomp == 4
        bic(4)=bic(4)+1;
    elseif estcomp == 5
        bic(5)=bic(5)+1;
    elseif estcomp == 6
        bic(6)=bic(6)+1;
    elseif estcomp == 7
        bic(7)=bic(7)+1;
    else estcomp == 8
        bic(8)=bic(8)+1;

    end

% CORRECTED BIC
[CR,I]=min(bic_cor_obj);
estcomp=I;
'CORRECTED BIC', estcomp


    if estcomp == 1
        bicc(1)=bicc(1)+1;
    elseif estcomp == 2
        bicc(2)=bicc(2)+1;
    elseif estcomp == 3
        bicc(3)=bicc(3)+1;
    elseif estcomp == 4
        bicc(4)=bicc(4)+1;
    elseif estcomp == 5
        bicc(5)=bicc(5)+1;
    elseif estcomp == 6
        bicc(6)=bicc(6)+1;
    elseif estcomp == 7
        bicc(7)=bicc(7)+1;
    else estcomp == 8
        bicc(8)=bicc(8)+1;

    end

% PAL
[CR,I]=min(pal_obj);
estcomp=I;
'PAL', estcomp


    if estcomp == 1
        pal(1)=pal(1)+1;
    elseif estcomp == 2
        pal(2)=pal(2)+1;
    elseif estcomp == 3
        pal(3)=pal(3)+1;
    elseif estcomp == 4
        pal(4)=pal(4)+1;
    elseif estcomp == 5
        pal(5)=pal(5)+1;
    elseif estcomp == 6
        pal(6)=pal(6)+1;
    elseif estcomp == 7
        pal(7)=pal(7)+1;
    else estcomp == 8
        pal(8)=pal(8)+1;

    end



% MAP
[CR,I]=min(map_obj);
estcomp=I;
'MAP', estcomp


    if estcomp == 1
        map(1)=map(1)+1;
    elseif estcomp == 2
        map(2)=map(2)+1;
    elseif estcomp == 3
        map(3)=map(3)+1;
    elseif estcomp == 4
        map(4)=map(4)+1;
    elseif estcomp == 5
        map(5)=map(5)+1;
    elseif estcomp == 6
        map(6)=map(6)+1;
    elseif estcomp == 7
        map(7)=map(7)+1;
    else estcomp == 8
        map(8)=map(8)+1;

    end
    
    
% USUAL AIC
[CR,I]=min(aic_obj);
estcomp=I;
'USUAL AIC', estcomp


    if estcomp == 1
        aic(1)=aic(1)+1;
    elseif estcomp == 2
        aic(2)=aic(2)+1;
    elseif estcomp == 3
        aic(3)=aic(3)+1;
    elseif estcomp == 4
        aic(4)=aic(4)+1;
    elseif estcomp == 5
        aic(5)=aic(5)+1;
    elseif estcomp == 6
        aic(6)=aic(6)+1;
    elseif estcomp == 7
        aic(7)=aic(7)+1;
    else estcomp == 8
        aic(8)=aic(8)+1;

    end


% USUAL CORRECTED AIC
[CR,I]=min(aicc_obj);
estcomp=I;
'USUAL CORRECTED AIC', estcomp


    if estcomp == 1
        aicc(1)=aicc(1)+1;
    elseif estcomp == 2
        aicc(2)=aicc(2)+1;
    elseif estcomp == 3
        aicc(3)=aicc(3)+1;
    elseif estcomp == 4
        aicc(4)=aicc(4)+1;
    elseif estcomp == 5
        aicc(5)=aicc(5)+1;
    elseif estcomp == 6
        aicc(6)=aicc(6)+1;
    elseif estcomp == 7
        aicc(7)=aicc(7)+1;
    else estcomp == 8
        aicc(8)=aicc(8)+1;

    end

% USUAL WIC
[CR,I]=min(wic_obj);
estcomp=I;
'USUAL WIC', estcomp


    if estcomp == 1
        wic(1)=wic(1)+1;
    elseif estcomp == 2
        wic(2)=wic(2)+1;
    elseif estcomp == 3
        wic(3)=wic(3)+1;
    elseif estcomp == 4
        wic(4)=wic(4)+1;
    elseif estcomp == 5
        wic(5)=wic(5)+1;
    elseif estcomp == 6
        wic(6)=wic(6)+1;
    elseif estcomp == 7
        wic(7)=wic(7)+1;
    else estcomp == 8
        wic(8)=wic(8)+1;

    end
    

% USUAL HANNAN QUINN
[CR,I]=min(hq_obj);
estcomp=I;
'USUAL HANNAN QUINN', estcomp


    if estcomp == 1
        hq(1)=hq(1)+1;
    elseif estcomp == 2
        hq(2)=hq(2)+1;
    elseif estcomp == 3
        hq(3)=hq(3)+1;
    elseif estcomp == 4
        hq(4)=hq(4)+1;
    elseif estcomp == 5
        hq(5)=hq(5)+1;
    elseif estcomp == 6
        hq(6)=hq(6)+1;
    elseif estcomp == 7
        hq(7)=hq(7)+1;
    else estcomp == 8
        hq(8)=hq(8)+1;

    end

result_all=[bic;bicc;pal;map;aic;aicc;wic;hq]    
end % Loop for number of simulations
'-----------------------------------------------------------------------',

'-----------------------------------------------------------------------',


'USUAL BIC'
bic,

'-----------------------------------------------------------------------',

'CORRECTED BIC'
bicc,

'-----------------------------------------------------------------------',

'-----------------------------------------------------------------------',

'PAL'
pal,

'-----------------------------------------------------------------------',


'MAP'
'-----------------------------------------------------------------------',

map,

'-----------------------------------------------------------------------',


'-----------------------------------------------------------------------',

'USUAL AIC'
aic,

'-----------------------------------------------------------------------',
'USUAL CORRECTED AIC'
aicc,

'-----------------------------------------------------------------------',
'USUAL WIC'
wic,

'-----------------------------------------------------------------------',
 'USUAL HANNAN QUINN'
 hq,

'-----------------------------------------------------------------------',
result_all=[bic;bicc;pal;map;aic;aicc;wic;hq]
%save 'freqL2_AR_no_out_sig_1txt' freq_L2 -ASCII;


end % loop for numsamp