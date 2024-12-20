function inits=per_maxm(yadj,numsamp,gridnum);
start_om=0;
end_om=pi;
for k=1:gridnum+1
    om_grid(k)=start_om+((end_om-start_om)/gridnum)*(k-1);
    per(k)=0.0;
    for kk=1:numsamp
        per(k)=per(k)+yadj(kk)*exp(-i*kk*om_grid(k));
    end
    per(k)=per(k)/numsamp;
    per(k)=abs(per(k))^2;
end
[Y,I]=max(per);
om=om_grid(I);

for kk=1:numsamp
    a_om(kk,1)=cos(om*kk);
    a_om(kk,2)=sin(om*kk);
    yvec(kk,1)=yadj(kk);
end
lin=(pinv(a_om'*a_om))*a_om'*yvec;
inits=[lin(1);lin(2);om];