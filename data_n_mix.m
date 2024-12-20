function y=data_n_mix(numsamp,numpar,a,b,om,stdev1,stdev2);

for k=1:numsamp
    x(k)=0.0;
    for j=1:numpar
      x(k)=x(k)+a(j)*cos(om(j)*k)+b(j)*sin(om(j)*k);
    end 
    bin=binornd(1,0.6,1,1); % previously 0.6 mixing
    z=bin.*normrnd(0,stdev1,1,1)+(1-bin).*normrnd(0,stdev2,1,1);
    y(k,1)=x(k)+z;
end
    