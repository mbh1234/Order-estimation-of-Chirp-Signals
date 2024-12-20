function y=data_nor(numsamp,numpar,a,b,om,stdev);
for k=1:numsamp
    x(k)=0.0;
    for j=1:numpar
      x(k)=x(k)+a(j)*cos(om(j)*k)+b(j)*sin(om(j)*k);
    end 
    y(k,1)=x(k)+normrnd(0,stdev);
end

    