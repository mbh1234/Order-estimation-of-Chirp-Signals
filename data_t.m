function y=data_t(numsamp,a,b,om,df,numpar);
for k=1:numsamp
    x(k)=0.0;
    for j=1:numpar
      x(k)=x(k)+a(j)*cos(om(j)*k)+b(j)*sin(om(j)*k);
    end 
      e=trnd(df);
      y(k,1)=x(k)+e;
end
   

