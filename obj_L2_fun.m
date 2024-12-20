function [obj_L2_fun] = fun(initial,y,numsamp,comp)

for pa=1:comp
    a(pa)=initial(pa,1);  
    b(pa)=initial(pa,2);
    om(pa)=initial(pa,3);
end
obj_L2_fun=0.0;
for ln=1:numsamp
    yp(ln)=0.0;
    for lp=1:comp
        yp(ln)=yp(ln)+(a(lp)*cos(om(lp)*ln)+b(lp)*sin(om(lp)*ln));
    end
    diff_sq(ln)=(y(ln)-yp(ln))^2;
        obj_L2_fun=obj_L2_fun+diff_sq(ln);
end

% CALCULATION OF VANDERMONDE MATRIX AND CONCENTRATED LIKELIHOOD
% for lr=1:numsamp
%  for lc=1:numpar
%   a(lr,lc)=cos(freqd(lc*exp(i*trfreq(lc)*lr);
%  end;
% end;
% CALCULATION OF OBJECTIVE FUNCTION
%obj=-(yd'*a*pinv(a'*a)*a'*yd);
