function xn = normalizedstd(x,missing)
if nargin<2
    temp = isnan(x);
else
    temp = (x==missing);
end
%空值置零
y = x;
y(temp) = 0;    
%减均值除以标准差
mn = sum(y)./sum(~temp);
sd = sqrt(sum((y-repmat(mn,size(y,1),1)).^2.*(~temp))./(sum(~temp)-1));
sd(sd==0) = eps;
xn = (x-repmat(mn,size(x,1),1))./repmat(sd,size(x,1),1);
% X=x;
% xn = (X-repmat(min(X),size(X,1),1))./(repmat(max(X),size(X,1),1)-repmat(min(X),size(X,1),1)+eps);

xn = (x-repmat(mn,size(x,1),1))./repmat(sd,size(x,1),1);

%missing说明用什么填充空值
if nargin==2
xn(temp) = missing;  
end

