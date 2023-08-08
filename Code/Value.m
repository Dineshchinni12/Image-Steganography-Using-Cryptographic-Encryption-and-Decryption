function [psn,See] = Value(x,p)

si=size(x);
m=si(1);
n=si(2);
x=double(x);
p=double(p);
Se=0;
for i=1:m
    for j=1:n
    Se=Se+(x(i,j)-p(i,j))^2;
    end
end
Se=Se/(m*n);
psn=20*log10((255^2)/Se);
See=Se/1000;
