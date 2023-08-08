
function [psn,Se] = Psnr_estimation(x,p)

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
psn=10*log10((255^2)/Se);
% Mxx=[Se psn];
% figure,plot(Mxx,'r*-');
% axis([0 5 0 50]);
% grid on;
% title('PSNR Plot');
