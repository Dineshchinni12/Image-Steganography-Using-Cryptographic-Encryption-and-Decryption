function [Encode, R, C, Cx,Ouu] = FEC_Encode(Image)
M=8;
[R, C, Cx] =size(Image);
V_im = zeros(1, R*C*Cx);
for i=1:Cx,
    for j=1:C,
        V_im((i-1)*R*C+(j-1)*R +1 : (i-1)*R*C+j*R) = Image(:,j,i);
    end
end

Tm = dec2base(V_im,M);
[row, col] = size(Tm);
for i =1:col,
    Encode(row*(i-1)+1:row*i)=str2num(Tm(:,i));  
end
Ouu=Image;
