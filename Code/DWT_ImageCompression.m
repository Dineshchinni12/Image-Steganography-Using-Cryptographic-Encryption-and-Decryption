function [im] = DWT_ImageCompression (infile)
coeff=2000;
    a = imread(infile);    
    if isa(a(:,:,1),'uint8')
        red = double(a(:,:,1));
        green = double(a(:,:,2));
        blue = double(a(:,:,3));
        
        red_dct=dct2(red);
        green_dct=dct2(green);
        blue_dct=dct2(blue);
        
        red_pow   = red_dct.^2;
        green_pow = green_dct.^2;
        blue_pow  = blue_dct.^2;
        
        red_pow=red_pow(:);
        green_pow=green_pow(:);
        blue_pow=blue_pow(:);
        
        [B_r,index_r]=sort(red_pow);
        [B_g,index_g]=sort(green_pow);
        [B_b,index_b]=sort(blue_pow);
        
        index_r=flipud(index_r);
        index_g=flipud(index_g);
        index_b=flipud(index_b);
        
        im_dct_r=zeros(size(red));
        im_dct_g=zeros(size(green));
        im_dct_b=zeros(size(blue));
        
        for ii=1:coeff
            im_dct_r(index_r(ii))=red_dct(index_r(ii));
            im_dct_g(index_g(ii))=green_dct(index_g(ii));
            im_dct_b(index_b(ii))=blue_dct(index_b(ii));
        end
        
        im_r=idct2(im_dct_r);
        im_g=idct2(im_dct_g);
        im_b=idct2(im_dct_b);
        
        im=zeros(size(red,1),size(red,2),3);
        im(:,:,1)=im_r;
        im(:,:,2)=im_g;
        im(:,:,3)=im_b;
        
        im=uint8(im);
        
%         
%         
%         figure('Name','Output image');
%         imshow(im);
        
        return;
    end
    
    if isa(a(:,:,1),'uint16')
        red = double(a(:,:,1));
        green = double(a(:,:,2));
        blue = double(a(:,:,3));
        
        red_dct=dct2(red);
        green_dct=dct2(green);
        blue_dct=dct2(blue);
        
        red_pow   = red_dct.^2;
        green_pow = green_dct.^2;
        blue_pow  = blue_dct.^2;
        
        red_pow=red_pow(:);
        green_pow=green_pow(:);
        blue_pow=blue_pow(:);
        
        [B_r,index_r]=sort(red_pow);
        [B_g,index_g]=sort(green_pow);
        [B_b,index_b]=sort(blue_pow);
        
        index_r=flipud(index_r);
        index_g=flipud(index_g);
        index_b=flipud(index_b);
        
        im_dct_r=zeros(size(red));
        im_dct_g=zeros(size(green));
        im_dct_b=zeros(size(blue));
        
        for ii=1:coeff
            im_dct_r(index_r(ii))=red_dct(index_r(ii));
            im_dct_g(index_g(ii))=green_dct(index_g(ii));
            im_dct_b(index_b(ii))=blue_dct(index_b(ii));
        end
        
        im_r=idct2(im_dct_r);
        im_g=idct2(im_dct_g);
        im_b=idct2(im_dct_b);
        
        im=zeros(size(red,1),size(red,2),3);
        im(:,:,1)=im_r;
        im(:,:,2)=im_g;
        im(:,:,3)=im_b;
        
        im=uint16(im);
        
%         imwrite(im, outfile);       
        
        
%         figure('Name','Output image');
%         imshow(im);
        
        
        return;
    end
    
    if isa(a(:,:,1),'double')
        red = double(a(:,:,1));
        green = double(a(:,:,2));
        blue = double(a(:,:,3));
        
        red_dct=dct2(red);
        green_dct=dct2(green);
        blue_dct=dct2(blue);
        
        red_pow   = red_dct.^2;
        green_pow = green_dct.^2;
        blue_pow  = blue_dct.^2;
        
        red_pow=red_pow(:);
        green_pow=green_pow(:);
        blue_pow=blue_pow(:);
        
        [B_r,index_r]=sort(red_pow);
        [B_g,index_g]=sort(green_pow);
        [B_b,index_b]=sort(blue_pow);
        
        index_r=flipud(index_r);
        index_g=flipud(index_g);
        index_b=flipud(index_b);
        
        im_dct_r=zeros(size(red));
        im_dct_g=zeros(size(green));
        im_dct_b=zeros(size(blue));
        
        for ii=1:coeff
            im_dct_r(index_r(ii))=red_dct(index_r(ii));
            im_dct_g(index_g(ii))=green_dct(index_g(ii));
            im_dct_b(index_b(ii))=blue_dct(index_b(ii));
        end
        
        im_r=idct2(im_dct_r);
        im_g=idct2(im_dct_g);
        im_b=idct2(im_dct_b);
        
        im=zeros(size(red,1),size(red,2),3);
        im(:,:,1)=im_r;
        im(:,:,2)=im_g;
        im(:,:,3)=im_b;  
        
        imshow(im);
        return;
    end
end





