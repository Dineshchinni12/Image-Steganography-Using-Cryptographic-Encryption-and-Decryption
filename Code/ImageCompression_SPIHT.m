function [out] = ImageCompression_SPIHT(infile,bx,by,outfile)


a=infile;
    double_a=double(a);
    ax=size(a,1)-mod(size(a,1),bx);
    ay=size(a,2)-mod(size(a,2),by);
    out_rgb=zeros(ax,ay,3);
    % -----------------------  RED component  ------------------   
    dvalue=double_a(:,:,1);
    
    
    dx=size(dvalue,1);
    dy=size(dvalue,2);
    % if input image size is not a multiple of block size image is resized
    modx=mod(dx,bx);
    mody=mod(dy,by);
    dvalue=dvalue(1:dx-modx,1:dy-mody);
    % the new input image dimensions (pixels)
    dx=dx-modx;
    dy=dy-mody;
    % number of non overlapping blocks required to cover 
    % the entire input image
    nbx=size(dvalue,1)/bx;
    nby=size(dvalue,2)/by;
    
    % the output compressed image
    matrice=zeros(bx,by);
    % the compressed data
    m_u=zeros(nbx,nby);
    m_l=zeros(nbx,nby);
    mat_log=logical(zeros(bx,by));
    
    posbx=1;
    for ii=1:bx:dx
        posby=1;
        for jj=1:by:dy
            % the current block
            blocco=dvalue(ii:ii+bx-1,jj:jj+by-1);
            % the average gray level of the current block
            m=mean(mean(blocco));
            % the logical matrix correspoending to the current block
            blocco_binario=(blocco>=m);
            % the number of pixel (of the current block) whose gray level
            % is greater than the average gray level of the current block
            K=sum(sum(double(blocco_binario)));
            % the average gray level of pixels whose level is GREATER than
            % the block average gray level
            mu=sum(sum(double(blocco_binario).*blocco))/K;
            % the average gray level of pixels whose level is SMALLER than
            % the block average gray level
            if K==bx*by
                ml=0;
            else
                ml=sum(sum(double(~blocco_binario).*blocco))/(bx*by-K);
            end
            % the COMPRESSED DATA which correspond to the input image
            m_u(posbx,posby)=mu;                            %---> the m_u matrix (see the cited reference)
            m_l(posbx,posby)=ml;                            %---> the m_l matrix 
            mat_log(ii:ii+bx-1,jj:jj+by-1)=blocco_binario;  %---> the logical matrix
            % the compressed image
            matrice(ii:ii+bx-1,jj:jj+by-1)=(double(blocco_binario).*mu)+(double(~blocco_binario).*ml);            
            posby=posby+1;
        end
        posbx=posbx+1;
    end
    out_rgb(:,:,1)=matrice;
    % -----------------------  GREEN component  ------------------   
    dvalue=double_a(:,:,2);
    
    
    dx=size(dvalue,1);
    dy=size(dvalue,2);
    % if input image size is not a multiple of block size image is resized
    modx=mod(dx,bx);
    mody=mod(dy,by);
    dvalue=dvalue(1:dx-modx,1:dy-mody);
    % the new input image dimensions (pixels)
    dx=dx-modx;
    dy=dy-mody;
    % number of non overlapping blocks required to cover 
    % the entire input image
    nbx=size(dvalue,1)/bx;
    nby=size(dvalue,2)/by;
    
    % the output compressed image
    matrice=zeros(bx,by);
    % the compressed data
    m_u=zeros(nbx,nby);
    m_l=zeros(nbx,nby);
    mat_log=logical(zeros(bx,by));
    
    posbx=1;
    for ii=1:bx:dx
        posby=1;
        for jj=1:by:dy
            % the current block
            blocco=dvalue(ii:ii+bx-1,jj:jj+by-1);
            % the average gray level of the current block
            m=mean(mean(blocco));
            % the logical matrix correspoending to the current block
            blocco_binario=(blocco>=m);
            % the number of pixel (of the current block) whose gray level
            % is greater than the average gray level of the current block
            K=sum(sum(double(blocco_binario)));
            % the average gray level of pixels whose level is GREATER than
            % the block average gray level
            mu=sum(sum(double(blocco_binario).*blocco))/K;
            % the average gray level of pixels whose level is SMALLER than
            % the block average gray level
            if K==bx*by
                ml=0;
            else
                ml=sum(sum(double(~blocco_binario).*blocco))/(bx*by-K);
            end
            % the COMPRESSED DATA which correspond to the input image
            m_u(posbx,posby)=mu;                            %---> the m_u matrix (see the cited reference)
            m_l(posbx,posby)=ml;                            %---> the m_l matrix 
            mat_log(ii:ii+bx-1,jj:jj+by-1)=blocco_binario;  %---> the logical matrix
            % the compressed image
            matrice(ii:ii+bx-1,jj:jj+by-1)=(double(blocco_binario).*mu)+(double(~blocco_binario).*ml);            
            posby=posby+1;
        end
        posbx=posbx+1;
    end
    out_rgb(:,:,2)=matrice;
    % -----------------------  BLUE component  ------------------   
    dvalue=double_a(:,:,3);
    
   
    dx=size(dvalue,1);
    dy=size(dvalue,2);
    % if input image size is not a multiple of block size image is resized
    modx=mod(dx,bx);
    mody=mod(dy,by);
    dvalue=dvalue(1:dx-modx,1:dy-mody);
    % the new input image dimensions (pixels)
    dx=dx-modx;
    dy=dy-mody;
    % number of non overlapping blocks required to cover 
    % the entire input image
    nbx=size(dvalue,1)/bx;
    nby=size(dvalue,2)/by;
    
    % the output compressed image
    matrice=zeros(bx,by);
    % the compressed data
    m_u=zeros(nbx,nby);
    m_l=zeros(nbx,nby);
    mat_log=logical(zeros(bx,by));
    
    posbx=1;
    for ii=1:bx:dx
        posby=1;
        for jj=1:by:dy
            % the current block
            blocco=dvalue(ii:ii+bx-1,jj:jj+by-1);
            % the average gray level of the current block
            m=mean(mean(blocco));
            % the logical matrix correspoending to the current block
            blocco_binario=(blocco>=m);
            % the number of pixel (of the current block) whose gray level
            % is greater than the average gray level of the current block
            K=sum(sum(double(blocco_binario)));
            % the average gray level of pixels whose level is GREATER than
            % the block average gray level
            mu=sum(sum(double(blocco_binario).*blocco))/K;
            % the average gray level of pixels whose level is SMALLER than
            % the block average gray level
            if K==bx*by
                ml=0;
            else
                ml=sum(sum(double(~blocco_binario).*blocco))/(bx*by-K);
            end
            % the COMPRESSED DATA which correspond to the input image
            m_u(posbx,posby)=mu;                            %---> the m_u matrix (see the cited reference)
            m_l(posbx,posby)=ml;                            %---> the m_l matrix 
            mat_log(ii:ii+bx-1,jj:jj+by-1)=blocco_binario;  %---> the logical matrix
            % the compressed image
            matrice(ii:ii+bx-1,jj:jj+by-1)=(double(blocco_binario).*mu)+(double(~blocco_binario).*ml);            
            posby=posby+1;
        end
        posbx=posbx+1;
    end
    out_rgb(:,:,3)=matrice;
    %------------------------------------------ -----------------
    %-----------------------------------------------------------
    %-----------------------------------------------------------
    if isa(a,'uint8')
        out=uint8(out_rgb);
        figure('Name','Compressed image');
        imshow(out);
        imwrite(out, outfile);
        return
    end
    
    if isa(a,'uint16')
        out=uint16(out_rgb);
        figure('Name','Compressed image');
        imshow(out);
        imwrite(out, outfile);
        return
    end
    
    if isa(a,'double')
        out=(out_rgb);
        figure('Name','Compressed image');
        imshow(out);
        imwrite(out, outfile);
        return
    end
    
end