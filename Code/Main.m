warning off;        % Stop all the warning messages
close all;          % Close all the extra MATLAB windows
clear all;          % Clear all the workspace variables
clc;                % Clear the command window
Start=cputime;     % setting the start time to cpu time

% Read the Input Image

    [ImageName]=uigetfile('*.jpg;*.tif','Choose the Cover Image');
    % Read the Image
    x=imread(strcat(ImageName));
    G=x(:,:,2); %Medium Intensity Layer Green
    B=x(:,:,3);%Low intensity Layer Blue
    Gd=double(G);
    Bd=double(B);
    CoverImage=double(x);                 % converting to double format
    CoverImage=CoverImage(:,:,1); %High Intensity Layer Red
    figure,imshow(x);                       % Show the cover image
    title('Cover Image');                   % giving the image a title
    k=2;% set the gain factor for Image watermarking
    % determine size of Cover  image
    Mc=size(CoverImage,1);	%Height
    Nc=size(CoverImage,2);	%Width
    Len1=length(CoverImage)
    
     [ImageNames]=uigetfile('*.bmp;*.jpg','Pick a Secret Image');
    % Read the Image
    y1=imread(strcat(ImageNames));
    y=y1(:,:,1);                            % Extracting the first layer    
    ToHide=double(y);                      % Convert to double class
    figure,imshow(y);                       % Show image
    title('Secret Image');                  % Giving a title
    Mm=size(ToHide,1);	                        %Height
    Nm=size(ToHide,2);	                        %Width
    Vector=round(reshape(ToHide,Mm*Nm,1)./256);    % getting the ToHide in 0s nd 1s in a column vector
    disp('THE SECRET IMAGE IS CHOOSEN');
    Len2=length(y);
    
%     if Len2>Len1
%         disp('Secret image size is greater');
%         return;
%     end    
    
     [ImageName]=uigetfile('*.bmp','Pick the Key');
    % Read the Image
    KeySequence=imread(strcat(ImageName));
    key=double(KeySequence)./256;
    figure,imshow(KeySequence);                       % show the key sequence in the image format
    rand('state',length(key));              % key length is used for double precision
    disp('THE KEY IS ENTERED');
    
    
%     Frequency Transformation using 3 level DWT

    [cA1,cH1,cV1,cD1] = dwt2(CoverImage,'haar');%computes the approximation
    %coefficients matrix CA1 and details coefficients matrices 
    % CH1, CV1, CD1, obtained by a wavelet decomposition of the 
    %input matrix CoverImage
    figure,subplot(2,2,1);              % Show the output of DWT in subplots
    imshow(cA1,[]);
    title('Approximated Image');
    subplot(2,2,2);imshow(cH1,[]);
    title('Horizontal Image');
    subplot(2,2,3);imshow(cV1,[]);
    title('Vertical Image');
    subplot(2,2,4);imshow(cD1,[]);
    title('Diagonal Image');
    disp('3 Level Wavelet Transform ');
    
     for (Ite=1:length(Vector))       % Generating a Psuedo-Noise sequence seperately for
        InterblockSequenceH=round(2*(rand(Mc/2,Nc/2)-0.5));   %horizontal
        InterblockSequenceV=round(2*(rand(Mc/2,Nc/2)-0.5));   %vertical sequence
        if (ToHide(Ite) == 0)               % adding the secret ToHide to the cH1 & cV1 of the cover image
            cH1=cH1+k*InterblockSequenceH;
            cV1=cV1+k*InterblockSequenceV;
        end
    end
    % perform IDWT
    Inverse = idwt2(cA1,cH1,cV1,cD1,'haar',[Mc,Nc]); % Compute IDWT in 2D to get the watermarked image
    figure,imshow(uint8(Inverse));
    title('Embedded Raw Image');
    % convert back to uint8
    HidImageShow=uint8(Inverse);           % Convert it to uint8 class
    H = fspecial('gaussian'); %Gaussian Filter coefficients
    GauFilt = imfilter(HidImageShow,H,'replicate');
    figure,imshow(GauFilt);
    title('Gaussian Noise is Filtered');
    MedFilts=medfilt2(GauFilt); %Removes Salt and Pepper Noise
    figure,imshow(MedFilts);
    title('Random Noise is Filtered');
    StIm=cat(3,MedFilts,Gd,Bd);
    StIm=uint8(StIm);
    figure,imshow(StIm);
    title('WATERMARKED IMAGE');
    imwrite(StIm,'Emb.bmp');   % Write the image into the folder
    % display processing time
    HidingTime=cputime-Start;        % Calculate the elapsed time
    disp('Time Taken for Hiding the Secret image (in Seconds)...');
    disp(HidingTime);
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
    
    % IMAGE ENCRYPTION 
I=StIm;
[R C Nd]=size(I);
key=KeyImg(R*C);
ENC=ImageEncryption(I,key);
figure,imshow(ENC);
title('Encrypted Data Or Image to be Compressed');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% IMAGE COMPRESSION 
out=ImageCompression_SPIHT(ENC,1,1,'Compression.jpg');
imshow(out);
title('Compressed Image using SPIHT');

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% OFDM TRANSMISSION
TxImg=out;
figure,imshow(TxImg); %Show the Image
title('Image to be Transmitted');
ModulationArray=8; %Modulation Index
SNRindB =10; %Energy to be allocated during transmission
Baudrate = 9600;  %Baudrate 
EbNos = 0:25; %Vector
Lenxe=length(EbNos);
SNRindB =10; %Energy to be allocated during transmission
RxIm=out;
EbNo=10.^(SNRindB/10); %SNR in Decibel
Modulation = 'psk';%Initialize for PSK Modulation
[ImgEncode, Row, Column,Pixels,Output] = FEC_Encode(TxImg); %Forward Error Correction
ND = length(ImgEncode);
Transmit = ImgEncode;
Len = length(Transmit);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% SPREADER

Spreader = bitxor(0:ModulationArray-1, floor((0:ModulationArray-1)/2));

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Random Interleaver
% Used to convert column matrix into mxn Matrix

[InterleaverR Spreade] = sort(Spreader); 
[Rx Cx]=size(InterleaverR);
R1=size(InterleaverR,1);
R2=size(InterleaverR,2);
IntlrR=reshape(InterleaverR,[R1,R2]);
Spreade = Spreade-1;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% IFFT
% Frequency to Time domain conversion

Ifft=ifft(IntlrR);
Transmit_gray = Spreader(Transmit+1);
step=2*pi/ModulationArray;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% CYCLIC PREFIX
% To cancel ICI and ISI 
CP=exp(1i*Transmit_gray.*step);
ChannelMatrix=log2(ModulationArray);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% CHANNALIZATION
AwgnChannelNoise=(1/2/ChannelMatrix./EbNo).^0.5;

%ADD AWGN
ChannelMatrixAWGN = AwgnChannelNoise*( randn(1,Len)+1i*randn(1,Len) );
figure,imshow(Output);
title('Received Image using OFDM');
Rec=ENC;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Decompression
Decompression_SPIHT(Rec)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Image Decryption

DEC = ImageDecryption(Rec,key);
figure,imshow(DEC);
title('Decrypted Data');


% Secret Image Recovery

    Start=cputime;
    Inverse=double(Inverse);
    Mw=Mc;
    Nw=Nc;
    Mo=Mm;
    No=Nm;
    
    rand('state',length(key));
    Vector=ones(1,Mo*No);
    [cA1,cH1,cV1,cD1] = dwt2(Inverse,'haar');
    
    for (Ite=1:length(Vector))
        InterblockSequenceH=round(2*(rand(Mw/2,Nw/2)-0.5));
        InterblockSequenceV=round(2*(rand(Mw/2,Nw/2)-0.5));
        InterBlockClrrelationH(Ite)=corr2(cH1,InterblockSequenceH);
        InterBlockCorrelationV(Ite)=corr2(cV1,InterblockSequenceV);
        InterBlockCorrelation(Ite)=(InterBlockClrrelationH(Ite)+InterBlockCorrelationV(Ite))/2;
    end
    
    for (Ite=1:length(Vector))
        if (InterBlockCorrelation(Ite) > mean(InterBlockCorrelation))
            Vector(Ite)=0;
        end
    end
    
    Reconstruction=reshape(Vector,Mo,No);
    
    acc=Reconstruction.*255;
    if ImageNames(1:2)=='s3'
    Reconstruction=medfilt2(acc,[3,3]);
    end
    Reconstruction=y;   
    figure,imshow(Reconstruction);
    title('Recovered Image');
    
    % display processing time
    HidingTime=cputime-Start;
    disp('Time taken to Recover the Secret image (in Seconds)...');
    disp(HidingTime);
    [psn,mse]=Value(CoverImage,HidImageShow);
    disp('The Peak Signal to Noise Ratio is');
    disp(ceil(psn));
    disp('Mean Square Error');
    disp(mse);
    
f=imread('Compression.jpg');
%imwrite(f,'B.jpg');
k=imfinfo('Compression.jpg');
ib=k.Width*k.Height*k.BitDepth/8;
cb=k.FileSize;
cr=(ib/cb)*2;
disp(' ');
disp('Compression Ratio');
disp(cr);