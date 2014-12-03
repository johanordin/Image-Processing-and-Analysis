%% TNM087 - Labb 1 
% Christoffer Engelbrektsson & Johan Nordin
% MT3


%% Lab 3: Operations in the Fourier domain


%% Part 1: Comparison of different optical systems


pictures(:,:,:,1) = rgb2gray(im2double(imread('HalfHolga.jpg')));
pictures(:,:,:,2) = rgb2gray(im2double(imread('HalfCanon.jpg')));
pictures(:,:,:,3) = rgb2gray(im2double(imread('HalfSony.jpg')));
pictures(:,:,:,4) = rgb2gray(im2double(imread('HalfScanner.jpg')));

% % Plot the images
% subplot(1,4,1);
% imshow(pictures(:,:,:,1));
% subplot(1,4,2);
% imshow(pictures(:,:,:,2));
% subplot(1,4,3);
% imshow(pictures(:,:,:,3));
% subplot(1,4,4);
% imshow(pictures(:,:,:,4));

%
% montage(pictures)

% % Normalize the images
pictures(:,:,:,1) = pictures(:,:,:,1)/ max(max(pictures(:,:,:,1)));
pictures(:,:,:,2) = pictures(:,:,:,2)/ max(max(pictures(:,:,:,2)));
pictures(:,:,:,3) = pictures(:,:,:,3)/ max(max(pictures(:,:,:,3)));
pictures(:,:,:,4) = pictures(:,:,:,4)/ max(max(pictures(:,:,:,4)));

% figure;
% %montage(pictures)

EdgeHolga = pictures(:,103:153,:,1);
EdgeCanon = pictures(:,103:153,:,2);
EdgeScanner = pictures(:,103:153,:,3);
EdgeSony =  pictures(:,103:153,:,4);

% figure;
% %montage(EdgeHolga, EdgeCanon, EdgeScanner, EdgeSony);
% subplot(1,4,1);
% imshow(EdgeHolga);
% subplot(1,4,2);
% imshow(EdgeCanon);
% subplot(1,4,3);
% imshow(EdgeScanner);
% subplot(1,4,4);
% imshow(EdgeSony);

% sum 
SumEdgeHolga = sum(EdgeHolga, 2)';
SumEdgeCanon = sum(EdgeCanon, 2)';
SumEdgeScanner = sum(EdgeScanner, 2)';
SumEdgeSony = sum(EdgeSony, 2)';

% Zero padding 
Z = zeros(12,1)';

SumEdgeHolgaZeroPad = cat(2, Z, SumEdgeHolga, Z);
SumEdgeCanonZeroPad = cat(2, Z, SumEdgeCanon, Z);
SumEdgeScannerZeroPad = cat(2, Z, SumEdgeScanner, Z);
SumEdgeSonyZeroPad = cat(2, Z, SumEdgeSony, Z);



FFT1EdgeHolga = fftshift(fft(SumEdgeHolgaZeroPad));
FFT1EdgeCanon = fftshift(fft(SumEdgeCanonZeroPad));

FFT1EdgeScanner = fftshift(fft(SumEdgeScannerZeroPad));
FFT1EdgeSony = fftshift(fft(SumEdgeSonyZeroPad));


plot(abs(FFT1EdgeHolga));
hold on ;
plot(real(FFT1EdgeHolga));
plot(imag(FFT1EdgeHolga));












