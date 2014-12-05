%% TNM087 - Labb 3
% Operations in the Fourier domain
% Christoffer Engelbrektsson & Johan Nordin
% MT3


%% Part 1: Comparison of different optical systems

%% A

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


% Normalize the images
pictures(:,:,:,1) = pictures(:,:,:,1)/ max(max(pictures(:,:,:,1)));
pictures(:,:,:,2) = pictures(:,:,:,2)/ max(max(pictures(:,:,:,2)));
pictures(:,:,:,3) = pictures(:,:,:,3)/ max(max(pictures(:,:,:,3)));
pictures(:,:,:,4) = pictures(:,:,:,4)/ max(max(pictures(:,:,:,4)));

% figure;
% %montage(pictures)

%% B

EdgeHolga = pictures(:,103:153,:,1);
EdgeCanon = pictures(:,103:153,:,2);
EdgeScanner = pictures(:,103:153,:,3);
EdgeSony =  pictures(:,103:153,:,4);

%% B-plot

figure;
%montage(EdgeHolga, EdgeCanon, EdgeScanner, EdgeSony);
subplot(1,4,1);
imshow(EdgeHolga);
subplot(1,4,2);
imshow(EdgeCanon);
subplot(1,4,3);
imshow(EdgeScanner);
subplot(1,4,4);
imshow(EdgeSony);

%% C

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

%% D

% FFT
FFT1EdgeHolga = fftshift(fft(SumEdgeHolgaZeroPad));
FFT1EdgeCanon = fftshift(fft(SumEdgeCanonZeroPad));
FFT1EdgeScanner = fftshift(fft(SumEdgeScannerZeroPad));
FFT1EdgeSony = fftshift(fft(SumEdgeSonyZeroPad));




%% D-plot

figure;
subplot(1,3,1);
plot(abs(FFT1EdgeHolga));
legend('Abs')
subplot(1,3,2);
plot(real(FFT1EdgeHolga));
legend('Real')
title('FFT1EdgeHolga')
subplot(1,3,3);
plot(imag(FFT1EdgeHolga));
legend('Imag')

figure;
subplot(1,3,1);
plot(abs(FFT1EdgeCanon));
legend('Abs')
subplot(1,3,2);
plot(real(FFT1EdgeCanon));
legend('Real')
title('FFT1EdgeCanon')
subplot(1,3,3);
plot(imag(FFT1EdgeCanon));
legend('Imag')

figure;
subplot(1,3,1);
plot(abs(FFT1EdgeScanner));
legend('Abs')
subplot(1,3,2);
plot(real(FFT1EdgeScanner));
legend('Real')
title('FFT1EdgeScanner')
subplot(1,3,3);
plot(imag(FFT1EdgeScanner));
legend('Imag')

figure;
subplot(1,3,1);
plot(abs(FFT1EdgeSony));
legend('Abs')
subplot(1,3,2);
plot(real(FFT1EdgeSony));
legend('Real')
title('FFT1EdgeSony')
subplot(1,3,3);
plot(imag(FFT1EdgeSony));
legend('Imag')


%% E

% DC komponenten är impulsen som sker vid frekvensen 0


%% F

% Varför är detta användbart?

NFFT1EdgeHolga = FFT1EdgeHolga/max(FFT1EdgeHolga);
NFFT1EdgeCanon = FFT1EdgeCanon/max(FFT1EdgeCanon);
NFFT1EdgeScanner = FFT1EdgeScanner/max(FFT1EdgeScanner);
NFFT1EdgeSony = FFT1EdgeSony/max(FFT1EdgeSony);

%% F-plot

figure;
subplot(1,3,1);
plot(abs(NFFT1EdgeHolga));
legend('Abs')
subplot(1,3,2);
plot(real(NFFT1EdgeHolga));
legend('Real')
title('FFT1EdgeHolga - Normerad')
subplot(1,3,3);
plot(imag(NFFT1EdgeHolga));
legend('Imag')

figure;
subplot(1,3,1);
plot(abs(NFFT1EdgeCanon));
legend('Abs')
subplot(1,3,2);
plot(real(NFFT1EdgeCanon));
legend('Real')
title('FFT1EdgeCanon - Normerad')
subplot(1,3,3);
plot(imag(NFFT1EdgeCanon));
legend('Imag')

figure;
subplot(1,3,1);
plot(abs(NFFT1EdgeScanner));
legend('Abs')
subplot(1,3,2);
plot(real(NFFT1EdgeScanner));
legend('Real')
title('FFT1EdgeScanner - Normerad')
subplot(1,3,3);
plot(imag(NFFT1EdgeScanner));
legend('Imag')

figure;
subplot(1,3,1);
plot(abs(NFFT1EdgeSony));
legend('Abs')
subplot(1,3,2);
plot(real(NFFT1EdgeSony));
legend('Real')
title('FFT1EdgeSony - Normerad')
subplot(1,3,3);
plot(imag(NFFT1EdgeSony));
legend('Imag')



%% G

%% H

%% I

%% J

%% K

%% L














