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

% Normalize the images
pictures(:,:,:,1) = pictures(:,:,:,1)/ max(max(pictures(:,:,:,1)));
pictures(:,:,:,2) = pictures(:,:,:,2)/ max(max(pictures(:,:,:,2)));
pictures(:,:,:,3) = pictures(:,:,:,3)/ max(max(pictures(:,:,:,3)));
pictures(:,:,:,4) = pictures(:,:,:,4)/ max(max(pictures(:,:,:,4)));


%% B

EdgeHolga = pictures(1:51,:,:,1);
EdgeCanon = pictures(1:51,:,:,2);
EdgeScanner = pictures(1:51,:,:,3);
EdgeSony =  pictures(1:51,:,:,4);


%% B-plot

figure;
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
SumEdgeHolga = sum(EdgeHolga, 1);
SumEdgeCanon = sum(EdgeCanon, 1);
SumEdgeScanner = sum(EdgeScanner, 1);
SumEdgeSony = sum(EdgeSony, 1);

% Zero padding 
Z = zeros(1,12);

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

% Viktfunktion
a =linspace(1,0,140);
b =linspace(0,1,140);
w = cat(2,a,b);

% Multiplicerar viktfunktionen med absolutbeloppet av skiftade
% fouriertransformen

weightHolga = abs(NFFT1EdgeHolga).*w;
weightCanon = abs(NFFT1EdgeCanon).*w;
weightScanner = abs(NFFT1EdgeScanner).*w;
weightSony = abs(NFFT1EdgeSony).*w;

disp('Sharpness of Holga:')
sum(weightHolga)
disp('Sharpness of Canon:')
sum(weightCanon)
disp('Sharpness of Scanner:')
sum(weightScanner)
disp('Sharpness of Sony:')
sum(weightSony)


%% H

% Adds padding aound the image
EdgeHolga2      = padarray(pictures(:,:,:,1),[128 128],0,'both');
EdgeCanon2      = padarray(pictures(:,:,:,2),[128 128],0,'both');
EdgeScanner2    = padarray(pictures(:,:,:,3),[128 128],0,'both');
EdgeSony2       = padarray(pictures(:,:,:,4),[128 128],0,'both');

% Compute fourier transform and shift
FFTEdgeHolga2 = fftshift(fft2(EdgeHolga2));
FFTEdgeCanon2 = fftshift(fft2(EdgeCanon2));
FFTEdgeScanner2 = fftshift(fft2(EdgeScanner2));
FFTEdgeSony2 = fftshift(fft2(EdgeSony2));


%% I

% Calculate abs of the images
AHolga = abs(FFTEdgeHolga2);
ACanon = abs(FFTEdgeCanon2);
AScanner = abs(FFTEdgeScanner2);
ASony = abs(FFTEdgeSony2);

HolgaDC = AHolga/AHolga(257,257); 
CanonDC = ACanon/ACanon(257,257);
ScannerDC = AScanner/AScanner(257,257);
SonyDC = ASony/ASony(257,257);


%% J

plot(HolgaDC(257,:));
figure;
plot(CanonDC(257,:));
figure;
plot(ScannerDC(257,:));
figure;
plot(SonyDC(257,:));


%% K

N = 512;
[X,Y] = meshgrid((1:N));
[T,R] = cart2pol(X-N/2,Y-N/2);

SR = R ./ R(512/2-1,1);

QR = round((SR*100));

average_Holga = zeros(100,1);
average_Canon = zeros(100,1);
average_Scanner = zeros(100,1);
average_Sony = zeros(100,1);

 
for m=1:100
    % Create logical matrix Maskm for each value m
    Maskm = QR == m;

    % Sum the pixelvalues for each m in QR
    sum_Holga = sum(sum(Maskm.*HolgaDC));
    sum_Canon = sum(sum(Maskm.*CanonDC));
    sum_Scanner = sum(sum(Maskm.*ScannerDC));
    sum_Sony = sum(sum(Maskm.*SonyDC));
    
    % Count the number of elements of m in QR 
    nr_objects = sum(sum(Maskm));

    average_Holga(m) = sum_Holga/nr_objects;
    average_Canon(m) = sum_Canon/nr_objects;
    average_Scanner(m) = sum_Scanner/nr_objects;
    average_Sony(m) = sum_Sony/nr_objects;

end

norm_average_Holga = average_Holga / max(average_Holga);
norm_average_Canon = average_Canon / max(average_Canon);
norm_average_Scanner = average_Scanner / max(average_Scanner);
norm_average_Sony = average_Sony / max(average_Sony);

%% K-plot

plot(norm_average_Holga)
title('Holga')

figure;
plot(norm_average_Canon)
title('Canon')

figure;
plot(norm_average_Scanner)
title('Scanner')

figure;
plot(norm_average_Sony)
title('Sony')


%% L

% Viktfunktion
a =linspace(1,0,50);
b =linspace(0,1,50);
w = cat(2,a,b);

% Multiplicerar viktfunktionen med absolutbeloppet av skiftade
% fouriertransformen

weightHolga = w'.*abs(norm_average_Holga);
weightCanon = w'.*abs(norm_average_Canon);
weightScanner = w'.*abs(norm_average_Scanner);
weightSony = w'.*abs(norm_average_Sony);

disp('Sharpness of Holga:')
sum(weightHolga)
disp('Sharpness of Canon:')
sum(weightCanon)
disp('Sharpness of Scanner:')
sum(weightScanner)
disp('Sharpness of Sony:')
sum(weightSony)


%% Part 2: Autofocus with Fourier Transforms

%% A and B

load('winsuint8.mat');

sharpness = zeros(192);

N = 64;
[X,Y] = meshgrid((1:N));
[T,R] = cart2pol(X-N/2,Y-N/2);

for i = 1:192
    
   % Adds padding aound the image
   patch = padarray(winsuint8(:,:,i),[16 16],0,'both');
   
   FFTpatch = abs(fftshift(fft2(patch)));
   
   FFTpatchDC =  FFTpatch / (FFTpatch(33,33));
   
   sharpness(i) = sum(sum(FFTpatchDC.*R));
    
end

plot(sharpness)

disp('Sharpest images is:')
find(sharpness == max(sharpness(:)))


