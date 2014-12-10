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

EdgeHolga = pictures(1:51,:,:,1);
EdgeCanon = pictures(1:51,:,:,2);
EdgeScanner = pictures(1:51,:,:,3);
EdgeSony =  pictures(1:51,:,:,4);

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

sharpHolga = sum(weightHolga)
sharpCanon = sum(weightCanon)
sharpScanner = sum(weightScanner)
sharpSony = sum(weightSony)


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

% % % normalisera --> fel
% % 
% % HolgaDC = HolgaDC/ (max(max(HolgaDC)));
% % CanonDC = CanonDC/ (max(max(CanonDC)));
% % ScannerDC = ScannerDC/ (max(max(ScannerDC)));
% % SonyDC = SonyDC/ (max(max(SonyDC)));
% % 

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

average = zeros(100,1);
average2 = zeros(100,1);
average3 = zeros(100,1);
average4 = zeros(100,1);

norm_average = zeros(100,1);
norm_average2 = zeros(100,1);
norm_average3 = zeros(100,1);
norm_average4 = zeros(100,1);
 

for m=1:100
    % Create logical matrix Maskm for each value m
    Maskm = QR == m;

    % Sum the pixelvalues for each m in QR
    sum_pv = sum(sum(Maskm.*HolgaDC));
    sum_pv1 = sum(sum(Maskm.*CanonDC));
    sum_pv2 = sum(sum(Maskm.*ScannerDC));
    sum_pv3 = sum(sum(Maskm.*SonyDC));
    
    % Count the number of elements of m in QR 
    nr_objects = sum(sum(Maskm));

    average(m) = sum_pv/nr_objects;
    average2(m) = sum_pv1/nr_objects;
    average3(m) = sum_pv2/nr_objects;
    average4(m) = sum_pv3/nr_objects;

end

norm_average = average / max(average);
norm_average2 = average2 / max(average2);
norm_average3 = average3 / max(average3);
norm_average4 = average4 / max(average4);

plot(norm_average)
figure;imshow(HolgaDC)

figure;plot(norm_average2)
figure;imshow(CanonDC)

figure;plot(norm_average3)
figure;imshow(ScannerDC)

figure;plot(norm_average4)
figure;imshow(SonyDC)
%% L


% Viktfunktion
a =linspace(1,0,50);
b =linspace(0,1,50);
w = cat(2,a,b);

% Multiplicerar viktfunktionen med absolutbeloppet av skiftade
% fouriertransformen

weightHolga = w'.*abs(norm_average);
weightCanon = w'.*abs(norm_average2);
weightScanner = w'.*abs(norm_average3);
weightSony = w'.*abs(norm_average4);

sharpHolga = sum(weightHolga)
sharpCanon = sum(weightCanon)
sharpScanner = sum(weightScanner)
sharpSony = sum(weightSony)














