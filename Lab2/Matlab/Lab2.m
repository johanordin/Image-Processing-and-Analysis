%% TNM087 - Labb 2: Operations in the image domain
% Christoffer Engelbrektsson & Johan Nordin
% MT3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Part 1: Vignette and uneven illumination
%% A

% Read in the pictures, convert to double and store in 4-dim matrix
pictures(:,:,:,1) = im2double(imread('CWhite1.jpg'));
pictures(:,:,:,2) = im2double(imread('HWhite1.jpg'));

% Resize the images to half the size
pictures = imresize(pictures, 0.5);

N = 512;
[X,Y] = meshgrid((1:N));
[T,R] = cart2pol(X-N/2,Y-N/2);

% Plot the vector R(N/2,:) and use imshow(R,[])
subplot(1,2,1)
plot(R(N/2,:));
subplot(1,2,2);
imshow(R, []);

%% B

SR = R ./ R(512/2-1,1);

QR = round((SR*100));

average = zeros(100,2);
norm_average = zeros(100,2);
 
for i=1:2
    for m=1:100
        % Create logical matrix Maskm for each value m
        Maskm = QR == m;

        % Sum the pixelvalues for each m in QR
        sum_pv = sum(sum(Maskm.*pictures(:,:,:,i)));

        % Count the number of elements of m in QR 
        nr_objects = sum(sum(Maskm));

        average(m,i) = sum_pv/nr_objects;
    end
    norm_average(:,i) = average(:,i) / max(average(:,i));
end

% Plot the images 
subplot(2,1,1);
plot(norm_average(:,1));
title('Canon lens');
ylabel('intensity');
subplot(2,1,2);
plot(norm_average(:,2));
title('Holga lens');
xlabel('radius');ylabel('intensity');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2: Correlation 

redEyes = im2double(imread('BoldRedEye.JPG'));

% Extract the read channel
redChannel = redEyes(:,:,1);

% Load in the mask
load('RedEyeMask');

sq_filter_32 = ones(48);

MFilterImage = imfilter(redChannel, sq_filter_32);
EyeFilterImage = imfilter(redChannel, RedEyeMask);

ratio = EyeFilterImage./MFilterImage;

% Normalize the ratio matrix
imshow(ratio / max(max(ratio)))
figure;

quant = quantile(quantile(ratio, 0.98), 0.98);

ratio = ratio >= quant;

imshow(ratio);
figure;

% return two-dimensional eight-connected neighborhood 
BW = ratio.*imregionalmax(ratio);

imshow(BW)
imshowpair(BW, redEyes)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 3: Registration
Gim = im2double(imread('GCPins512.jpg'));
Him = im2double(imread('GHPins512.jpg'));

% Minsta antalet punkter som kravs for att rotera, skala och translatera.
M = 3;
k = 1:3;
numstr_1 = {'G1','G2','G3'};
imshow(Gim);

[GX,GY] = ginput(M);
   
hold on;
text(GX(k),GY(k), numstr_1(k));
close

% Canon
GC = zeros(M,2);

GC(:,1) = GX;
GC(:,2) = GY;

GC1 = GC;
GC1(:,3) = ones(3,1);

% Holga
figure
imshow(Him);
numstr_2 = {'H1','H2','H3'};
[HX,HY] = ginput(M);

hold on;
text(HX(k),HY(k), numstr_2(k));
close

HC = zeros(M,2);

HC(:,1) = HX;
HC(:,2) = HY;

HC1 = HC;
HC1(:,3) = ones(3,1);

% GC1*A = HC1 --> mldivide to sovle A = HC1 / GC1

A = GC1 \ HC1  ;
%A = mldivide(HC1,GC1);

% Create a meshgrid
[X, Y] = meshgrid(1:500, 1:512);

% transform the coordinates
new = [X(:), Y(:), ones(256000,1) ]*A;

% interpolate
Vq = interp2(Him, new(:,1), new(:,2) );

% reshape back to image 
B = reshape(Vq, [512 500]);

% show the image
%imshow(B);

imshowpair(Gim, B);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Losning med 2 for loopar, som loopar igenom alla pixlar pa bilden.
% http://se.mathworks.com/help/images/performing-general-2-d-spatial-transformations.html#f12-31921

% tform = affine2d(A);
% cb_rgb = imwarp(Him,tform);
% 
% figure;
% imshow(cb_rgb)
% 
% imshowpair(Him, cb_rgb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% forsoker loopa med ett linjart index
% ska en kolumn vektor

% % Vektor med de orginal koordinater
% q = [0 0 0]';
% % Vektor med de nya transformerade koordinaterna
% p = [0 0 0]';
% % Tom bilder matris som ska innehalla den nya transformerade bilden
% result = Gim*0;

% for i = 1: 256000 %256000
%     
%     %konvertarar linjart index till subscript
%     [r,c] = ind2sub(size(Gim), i);
%    
%     % Nya icke-linjara koordinater
%     q(1) = r;
%     q(2) = c;
%     
%     % Multiplicera koordinaterna med matrisen A
%     p = A*q;
%     % Avrunda till heltal
%     p = round(p);
%      
%     % Vilkor for att transformerade koordinater ska ligga inom bilden
%     if (p(2) <= 0)
%         p(2) = 1; 
%     end
%     if (p(1) <= 0)
%         p(1) = 1;
%     end
%     
%     if ((p(2) > 512))
%         p(2) = 512;
%     end
%     if ((p(1) > 500))
%         p(1) = 500;
%     end
%     
%     % Konvertera tillbaka till linjara koordinater 
%     index = sub2ind(size(Gim), p(1), p(2));
%     
%     % Satt in motsvarande pixelvarde i den nya tranformerade bilder.
%     result(index) = Gim(i);
%     
% end

% imshow(result);
% figure;
% imshowpair(Him, result);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%