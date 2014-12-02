

%%Lab 2: Operations in the image domain

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

SR = R;
SR = R ./ R(512/2-1,1);

QR = round((SR*100));

average = zeros(141,2);
norm_average = zeros(141,2);
 
for i=1:2
    for m=1:100
        % Create logical matrix Maskm for each value m
        Maskm = QR == m;

        % Sum the pixelvalues for each m in QR
        sum_pv = sum(sum(Maskm.*pictures(:,:,:,i)));

        % Count the number of elements of m in QR // samma sak som sum(sum(QR
        % == 1)) //
        nr_objects = sum(sum(Maskm));

        average(m,i) = sum_pv/nr_objects;
    end
    norm_average(:,i) = average(:,i) / max(average(:,i));
end

subplot(2,1,1);
plot(norm_average(:,1));
subplot(2,1,2);
plot(norm_average(:,2));

%pictures(:,:,:,2)


%% Part 2: Correlation 

redEyes = im2double(imread('BoldRedEye.JPG'));

% Extract the read channel
redChannel = redEyes(:,:,1);

% Load in the mask
load('RedEyeMask');

sq_filter_32 = ones(32);

MFilterImage = imfilter(redChannel, sq_filter_32);
EyeFilterImage = imfilter(redChannel, RedEyeMask);

ratio = EyeFilterImage./MFilterImage;

% Normalize the ratio matrix
imshow(ratio / max(max(ratio)))
figure;

quant = quantile(quantile(ratio, 0.90), 0.90);

ratio = ratio >= quant;

imshow(ratio);
figure;

% return two-dimensional eight-connected neighborhood 
BW = ratio.*imregionalmax(ratio);

imshow(BW)


%% Part 3: Registration
































