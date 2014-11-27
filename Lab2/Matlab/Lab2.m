

%%Lab 2: Operations in the image domain

%%Part 1: Vignette and uneven illumination
%% A



pictures(:,:,:,1) = im2double(imread('CWhite1.jpg'));

pictures(:,:,:,2) = im2double(imread('HWhite1.jpg'));

pictures = imresize(pictures, 0.5);

N = 512;

[X,Y] = meshgrid((1:N));
[T,R] = cart2pol(X-N/2,Y-N/2);


%Plot the vector R(N/2,:) and use imshow(R,[])

subplot(1,2,1)
plot(R(N/2,:));
subplot(1,2,2);
imshow(R, []);

%% B

RS = R;

%for i=1:512*512
    RS = R ./ R(512/2-1,1);
%end
