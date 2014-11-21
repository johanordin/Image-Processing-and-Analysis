%% Labb 1 - 
% - TNM087
%%%%%%%
%%%%%%%
%
clear all;


%% Dot 1
% Read in and plot the gfun

makegfun;
load('gfun.mat');
plot(gfun);
hold on
plot(2.^gfun)


%% Dot 2
%

%Kontroll över index
%[x,y] = ind2sub(size(bild1), index_M)

ImgCell = {'Img1.tiff', 'Img2.tiff','Img3.tiff','Img4.tiff','Img5.tiff','Img6.tiff','Img7.tiff','Img8.tiff','Img9.tiff','Img10.tiff','Img11.tiff','Img12.tiff','Img13.tiff','Img14.tiff'};

%%%%
%Spara bilderna i en 4-dim matris
imagesCell = {imread('Img1.tiff'), imread('Img2.tiff'),imread('Img3.tiff'),imread('Img4.tiff'),imread('Img5.tiff'),imread('Img6.tiff'),imread('Img7.tiff'),imread('Img8.tiff'),imread('Img9.tiff'),imread('Img10.tiff'),imread('Img11.tiff'),imread('Img12.tiff'),imread('Img13.tiff'),imread('Img14.tiff')};
pics4dimarray  = imread('Img1.tiff');

for i =2:14
    pics4dimarray(:,:,:,i) = imagesCell{i};
end

%%%
% maxvärde av en av de första bilderna
image_dark=imread('Img1.tiff');
M = max(max(image_dark(:,:,1)));
index_M = find(image_dark==M);

% minvärde av en av de sista bilderna
image_light=imread('Img14.tiff');
m = min(min(image_light(:,:,1)));
index_m = find(image_light==m,1);

% Medianvärde
image_median=imread('Img9.tiff');
med = median(median(image_median(:,:,1)));
index_med = find(image_median==med,1);

% Allokerar minne
arr1 = zeros(14,1);
arr2 = zeros(14,1);
arr3 = zeros(14,1);

for i=1:length(ImgCell)
    
    %temp = ImgCell{i};
    %bild = imread(temp);
    
    arr1(i) = imagesCell(index_M);
    arr2(i) = imagesCell(index_m);
    arr3(i) = imagesCell(index_med);
end


% Plottar de pixelvärderna för de olika bilderna
plot(arr1);
plot(arr2);
plot(arr3);



%% Dot 3

% Set some time constant
time = 1/128;

%E = zeros(size(bild,1)*size(bild,2),1);

% Read in the picture --> storage for the new pictures
newimagesCell = {imread('Img1.tiff'), imread('Img2.tiff'),imread('Img3.tiff'),imread('Img4.tiff'),imread('Img5.tiff'),imread('Img6.tiff'),imread('Img7.tiff'),imread('Img8.tiff'),imread('Img9.tiff'),imread('Img10.tiff'),imread('Img11.tiff'),imread('Img12.tiff'),imread('Img13.tiff'),imread('Img14.tiff')};
%sample = uint(zeros(683, 1024, 3));   %# Creates a 20x10x3 matrix
%newimagesCell = {sample,sample,sample,sample,sample,sample,sample,sample,sample,sample,sample,sample,sample,sample};


%%%
%index i --> picture in the picture array
%index j --> every pixel in the picture

for i = 1:14
    %temp = ImgCell{i};
    %bild = imread(temp);
    for j = 1:683*1024*1
        
        z = imagesCell{i}(j)+1;
        newimagesCell{i}(j) = gfun(z) - log2(time);
    end
       time = time*2;
end


imshow(newimagesCell{13})


%%

%imagesCell{i}(j)+1

for i = 1:14
    for j = 1:683*1024*2
        
        z = imagesCell{i}(j)+1;
        
        if z <= 349697
            newimagesCell{i}(j) = z - imagesCell{i}(1)+1;
        else
            newimagesCell{i}(j) = imagesCell{i}(end)+1 - z ;
        end
    end
       time = time*2;
end



%%

% % % % C = permute(B,2)
% % % % C = permute(B,1)
% % % % C = permute(B,1)C = zeros(4, 4)
% % % % C = zeros(4, 4)
% % % % B = cat(dim, A1, A2...)
% % % % B = cat(3, A, A, A)
% % % % B = cat(3, A, A, A);
% % % % B = cat(2, A, A, A);
% % % % 


%







