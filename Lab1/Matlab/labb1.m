%% TNM087 - Labb 1 
% Christoffer Engelbrektsson & Johan Nordin
% MT3



%% Dot 1
% Read in and plot the gfun

clear all;
makegfun;
load('gfun.mat');
plot(gfun);
hold on
plot(2.^gfun)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spara bilderna i en 4-dim matris

load('gfun.mat');
% Ladda in bilderna i en cell array
imagesCell = {imread('Img1.tiff'), imread('Img2.tiff'),imread('Img3.tiff'),imread('Img4.tiff'),imread('Img5.tiff'),imread('Img6.tiff'),imread('Img7.tiff'),imread('Img8.tiff'),imread('Img9.tiff'),imread('Img10.tiff'),imread('Img11.tiff'),imread('Img12.tiff'),imread('Img13.tiff'),imread('Img14.tiff')};

% Skapa en ny matris fr den frsta bilden fr att alokera minne
 pics4dimarray  = im2double(imread('Img1.tiff'));

% Skapa den 4-dimensionella matrisen med resterande bilderna
for i = 2:14
    pics4dimarray(:,:,:,i) = im2double(imagesCell{i});
end
 
% Creating a copy of the 4dim array index with zeros
pics4dimarrayNew = pics4dimarray.*0;
weightfunc = pics4dimarray.*0;
finalMatrix = pics4dimarray*0;
finalWeight = pics4dimarray*0;

montage(pics4dimarray);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% "Preperations assigments"

% maxvärde av en av de första bilderna
image_dark=imread('Img1.tiff');
M = max(max(image_dark(:,:,1)));
index_M = find(image_dark == M);

% minvärde av en av de sista bilderna
image_light=imread('Img14.tiff');
m = min(min(image_light(:,:,1)));
index_m = find(image_light == m,1);

% Medianvärde
image_median=imread('Img9.tiff');
med = median(median(image_median(:,:,1)));
index_med = find(image_median == med,1);

% Allokerar minne
arr1 = zeros(14,1);
arr2 = zeros(14,1);
arr3 = zeros(14,1);

for i = 1:length(ImgCell)
    
    arr1(i) = imagesCell(index_M);
    arr2(i) = imagesCell(index_m);
    arr3(i) = imagesCell(index_med);
    
end

% Plottar de pixelvärderna frår de olika bilderna
plot(arr1);
plot(arr2);
plot(arr3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Logaritm images to exposure images

for i=1:3
    for values = 0:255
        index = 255.*pics4dimarray(:,:,i,:) == values;
        pics4dimarrayNew(index) = gfun(values+1,i);
    end
end

time = 1;

for i=1:14
    pics4dimarrayNew(:,:,:,i) = pics4dimarrayNew(:,:,:,i) - log2(time);
    time = time*2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Weightfunction

maximum = 1;
minimum = 0;

condition = (maximum+minimum)/2;

d = (pics4dimarray <= condition);

w1 = d.* (pics4dimarray);

e = (pics4dimarray > condition);

w2 = e.* (maximum - pics4dimarray );

weightfunc = w1+w2;

montage(weightfunc)
        
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Summerar alla bilder till en HDR-bild

finalMatrix = sum(weightfunc.*pics4dimarrayNew, 4);

finalWeight = sum(weightfunc, 4);

finalMatrix = 2.^((finalMatrix) ./ finalWeight);

imshow(tonemap(finalMatrix));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sample of small 4-dim array
% % % %intensity in picture --> 0 - 255
% % % %gfun values maps --> 1 - 256 dvs 1 -> 0 och 256->255
% for values= 0:255
%     index = find(A(:) == values);
%     A(index) = gfun(values+1);
% end