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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spara bilderna i en 4-dim matris
% Ladda in bilderna i en cell array
imagesCell = {imread('Img1.tiff'), imread('Img2.tiff'),imread('Img3.tiff'),imread('Img4.tiff'),imread('Img5.tiff'),imread('Img6.tiff'),imread('Img7.tiff'),imread('Img8.tiff'),imread('Img9.tiff'),imread('Img10.tiff'),imread('Img11.tiff'),imread('Img12.tiff'),imread('Img13.tiff'),imread('Img14.tiff')};
% Skapa en ny matris för den första bilden
pics4dimarray  = imread('Img1.tiff');

% Skapa den 4-dimensionella matrisen med alla bilderna
for i =2:14
    pics4dimarray(:,:,:,i) = imagesCell{i};
end

% Creating a copy of the 4dim array index with zeros
pics4dimarrayNew = pics4dimarray.*0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
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

%%%
%index i --> picture in the picture array
%index j --> every pixel in the picture

for i = 1:14
    for j = 1:683*1024*1
        z = imagesCell{i}(j)+1;
        newimagesCell{i}(j) = gfun(z) - log2(time);
    end
       time = time*2;
end


imshow(newimagesCell{13})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% % % %intensity in picture --> 0 - 255
% % % %gfun values maps --> 1 - 256 dvs 1 -> 0 och 256->255

time = 1/128;

for i=1:14
    for values= 0:255
        index = find(pics4dimarray(:,:,:,i) == values);
        pics4dimarrayNew(index) = gfun(values+1) - log2(time);
        time = time*2;
    end
end

% Plot the exposure images
for i=1:14
    imshow(pics4dimarrayNew(:,:,:,i))
    figure;
end


%% Sample of small 4-dim array
% % % %intensity in picture --> 0 - 255
% % % %gfun values maps --> 1 - 256 dvs 1 -> 0 och 256->255
% for values= 0:255
%     index = find(A(:) == values);
%     A(index) = gfun(values+1);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% Set some time constant
time = 1/128;
time2 = 1/128;


storage = zeros(14,1);
storage2 = zeros(14,1);
for x=1:14

    storage(x) = log2(time);
    storage2(x) = log2(time2);
    time = time*2;
end


plot(storage);
hold on ;
plot(storage2);
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Skapa en 3x3 matris
% A = [1, 2, 1;2 4 5; 3 5 6]
% % skapa en 3x3x2 matris (3-dim)
% A(:,:,2) = [1 2 2;7 6 5; 2 5 2]
% % hitta värdet 7
% index = find(A(:,:,:) == 7)
% A(index) = 16
% 
% index = find(A(:) == 1)
% A(:,:,3) = [1 2 2;7 6 5; 2 5 2]
% 
% %Skapa en 3x3x3x2 matris (4-dim)
% A(:,:,:,2) = A
% index = find(A(:) == 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






