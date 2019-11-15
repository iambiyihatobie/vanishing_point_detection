%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Vanishing Point Estimation from Image %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

%% Import Tset Image
Iorig = imread('Sample_1_Buildings.jpg');
figure; imshow(Iorig, [])
title('Original Image')

%% Important Parameters
N = 500;          % Number of Pixels added to the image
CT = 0.6;         % Canny Edge Detection Threshold
HT = 60;          % Line length

%% Surround the image with white pixels N = number of pixels added
I = im2double(Iorig);
[nrows,ncols,~] = size(I);
Im = ones(nrows+N,ncols+N,3);
[n1rows,n1cols,~] = size(Im);
r = (n1rows - nrows) / 2;
L = N - r + 1;
Im(r:end-L, r:end-L,:) = I;

%% Edge Detection
Ig = rgb2gray(I);
BW = edge(Ig,'Canny', CT);
BW = padarray(BW, [r, r]);
%figure; imshow(BW, []);

%% Hough Transform
[H,theta,rho] = hough(BW);

%% Hough transform matrix peak
P = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));

%% Hough LInes
lines = houghlines(BW,theta,rho,P,'FillGap',5,'MinLength',HT);

%% Drawing the lines using Hough LInes
figure, imshow(Im), hold on
max_len = 0;
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   slope = @(line) (line(2,2) - line(1,2))/(line(2,1) - line(1,1));
   intercept = @(line,m) line(1,2) - m*line(1,1);
   
   m = slope(xy);
   b = intercept(xy, m);
   x = linspace(0, length(BW));
   y = m*x + b;
   
   %%%plot extended lines
   plot(x,y,'LineWidth',1.5,'Color','red'); 
   
   %%%Plot beginnings and ends of lines of original lines
   % plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','red'); %view the real lines
   % plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   % plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','blue');
end

%% Intersection and Vanishing Point Estimation
interX= [];
for i = 1:length(lines)-1
  for j = i+1:length(lines)
    xy1 = [lines(i).point1; lines(i).point2];
    xy2 = [lines(j).point1; lines(j).point2];
    m1 = slope(xy1);
    m2 = slope(xy2);
    b1 = intercept(xy1,m1);
    b2 = intercept(xy2,m2);
    xintersect = (b2-b1)/(m1-m2);
    yintersect = m1*xintersect + b1;
    if lines(i).theta ~= lines(j).theta
        interX = [interX;[xintersect yintersect]];
        %plot(xintersect, yintersect,'g.','MarkerSize', 50)
    end
  end
end

% Vanishing Point Estimation
vanishingPoint = mostfreq(interX);
interX_rd = round(interX, -1);
plot(vanishingPoint(1), vanishingPoint(2), 'g.','MarkerSize', 50)
title('Vanishing Point Detection')
disp("Vanishing Point Coordinate")
disp(vanishingPoint)
hold off

%% Voting Function [ Intersection point with most occurence get sorted out]
function t = mostfreq(m)
m = round(m, -1);     % roundup the intersection points
[Au,~,ic] = unique(m, 'rows');
tally = accumarray(ic, 1);
n = max(tally);

for i=1:length(tally)
    if tally(i)== n
        index = i;
        break;
    end
end
t = Au(index,:);
end
%%