%hough transform
clc;
close all;

img_1=imread('HT.jpg');
%img_1 = rgb2gray(img_1);
subplot(2,3,1);
imshow(img_1);
title('Input image');
input_img = im2double(img_1);
[m,n]=size(input_img); %size of input image
z=ones(3); %mask size is 3x3
[p,q]=size(z);
smooth = zeros(m,n);

%Padding with zeroes
w=1:p;
x=round(median(w));
anz=zeros(m+2*(x-1),n+2*(x-1));


for i=x:(m+(x-1))
    for j=x:(n+(x-1))
        anz(i,j)=img_1(i-(x-1),j-(x-1));
    end
end


% smoothing of image 
sum=0;
x=0;
y=0;
for i=1:m
    for j=1:n
        for k=1:p
            for l=1:q 
                sum= sum+anz(i+x,j+y)*z(k,l);
                y=y+1;
            end
            y=0;
            x=x+1;
        end
        x=0;
        smooth(i,j)=(1/(p*q))*(sum);
        sum=0;
    end
end
subplot(2,3,2);
imshow(uint8(smooth));
title('smoothed image');

C=double(smooth);
treshold1 = [0.06,0.5];
edge1 = edge(img_1,'canny',treshold1);

% % Edge detection
%   dx = zeros(size(C,1),size(C,2));
%   dy = zeros(size(C,1),size(C,2));
% 
% 3x1 i.e. [-1;0;1] mask for x-derivative:
% for i=1:size(C,1)
%     for j=2:size(C,2)-1
%           dx(i,j) = (C(i,j+1)- C(i,j-1))/2;
%     end
% end
% 
% 1x3 i.e. [-1,0,1] mask for y-derivative:
% for i=2:size(C,1)-1
%     for j=1:size(C,2)
%           dy(i,j) = (C(i+1,j)- C(i-1,j))/2;
%     end
% end
% 
% The magnitude of the gradient of the image
% edge =sqrt(dx.^2+dy.^2);
% Direction of the gradient of image 
%  dir = atan2(dy,dx);
%  
%  edge1 = uint8(255 * mat2gray(edge));
% [H,T,R] = hough(edge1);
% subplot(2,3,6);
% imshow(uint8(255 * mat2gray(H)),[],'XData',T,'YData',R,...
%             'InitialMagnification','fit');
% xlabel('\theta'), ylabel('\rho');
% axis on, axis normal, hold on;
% rgb =im2uint8(imoverlay(smooth,dx,[1 0 0]));
% subplot(2,3,2);
% imshow(uint8(255 * mat2gray(dx))); 
% title('x-derivative');
% 
% subplot(2,3,3);
% imshow(uint8(255 * mat2gray(dy))); 
% title('y-derivative');
% 
subplot(2,3,3);
imshow(uint8(255 * mat2gray(edge1)));
title('Edge Magnitude');
% edge1 = uint8(255 * mat2gray(edge));
% subplot(2,3,6);
% imshow(uint8(255 * mat2gray(dir)));
% title('Gradient Direction');

%HT
% set up variables for hough transform
theta_frequency = 0.01;                                             
[x, y] = size(edge1);
rho_limit = norm([x y]);                                                
rho = (-rho_limit:1:rho_limit);
theta = ((-pi/2):theta_frequency:(pi/2));
num_of_thetas = numel(theta);
num_of_rhos = numel(rho);
hough_space = zeros(num_of_rhos, num_of_thetas);


%perform hough transform
for xi = 1:x
    for yj = 1:y
        if edge1(xi, yj) ==1
            for theta_index = 1:num_of_thetas
                th = theta(theta_index);
                r  = xi * cos(th) + yj * sin(th);
                rho_index = round(r + num_of_rhos/2);                      
                hough_space(rho_index, theta_index) = hough_space(rho_index, theta_index) + 1;
            end
        end
    end
end  

%show hough transform
% mat2gray(), which converts any matrix to a 0-1 range.
%Then multiply by 255 and cast to uint8 if you want a uint8 data type.
subplot(2,3,4);
imagesc(theta, rho, hough_space);
title('Peaks detected on Hough Transform');
xlabel('Theta (radians)');
ylabel('Rho (pixels)');
colormap(gca,hot);


%detect peaks in hough transform 
r_peak = [];
c_peak = [];
%max(hough_space): is a row vector containing the maximum element from each column
%     It returns the a row vector 'max_in_col' that contains the maximum value 
%     of each column of 'hough_space', 
%     'row_number' is a row vector containing the row positions of maximum

[max_in_col, row_number] = max(hough_space);
[rows, cols] = size(edge1);
mxval = max(max(hough_space));
thresh = mxval/2;
%After peak detection write the parameters (theta,rho,votes) Mymatrix.txt file
file_id = fopen('Mymatrix.txt','wt'); 
for i = 1:size(max_in_col, 2)
   if max_in_col(i) > thresh
       c_peak(end + 1) = i;
       r_peak(end + 1) = row_number(i);
       fprintf(file_id,'%g\t',theta(i));
        fprintf(file_id,'%g\t',rho(row_number(i)));
        fprintf(file_id,'%g\t',max_in_col(i));
        fprintf(file_id,'\n');
       
   end
end
fclose(file_id);

% plot all the detected peaks on hough transform image
hold on;
plot(theta(c_peak), rho(r_peak),'gx');
hold off;



% plot the detected line superimposed on the original image
subplot(2,3,5)
imagesc(input_img);
title('Straight lines highlighted in original image');
colormap(gray);
hold on;

for i = 1:size(c_peak,2)
    th = theta(c_peak(i));
    rh = rho(r_peak(i));
    m = -(cos(th)/sin(th));
    b = rh/sin(th);
    x = 1:cols;
    hold on;
    plot(m*x+b, x,'g');
    hold off;
end