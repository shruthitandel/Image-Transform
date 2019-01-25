%1. Affine Image Transformation
%1.1 Resampling

clc;
close all;
clear all;

%Reading source image
src_img1=imread('edges-lines-orig.png');
src_img = im2double(src_img1);
figure;
imshow(src_img1);
title('Original Image');
axis on;

%Initialize target image
target = zeros([size(src_img1,1),size(src_img1,2)]);
%target1 = zeros([size(src_img1,1),size(src_img1,2)]);

%Translation
tx= input('Enter Tx for Translation');
ty= input('Enter Ty for Translation');
transfo_mat4 = [1,0,tx;0,1,ty;0,0,1];
% 
% % Interpolation using Nearest Neighbour 
% translated_nn = nearest(src_img,target,transfo_mat1);
% subplot(1,3,2);
% imshow(translated_nn);
% title('Translated Image using Nearest Neighbor');
% axis on;
% 
% 
% % Interpolation using Bilinear Interpolation
% translated_biln = bilinear(src_img,target,transfo_mat1);
% subplot(1,3,3);
% imshow(translated_biln);
% title('Translated Image using Bilinear Interpolation');
% axis on;

%Scaling
  degree=input('enter the rotating angle   ');
  transfo_mat = [cos(deg2rad(degree)),-sin(deg2rad(degree)),0;sin(deg2rad(degree)),cos(deg2rad(degree)),0;0,0,1];
  sx= input('Enter Sx for Scaling');
  sy= input('Enter Sy for Scaling');
  transfo_mat2 = [sx,0,0;0,sy,0;0,0,1];
  shx= input('Enter Sx for Shearing');
  shy= input('Enter Sy for Shearing');
  transfo_mat3 = [1,0,shy;shx,1,0;0,0,1];
  transfo_mat1 = transfo_mat4 * transfo_mat2 * transfo_mat;

% Interpolation using Nearest Neighbour 
scaled_nn = nearest(src_img,target,transfo_mat1);
figure;
imshow(scaled_nn);
title('Rotated,scaled and translated using Nearest Neighbor');
axis on;
imwrite(scaled_nn,'scaled_nn.png');


% Interpolation using Bilinear Interpolation
scaled_biln = bilinear(src_img,target,transfo_mat1);
figure;
imshow(scaled_biln);
title('Rotated,scaled and translated using Bilinear Interpolation');
axis on;
imwrite(scaled_biln,'scaled_biln.png');




 function[result_img] = nearest(src_img,target,transfo_mat1)
    result_img = zeros([size(src_img,1),size(src_img,2)]);
    for i=1:size(target,1)
        for j=1:size(target,2)
            current_v=[i;j;1];
            new_v= mtimes(inv(transfo_mat1),current_v);
            if(new_v(1)>=1 && new_v(2)>=1 && new_v(1)<=size(src_img,1) && new_v(2)<=size(src_img,2))
                intensity = src_img(round(new_v(1)),round(new_v(2)));
                result_img(i,j) = intensity;
            end
        end
    end
end

function[result2_img] = bilinear(src_img,target,transfo_mat1)
    result2_img = zeros([size(src_img,1),size(src_img,2)]);
    for i=1:size(target,1)
        for j=1:size(target,2)
            current_v=[i;j;1];
            new_v= mtimes(inv(transfo_mat1),current_v);
             if(new_v(1)>=1 && new_v(2)>=1 && new_v(1)<=size(src_img,1) && new_v(2)<=size(src_img,2))
                  x1 =  floor(new_v(1));
                  x2 = ceil(new_v(1));
                  y1 =  floor(new_v(2));
                  y2 = ceil(new_v(2));
                   if x1==x2 
                       if (x2 + 1)<=size(src_img,1)
                         x2= x2 + 1;
                       end
                   end
                   if y1 == y2
                       if (y2 + 1)<=size(src_img,2)
                        y2= y2 + 1;
                       end
                   end
                  
             neighbor=[x1, x2, y1, y2];
             
             % Compute coefficients
              b1 = src_img(neighbor(1),neighbor(3));
              b2 = src_img(neighbor(2),neighbor(3));
              b3 = src_img(neighbor(1),neighbor(4));
              b4 = src_img(neighbor(2),neighbor(4));

              % compute new intensity
              newint = (b1*(neighbor(2)-new_v(1))*(neighbor(4)-new_v(2))) ...
              +(b2*(new_v(1)-neighbor(1))*(neighbor(4)-new_v(2))) ...
              +(b3*(neighbor(2)-new_v(1))*(new_v(2)-neighbor(3))) ...
              +(b4*(new_v(1)-neighbor(1))*(new_v(2)-neighbor(3))) ...
              /((neighbor(2)-neighbor(1))*(neighbor(4)-neighbor(3)));
              result2_img(i,j) = newint;
            end
        end
    end
end


