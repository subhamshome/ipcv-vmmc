clc;
close all;
clear;
warning off;
%% Reading user points once marked and saved
load userpoints_part1_image1.mat;

%% Own camera calibration
num_points = 9; % Total number of points
width_checkerboard = 90;
[coords, ima_pattern]= get_real_points_checkerboard_vmmc(num_points, width_checkerboard, 1);

%% Reading the images
im1 = imread('input_images\1080p\image1.jpg');
im2 = imread('input_images\1080p\image2.jpg');
im3 = imread('input_images\1080p\image3.jpg');
im4 = imread('input_images\1080p\image4.jpg');
im5 = imread('input_images\1080p\image5.jpg');
im6 = imread('input_images\1080p\image6.jpg');
im7 = imread('input_images\1080p\image7.jpg');
im8 = imread('input_images\1080p\image8.jpg');
im9 = imread('input_images\1080p\image9.jpg');

% points1 = get_user_points_vmmc(im1); 
% points2 = get_user_points_vmmc(im2);
% points3 = get_user_points_vmmc(im3);
% points4 = get_user_points_vmmc(im4);
% points5 = get_user_points_vmmc(im5);
% points6 = get_user_points_vmmc(im6);
% points7 = get_user_points_vmmc(im7);
% points8 = get_user_points_vmmc(im8);
% points9 = get_user_points_vmmc(im9);

%% Getting the Homographies
figure;
sgtitle("Homography with Zhang's Method");
set(gcf, 'WindowState', 'maximized');

for i = 1:num_points

    homo = homography_solve_vmmc(coords', eval(['points' num2str(i) '(:,1:9)']));
    [homo_ref, r] = homography_refine_vmmc(coords', eval(['points' num2str(i) '(:,1:9)']), homo);
    H_ref{i} = homo_ref;
    tform = maketform('projective', homo_ref');
    J_homo = imtransform(ima_pattern, tform, 'XData', ...
        [1 size(eval(['im' num2str(i)]), 2)], 'YData', [1 size(eval(['im' num2str(i)]), 1)]);

    % Subplots
    subplot(2, num_points, i);
    imshow(eval(['im' num2str(i)]));
    subplot(2, num_points, i + num_points);
    imshow(J_homo);
end

%% Calculation of internal parameters 

A_ref = internal_parameters_solve_vmmc(H_ref);
[m, n] = size(A_ref);

disp("Internal Parameters Matrix: ");
for i = 1:m
    for j = 1:n
        fprintf('%f     ',A_ref(i, j));
    end
    fprintf('\n');
end

