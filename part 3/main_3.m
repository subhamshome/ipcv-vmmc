%% Part-3
clc;
clear;
close all;
warning off;

ACT_path = 'ACT_lite/';
addpath(genpath(ACT_path));
% include extra funs
extra_funs_path = 'extra_funs/';
addpath(genpath(extra_funs_path));
allfns = 'allfns/';
addpath(genpath(allfns));

%% Loading the features
part2_data = "part2_data.mat";
features=load(part2_data).features;
ima=load(part2_data).ima;
points=load(part2_data).points;
load('part1_data.mat');

%%Matching parameters (see doc matchFeatures for explanation and additional parameters)
MaxRatio     =  0.6;  
Metric       =  'SSD';
 
n=(3:10); % (3:10); (1:8); (5:12);
ncam=length(n);


% extracting selected features.
for i=1:ncam
    features_{i}=features{n(i)};
    points_{i}=points{n(i)};
    ima_{i}=ima{n(i)};
end

%% 

%generating the data points.
q_data = n_view_matching(points_, features_, ima_, MaxRatio, Metric);
q_data = homogenize_coords(q_data);

%Scaling factor
imscale = 0.5;

disp('************************************* START *************************************');

K = A_ref;
K_est=K.*imscale; K_est(3,3)=1;

K_=zeros(3, 3, 2);
K_(:, :, 1)=K_est;
K_(:, :, 2)=K_est;


% ------------------------------------------------------------------------
% 2. Compute the fundamental matrix using the first and last cameras
% of the camera set (N cameras)
% ------------------------------------------------------------------------
% initial reconstruction

q_2cams(:,:,1)=q_data(:,:,1); 
q_2cams(:,:,2)=q_data(:,:,ncam);

[F, P_2cam_est,Q_2cam_est,q_2cam_est] = MatFunProjectiveCalib(q_2cams);

% for part---2
disp(['Residual reprojection error. 8 point algorithm   = ' num2str( ErrorRetroproy(q_2cams,P_2cam_est,Q_2cam_est)/2 )]);
draw_reproj_error(q_2cams,P_2cam_est,Q_2cam_est);

%%
% ------------------------------------------------------------------------
% 3. Resection. Obtain the projection matrices of the rest of cameras using the PDLT_NA function 
% ------------------------------------------------------------------------
% ...

P_cams(:,:,:)=zeros(3,4,ncam);
P_cams(:,:,1)=P_2cam_est(:,:,1);
P_cams(:,:,ncam)=P_2cam_est(:,:,2);

for i= 2:ncam-1
    P_cams(:,:,i)=PDLT_NA(q_data(:,:,i),Q_2cam_est);
end


% ------------------------------------------------------------------------
% 4. Compute the statistics of the reprojection error for the initial projective reconstruction
% ------------------------------------------------------------------------
disp(['Residual reprojection error, After resectioning  = ' num2str( ErrorRetroproy(q_data,P_cams,Q_2cam_est)/2 )]);
draw_reproj_error(q_data,P_cams,Q_2cam_est); %q_data

% ------------------------------------------------------------------------
% 5a. Projective Bundle Adjustment. Use BAProjectiveCalib function
% Coordinates of 3D and 2D points are given in homogeneus coordinates
% ------------------------------------------------------------------------
% auxiliary matrix that indicates that all points are visible in all the cameras

npoints=size(q_data,2);
vp = ones(npoints,ncam);
% ...
[P_bundle,Q_bundle,q_bundle]=BAProjectiveCalib(q_data,P_cams,Q_2cam_est,vp);

% ------------------------------------------------------------------------
% 5b. Compute the statistics of the reprojection error for the improved projective reconstruction
% ------------------------------------------------------------------------
% ...
disp(['Reprojection error, After Bundle Adjustment  = ' num2str( ErrorRetroproy(q_data,P_bundle,Q_bundle)/2 )]); %q_data
draw_reproj_error(q_data,P_bundle,Q_bundle);


% 6. Obtain the essential matrix (E) from the fundamental matrix (F) and the
% intrinsic parameter matrices (K).
% ------------------------------------------------------------------------
%E=K1'FK;

% estimating F from P_BA (P_bundle)
F_bundle=vgg_F_from_P({P_bundle(:,:,1),P_bundle(:,:,end)});

%estimating External matrix E
E = K_(:,:,1)'*F_bundle*K_(:,:,2);

% calculate the T & R matrices
[R_est,T_est]=factorize_E(E);

q_bundle_2cam(:,:,1)=q_bundle(:,:,1);
q_bundle_2cam(:,:,2)=q_bundle(:,:,end);

% ------------------------------------------------------------------------
% Save the 4 solutions (R,t) in the structures Rcam(3,3,cam,sol), T(3,cam,sol),
% where cam indicates the camera number and sol indicates the solution number (1, 2, 3 or 4).
% ------------------------------------------------------------------------
Rcam = zeros(3,3,2,4);
Tcam = zeros(3,2,4);

for i=1:size(Rcam,4)
    Rcam(:,:,1,i)=eye(3,3);
    if i<=2
        Rcam(:,:,2,i) = R_est(:,:,1);
    else
        Rcam(:,:,2,i) = R_est(:,:,2);
    end
end


Tcam(:,1,:)=0;
for i=1:size(Tcam,3)
    if mod(i,2) ~= 0
        Tcam(:,2,i) = T_est;
    else
        Tcam(:,2,i) = -T_est;
    end
end


% ------------------------------------------------------------------------
% 5. For each solution we obtain an Euclidean solution and we visualize it.
% ------------------------------------------------------------------------
%npoints = size(q_data,2);
Q_euc = zeros(4,npoints,2); % Variable for recontructed points
P_euc = zeros(3,4,2);       % Variable for projection matrices
figNo=figure;


%%
for sol=1:4
    % Euclidean triangulation to obtain the 3D points (use TriangEuc)
     Q_euc = TriangEuc(Rcam(:,:,2,sol),Tcam(:,2,sol),K_,q_bundle_2cam);
       
    % visualize 3D reconstruction
    figure();
    draw_scene(Q_euc, K_, Rcam(:,:,:,sol), Tcam(:,:,sol));
    title(sprintf('Solution %d', sol));
     
    % Compute the projection matrices from K, Rcam, Tcam
    for k=1:2
        rcam=Rcam(:,:,k,sol);
        tcam=Tcam(:,k,sol);
        P_euc(:,:,k) = K_(:,:,k)*[rcam, -1.*rcam*tcam] ; 
        %P_euc(:,:,k) = K(:,:,k)*[rcam -1*rcam*tcam] ;
    end
    
    % Obtain the re-projected points q_rep
     
     for j=1:2
        q_rep(:,:,j) = P_euc(:,:,j)*Q_euc;
    %Q_euc
     end
    
    % Visualize reprojectd points to check that all solutions correspond to
    % the projected images
    q_rep = un_homogenize_coords(q_rep);
    for k=1:2
      figure(figNo); 
      subplot(4,2,2*(sol-1)+k); scatter(q_rep(1,:,k),q_rep(2,:,k),30,[1,0,0]);
      title(sprintf('Reprojection %d, image %d', sol, k));
      daspect([1, 1, 1]);
      pbaspect([1, 1, 1]);
      axis([-1000, 1000, -1000, 1000]);
    end
end

disp(['Reprojection error, Euclidean Re-construction  = ' num2str( ErrorRetroproy(q_bundle_2cam,P_euc,Q_euc)/2 )]); %q_data
draw_reproj_error(q_bundle_2cam,P_euc,Q_euc);

disp('************************************* END *************************************');

save('final_results.mat');