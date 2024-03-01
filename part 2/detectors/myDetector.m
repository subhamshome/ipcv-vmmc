%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VIDEO FOR MULTIPLE AND MOVING CAMERAS (VMMC)
% IPCV @ UAM
% Marcos Escudero-Viñolo (marcos.escudero@uam.es)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [points,decomposition] =  myDetector(gima,params)
    %% get detections:
    switch params.detector

         case 'LoG_ss'

              SS                  =   doScaleSpaceGivenSigmas(gima,params.sigmas);
              LoG                 = extractLaplacianGivenSigmas(SS,params.sigmas);
              [M,MaxPos,m,MinPos] = MinimaMaxima3D(LoG,1,0,params.npoints,params.npoints);

    
              points.Metric          = abs(m);
              points.Location        = [MinPos(:,2),MinPos(:,1)];
              points.SignOfLaplacian = -1;

              points.sid             = MinPos(:,3);

              points.Metric          = [points.Metric;M];
              points.Location        = [points.Location;[MaxPos(:,2),MaxPos(:,1)]];
              points.SignOfLaplacian = 1;
              points.sid             = [points.sid;MaxPos(:,3)];
              
              points.Scale           = params.sigmas(points.sid);
              
              decomposition          = SS;
        case 'DoH'

              SS                  = doScaleSpaceGivenSigmas(gima,params.sigmas);
              LoG                 = extractDeterminantOfHessianGivenSigmas(SS,params.sigmas);
              [M,MaxPos,m,MinPos] = MinimaMaxima3D(LoG,1,0,params.npoints,params.npoints);

    
              points.Metric          = abs(m);
              points.Location        = [MinPos(:,2),MinPos(:,1)];
              points.SignOfLaplacian = -1;

              points.sid             = MinPos(:,3);

              points.Metric          = [points.Metric;M];
              points.Location        = [points.Location;[MaxPos(:,2),MaxPos(:,1)]];
              points.SignOfLaplacian = 1;
              points.sid             = [points.sid;MaxPos(:,3)];
              
              points.Scale           = params.sigmas(points.sid);
              
              decomposition          = SS;
              
         case 'SURF'
             
             points                 =  detectSURFFeatures(gima);
             
             decomposition{1}       = gima;
            
        case 'KAZE'
             
             points                 =  detectKAZEFeatures(gima);
             
             decomposition{1}       = gima;
             
        case 'SIFT'
             % see additional detection methods and options provided by vlfeat in help vl_covdet
              [frame]= ...
              vl_sift(single(gima),'Octaves',params.noctaves,'Levels',params.nscales); 
            % parse to common structure
               points.Location        = [frame(1,:)',frame(2,:)'];
               points.Scale           =  frame(3,:)';
               points.Orientation     =  frame(4,:)';
               points.Metric          =   0.5.*ones(size( frame(4,:)'));% unknown
               points.SignOfLaplacian =        zeros(size(frame(4,:)'));     % unknown
               
               decomposition{1}       = gima;    

         case 'HARRIS'
             points = detectHarrisFeatures(gima);
             decomposition{1}       = gima; 
               ...

        otherwise  % SURF

              points                 =  detectSURFFeatures(gima);
             
              decomposition{1}       = gima;

    end
% note that matlab's functions: detectKAZEFeatures, detectSURFFeatures,
% detect<<XXX>>Features,.. allow you to detect using several of the
% handcrafted detectors studied in class.
end