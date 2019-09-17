%
%	CT project
%   
%   This template is provided to guide you through the CT project
%   Find and replace all question marks (???)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   CREATE DISK PHANTOM AND SINOGRAM
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Parameters for the disk phantom are
%	        x-center, y-center, radius, attenuation coefficient
%   The data is stored in MATLAB variable "phantom"
%   The sinogram is stored in corresponding MATLAB variable named "sg1" 
%
circ = [  0   0 110  2; 
        -65   0  20  1; 
          0   0  35  0; 
         65 -25  25  4
         50  50  7   8];
%
%	Image parameters: number of pixels, size, etc.
%
nx = 128; ny = 128;
dx = 2;		                        % 2 mm / pixel
x = dx * ([1:nx]'-(nx+1)/2);
y = -dx * ([1:ny]'-(ny+1)/2);
xx = x(:,ones(1,ny));
yy = y(:,ones(1,ny))';
%
%	Generate data for disk phantom
%
  phantom = zeros(nx,ny);
  for ii=1:size(circ,1)
    cx = circ(ii,1); cy = circ(ii,2); rad = circ(ii,3); amp = circ(ii,4);
    t = find( ((xx-cx)/rad).^2 + ((yy-cy)/rad).^2 <= 1 );
    phantom(t) = amp * ones(size(t));
  end
%
% Image the phantom
%
  figure(1)
    imagesc(x, y, phantom')               % NOTE the transpose (') here and the x and y values
    colormap('gray')
    axis('square')
    title('Disk Phantom')
    xlabel('Position')
    ylabel('Position')
%
%	Geometry parameters
%
nr = 128;	dr = 2;		            % number of radial samples and ray spacing
na = nr*2;          	            % number of angular samples
r = dr * ([1:nr]'-(nr+1)/2);	    % radial sample positions
angle = [0:(na-1)]'/na * pi;	    % angular sample positions
%
%	Compute sinogram for the phantom
%
     rr = r(:,ones(1,na));
     sg1 = zeros(nr, na);
  for ii=1:size(circ,1)
    cx = circ(ii,1); cy = circ(ii,2); rad = circ(ii,3); amp = circ(ii,4);
    tau = cx * cos(angle) + cy * sin(angle);
    tau = tau(:,ones(1,nr))';
    t = find( (rr-tau).^2 <= rad.^2 );
    if ii > 1, amp = amp - circ(1,4);, end	% small disks embedded
    sg1(t) = sg1(t)+amp*2*sqrt(rad^2-(rr(t)-tau(t)).^2);
  end
%
% Sinogram of the phantom
%
  figure(2)
    imagesc(r, angle/pi*180, sg1')   % NOTE the transpose (') here and
    colormap('gray')                 % the fact that angle is displayed in degrees
    title('Sinogram: Disk Phantom')
    xlabel('Position (i.e., Rays)')
    ylabel('Angle (i.e., Views)')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   YOUR MODIFICATIONS BEGIN HERE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	Let's make a common variable (sinogram) so that your code is not linked
%	to any specific sinogram - later you should be able to cut and paste
%	your code and reconstruct an unknown object
%
sinogram = sg1;                         % disk phantom sinogram
%
%	Since different size sinograms are used
%
disp(sprintf('number of rays = %g', nr))
disp(sprintf('number of views = %g', na))
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   QUESTION A - COMPUTE THE 0th MOMENT (a.k.a. AREA UNDER THE CURVE) 
%                AS A FUNCTION OF THE PROJECTION ANGLE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
 disp('Question A: computing 0th moment')
 theta = 45;
%  Converting angle to y-axis index
[~, angle_index] = min(abs(angle-deg2rad(theta)));

projection = sg1(:, angle_index-1);
%
% Plot Oth moment of projections
%
figure(3)
   projection_max=max(projection);
   plot(projection,'k');
   axis([0,na,0,2*projection_max]);
%    Remove 'interpreter', 'latex' if computer does not support latex
   title('\textbf{$0^{th}$ moment of each projection at $45^\circ$}', 'interpreter', 'latex')
   xlabel('Position (i.e., Rays)')
   ylabel('Area Under the curve')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   QUESTION B - IMPLEMENT SIMPLE BACKPROJECTION, i.e. 
%                PRODUCE LAMINOGRAMS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
 disp('Question B: simple backprojection')
if 1~=exist('lamin')  % Check for variable ..., if exist - do not compute it gain
                      % You may want to comment out these "if-end" construct if you like  
%  
  lamin = zeros(nx,ny);
  mid_pt = nx/2;
  for ia = 1:na
    disp(sprintf('angle %g of %g', ia, na));
%
% first backproject at theta = 0 (copy/smear the contents of the projection
    % throughout the image)
    temp = zeros(nx); % Creating stencil for rotation
    for y_pos = 1:1:nx
        temp(y_pos, :) = sg1(:, ia)';
    end
    rot_angle = ((ia-1)*180)/na; % in degree
    %    
    % now rotate the projection
    % !!!Warning - imrotate expects angles in degrees!!!
    backproject_theta = imrotate(temp, rot_angle, 'bicubic','crop'); 
    % 
    lamin = lamin + backproject_theta;
  end
%
% Display Image
%
  figure(4)
  imagesc(x, y, lamin); colormap('gray'); axis('image')  
  title('Simple Backprojection Image')
  xlabel('Position')
  ylabel('Position')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   QUESTION C1 - FILTER PROJECTIONS 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
 disp('Question C, part 1: filter the sinogram')
%  if 1~=exist('sinogramfilt')
  % zero pad sinogram (if using Fourier methods)
  % Implementing the filter
  ramlak = abs(linspace(-1,1,nr).');
  ramlak = repmat(ramlak, [1 na]);
  
  sinogrampad = sg1;%padarray(sg1, nr/2, 0);
  sinogramfft = fft(sinogrampad, [], 1);
  sinogramfft = fftshift(sinogramfft, 1);
  sinogramfiltfreq = ifftshift(sinogramfft.*ramlak, 1);
  % 
  % filter the sinogram
  % !!!Warning - fftshift works differently for vectors and matrices!!!
  sinogramfilt = real(ifft(sinogramfiltfreq, nr, 1));
  %sinogramfilt((nr/2 + 1):(end-(nr/2)), :);
%
% Plot Filtered Sinogram
%
  figure(5)
  plot(r, sinogram(:,64)./max(sinogram(:,64)), '-',...
       r, sinogramfilt(:,64)./max(sinogramfilt(:,64)),':');
  title('Filtered sinogram')
  legend('Original Sinogram', 'Filtered SInogram')
  xlabel('Position');
  ylabel('Projection');
%
% Display Filtered Sinogram
% 
  figure(6)
  imagesc(r, angle/pi*180, sinogramfilt'); colormap('gray'); axis('image')  
  title('Filtered sinogram')
  xlabel('angle')
  ylabel('Position(rays)')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   QUESTION C2 - BACKPROJECT THE FILTERED SINOGRAMS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
 disp('Question C, part 2: Filtered backprojection')
%  if 1~=exist('bpf_recon')  % Just checking ... 
%  
  bpf_recon = zeros(nx,ny);
  for ia = 1:na
    disp(sprintf('angle %g of %g', ia, na))
%
    % first backproject at theta = 0 (copy the contents of the projection
    % through the image)
    temp = zeros(nx); % Creating stencil for rotation
    for y_pos = 1:1:nx
        temp(y_pos, :) = sinogramfilt(:, ia)';
    end
    rot_angle = ((ia-1)*180)/na; % in degree
    %    
    % now rotate the projection
    % !!!Warning - imrotate expects angles in degrees!!!
    backproject_theta = imrotate(temp, rot_angle, 'bicubic','crop'); 
    % 
    bpf_recon = bpf_recon + backproject_theta;
  end
%
% Display Reconstructed Image with Negative Values Set to Zero
%
  figure(7)
  imagesc(x, y, max(bpf_recon,0)); colormap('gray'); axis('image')  
  title('Filtered Backprojection Image')
  xlabel('Position')
  ylabel('Position')
% end
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %   QUESTION D - COMPARE RECONSTRUCTIONS, PLOT PROFILES
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
 disp('Question D: compare reconstructions')
%  Prepare the lines for comparison
%   y = -127:2:127;
index = find(y==-7);
% Prepare lines for comparison 
% Warning: you may want to check that your recons are not rotated
  line_phantom = phantom(:, index);
  line_lamin = lamin(index,:);
  line_bpf = bpf_recon(index, :);
%
% Plot Reconstructed Values (normalized) 
%
  figure(8)
  plot(y, line_phantom/max(line_phantom),'r--', ...
      y, line_bpf/max(line_bpf),'b-', ...
      y, line_lamin/max(line_lamin), 'k-');
  legend('Original', 'Filtered BP', 'Laminogram')
  title('Reconstruction Comparison')
  xlabel('Horizontal position (mm)');
  ylabel('Normalized attenuation');
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %   QUESTION E - RECONSTRUCTION USING SUBSAMPLED SINOGRAM
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% 
%  Downsample/subsample the original sinogram by removing 7 out of 8 projections   
    sinogram8 = sinogram(:,1:8:end); 
%  Using the result of Question C, reconstruct the image
    % Developing the filter
    [dim1 dim2] = size(sinogram8);
    ramlak = abs(linspace(-1,1,nr).');
    ramlak = repmat(ramlak, [1 dim2]);
  
    sinogram8fft = fft(sinogram8, [], 1);
    sinogram8fft = fftshift(sinogram8fft, 1);
    sinogramfilt8freq = ifftshift(sinogram8fft.*ramlak, 1);
    sinogramfilt8 = real(ifft(sinogramfilt8freq, nr, 1));
    bpf_recon8 = zeros(nx,ny);
    for ia = 1:na/8
        temp = zeros(nx); % Creating stencil for rotation
        for y_pos = 1:1:nx
            temp(y_pos, :) = sinogramfilt(:, ia)';
        end
        rot_angle = ((ia-1)*180)/(na/8); % in degree
        backproject_theta = imrotate(temp, rot_angle, 'bicubic','crop'); 
        bpf_recon8 = bpf_recon8 + backproject_theta;     
    end   
%
% Display Reconstructed Image with Negative Values Set to Zero
%
  figure(9)
  imagesc(x, y, max(bpf_recon8,0)); colormap('gray'); axis('image')  
  title('Subsampled Filtered Backprojection Image')
  xlabel('Position')
  ylabel('Position')
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %   QUESTION F - MYSTERY OBJECT
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  MISTERY SINOGRAM
%
load sinogram2; sinogram = sg2; 
%
%	We must redefine some of the parameters again
%
 dr=1;
 [nr, na] = size(sinogram);
 nx=nr; ny=nr; dx=dr;  % 
 angle = [0:(na-1)]'/na * pi;
 r = dr * ([1:nr]'-(nr+1)/2);
    disp(sprintf('number of rays = %g', nr))
    disp(sprintf('number of views = %g', na))
%
x = dx * ([1:nx]'-(nx+1)/2);
y = -dx * ([1:ny]'-(ny+1)/2);
xx = x(:,ones(1,ny));
yy = y(:,ones(1,ny))';
figure(10)
    imagesc(r, angle/pi*180, sg2')   % NOTE the transpose (') here and
    colormap('gray')                 % the fact that angle is displayed in degrees
    title('Sinogram: Unknown')
    xlabel('Position (i.e., Rays)')
    ylabel('Angle (i.e., Views)')
%
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %   QUESTION F2 - IMPLEMENT SIMPLE BACKPROJECTION, i.e. 
% %                PRODUCE LAMINOGRAMS
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% you may want to cut-and-paste your code from Question B here
   lamin = zeros(nx,ny);
   mid_pt = nx/2;
   for ia = 1:na
     disp(sprintf('angle %g of %g', ia, na));
     temp = zeros(nx); % Creating stencil for rotation
     for y_pos = 1:1:nx
         temp(y_pos, :) = sg2(:, ia)';
     end
     rot_angle = ((ia-1)*180)/na; % in degree
     backproject_theta = imrotate(temp, rot_angle, 'bicubic','crop'); 
     lamin = lamin + backproject_theta;
   end
%
% Display Image
%
  figure(12)
  imagesc(x, y, lamin); colormap('gray'); axis('image') 
  title('Simple Backprojection Image')
  xlabel('Position')
  ylabel('Position')
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %   QUESTION F2.1 - FILTER PROJECTIONS 
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% you may want to cut-and-paste your code from Question C1 here
   ramlak = abs(linspace(-1,1,nr).');
   ramlak = repmat(ramlak, [1 na]);
   sinogram2pad = sg2;
   sinogram2fft = fft(sinogram2pad, [], 1);
   sinogram2fft = fftshift(sinogram2fft, 1);
   sinogram2filtfreq = ifftshift(sinogram2fft.*ramlak, 1);
   sinogramfilt2 = real(ifft(sinogram2filtfreq, nr, 1));
  
%
% Plot Filtered Sinogram
%
  figure(13)
  plot(r, sinogram(:,64)./max(sinogram(:,64)), '-',...
       r, sinogramfilt2(:,64)./max(sinogramfilt2(:,64)),':');
  xlabel('Position');
  ylabel('Projection');
  
  figure(14)
  imagesc(r, angle/pi*180, sinogramfilt2'); colormap('gray'); axis('image')  
  title('Filtered sinogram')
  xlabel('Position(Rays)')
  ylabel('Angle')
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %   QUESTION F2.2 - BACKPROJECT THE FILTERED SINOGRAMS
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% you may want to cut-and-paste your code from Question C2 here
    bpf_recon2 = zeros(nx,ny);
    for ia = 1:na
        disp(sprintf('angle %g of %g', ia, na))
        temp = zeros(nx); % Creating stencil for rotation
        for y_pos = 1:1:nx
            temp(y_pos, :) = sinogramfilt2(:, ia)';
        end
        rot_angle = ((ia-1)*180)/na; % in degree
        backproject_theta = imrotate(temp, rot_angle, 'bicubic','crop'); 
        bpf_recon2 = bpf_recon2 + backproject_theta;
    end
%
% Display Reconstructed Image with Negative Values Set to Zero
%
  figure(15)
  imagesc(x, y, max(bpf_recon2,0)); colormap('gray'); axis('image')
  title('Filtered Backprojection Image')
  xlabel('Position')
  ylabel('Position')
% %
