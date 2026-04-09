function helfrich_fitting_model()
% =========================================================================
% Script Name: helfrich_fitting_model
% Author: Hongyan Yuan (SUSTech) hongyan6@outlook.com         
%         Xin Zhang (SUSTech) 12432359@mail.sustech.edu.cn 
% Description: This script calculates the bending rigidity (kappa) and 
%              surface tension (sigma) of a lipid membrane from LAMMPS 
%              trajectories using the Helfrich fluctuation model.
% =========================================================================
%                               INPUT PARAMETERS
% =========================================================================
% Lx, Ly    - System dimensions in X and Y directions 
% n         - Grid resolution for surface interpolation
% kbT       - Thermal energy (consistent with simulation units)
% filename  - Path to the LAMMPS trajectory file (.lammpstrj)
% begini    - Starting timestep for analysis
% endi      - Ending timestep for analysis
% interval  - Timestep frequency for frame sampling
% =========================================================================
    Lx = 98;
    Ly = 101;
    n = 80;
    kbT = 0.23; 
    prefix = 'dump_4_onlymembrane';
    filename = [prefix '.lammpstrj'];
    begini = 100000;
    endi   = 200000;
    interval = 5000;

    xi = linspace(0, Lx, n+1);
    yi = linspace(0, Ly, n+1);
    xi = xi(1:n) + (Lx/n)/2;
    yi = yi(1:n) + (Ly/n)/2;
    [XI, YI] = meshgrid(xi, yi);
    
    hq2_sum = zeros(n, n); 
    fid = fopen(filename, 'r'); 
    current_frame_count = 0; 
    
    while ~feof(fid)
        tline = fgetl(fid);
        if ~ischar(tline), break; end
        if strncmp(tline, 'ITEM: TIMESTEP', 14) 
            timestep = str2double(fgetl(fid)); 
            fgetl(fid); 
            natoms_current = str2double(fgetl(fid));
            
            if timestep >= begini && timestep <= endi && mod(timestep - begini, interval) == 0
                fgetl(fid); % BOX BOUNDS
                box_x = sscanf(fgetl(fid), '%f %f');
                box_y = sscanf(fgetl(fid), '%f %f');
                box_z = sscanf(fgetl(fid), '%f %f');
                fgetl(fid); % ATOMS
                
                data_cell = textscan(fid, '%f %f %f %f %f', natoms_current);
                xyz = [data_cell{3}, data_cell{4}, data_cell{5}];
                
                xyz(:,1) = xyz(:,1) * (box_x(2) - box_x(1)) + box_x(1);
                xyz(:,2) = xyz(:,2) * (box_y(2) - box_y(1)) + box_y(1);
                xyz(:,3) = xyz(:,3) * (box_z(2) - box_z(1)) + box_z(1);
                
                ZI = griddata(xyz(:,1), xyz(:,2), xyz(:,3), XI, YI);
                ZI(isnan(ZI)) = nanmean(ZI(:)); 
                ZI = ZI - mean(ZI(:)); 
                [PX, PY] = meshgrid(xi, yi);
                A_mat = [PX(:), PY(:), ones(numel(PX), 1)];
                coeffs = A_mat \ ZI(:); 
                ZI = ZI - (coeffs(1)*PX + coeffs(2)*PY + coeffs(3)); 
                
                if current_frame_count == 0
                    figure('Color', 'w', 'Name', 'Real Space Membrane Surface');
                    surf(XI, YI, ZI, 'EdgeColor', 'none');
                    colormap jet; 
                    colorbar;
                    title(sprintf('Membrane Surface (Tilt Corrected) at Timestep %d', timestep));
                    xlabel('X'); ylabel('Y'); zlabel('Z (Height Fluctuation)');
                    
                    view(-30, 45); 
                    camlight; lighting gouraud; 
                    drawnow;
                end
                dx = Lx/n; dy = Ly/n;
                hq_t = fft2(ZI) * (dx * dy) / (Lx * Ly); 
                
                hq2_sum = hq2_sum + abs(hq_t).^2;
                current_frame_count = current_frame_count + 1;
                disp(['Timestep: ', num2str(timestep), ' OK']);
            else
                for sk=1:5, fgetl(fid); end
                textscan(fid, '%*s', natoms_current, 'Delimiter', '\n');
            end
        end
    end
    fclose(fid);
    
    hq2_avg = hq2_sum ./ current_frame_count;

    m = floor(n/2);
    q_list = []; hq2_list = [];
    
    for i = 0:m-1
        for j = 0:m-1
            if i==0 && j==0, continue; end
            qi = i * (2*pi/Lx);
            qj = j * (2*pi/Ly);
            q_val = sqrt(qi^2 + qj^2);
            q_list = [q_list; q_val];
            hq2_list = [hq2_list; hq2_avg(i+1, j+1)];
        end
    end


    m = floor(n/2);
    [qx, qy] = meshgrid((0:m-1)*(2*pi/Lx), (0:m-1)*(2*pi/Ly));
    q_map = sqrt(qx.^2 + qy.^2);

    S_map = (Lx * Ly) * hq2_avg(1:m, 1:m); 

    q_flat = q_map(:);
    S_flat = S_map(:);
    
    valid_idx = q_flat > 0;
    q_flat = q_flat(valid_idx);
    S_flat = S_flat(valid_idx);

    dq = min(2*pi/Lx, 2*pi/Ly); 
    max_q = max(q_flat);
    
    edges = 0 : dq : (max_q + dq);  
    
    bin_indices = discretize(q_flat, edges);
    
    num_bins = max(bin_indices);
    q_binned = zeros(num_bins, 1);
    S_binned = zeros(num_bins, 1);

    for k = 1:num_bins
        mask = (bin_indices == k); 
        if any(mask)
            q_binned(k) = mean(q_flat(mask));
            S_binned(k) = mean(S_flat(mask)); 
        end
    end

    valid_bins = S_binned > 0;

    q = q_binned(valid_bins);
    S_q = S_binned(valid_bins);

    logq = log(q); logS = log(S_q);
    
    %% %% Please enter the cutoff wavenumber and the number of high-frequency points to be omitted
    cutoff = 3;
    omitted = 40;
    idx_sigma = 1:cutoff;              
    idx_kappa = cutoff:length(q)-omitted;     

    b1 = mean(logS(idx_sigma) + 2*logq(idx_sigma));
    sigma = kbT / exp(b1);

    b2 = mean(logS(idx_kappa) + 4*logq(idx_kappa));
    kappa = kbT / exp(b2);

    fprintf('\n====================================\n');
    fprintf(' sigma = %.2f kbT (%.4e)\n', sigma/kbT, sigma);
    fprintf(' kappa = %.2f kbT (%.4e)\n', kappa/kbT, kappa);
    fprintf('====================================\n');

    figure('Color','w','Name','Membrane Fluctuation Spectrum');
    loglog(q, S_q, 'bo', 'MarkerSize', 4, 'DisplayName', 'MD Data'); hold on;
    q_low = q(idx_sigma);
    plot(q_low, exp(-2*log(q_low) + b1), 'r-', 'LineWidth', 2, 'DisplayName', '\sim q^{-2} (Sigma)');
    q_high = q(idx_kappa);
    plot(q_high, exp(-4*log(q_high) + b2), 'g-', 'LineWidth', 2, 'DisplayName', '\sim q^{-4} (Kappa)');
    
    xlabel('Wave vector q'); ylabel('S(q) = A \langle|h_q|^2\rangle');
    grid on; legend('Location', 'southwest');
    title('Helfrich Model Fitting');
end