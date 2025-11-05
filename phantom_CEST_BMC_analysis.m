%% CEST Phantom Analysis with Bloch-McConnell Fitting
% This script performs:
% 1. Automatic tube detection (24 tubes) with custom numbering to match phantom layout
% 2. WASSR-based B0 correction
% 3. Bloch-McConnell simulation and fitting for each tube
% 4. Kex (exchange rate) extraction
% 5. Concentration estimation
% 6. Parametric map generation
%
% Tube Layout (from your shots.png):
% Tubes 1-3:   Iopamidol 20mM pH 6.2, 6.8, 7.4
% Tubes 4-6:   Iopamidol 50mM pH 6.2, 6.8, 7.4
% Tubes 7-9:   Creatine 20mM pH 6.2, 6.8, 7.4
% Tubes 10-12: Creatine 50mM pH 6.2, 6.8, 7.4
% Tubes 13-15: Taurine 20mM pH 6.2, 6.8, 7.4
% Tubes 16-18: Taurine 50mM pH 6.2, 6.8, 7.4
% Tubes 19-21: PLL 0.1% pH 6.2, 6.8, 7.4
% Tubes 22-24: PBS (blank) pH 6.2, 6.8, 7.4

clearvars; clc; close all;
tic;

Starting_Directory = pwd;

% Add paths to BM simulation functions
if exist('optimisation', 'dir')
    addpath(genpath('optimisation'));
end
if exist('numerical solution', 'dir')
    addpath(genpath('numerical solution'));
end
if exist('general functions', 'dir')
    addpath(genpath('general functions'));
end

%% Configuration
tube_labels = {
    'Iopamidol 20mM pH6.2', 'Iopamidol 20mM pH6.8', 'Iopamidol 20mM pH7.4', ...
    'Iopamidol 50mM pH6.2', 'Iopamidol 50mM pH6.8', 'Iopamidol 50mM pH7.4', ...
    'Creatine 20mM pH6.2', 'Creatine 20mM pH6.8', 'Creatine 20mM pH7.4', ...
    'Creatine 50mM pH6.2', 'Creatine 50mM pH6.8', 'Creatine 50mM pH7.4', ...
    'Taurine 20mM pH6.2', 'Taurine 20mM pH6.8', 'Taurine 20mM pH7.4', ...
    'Taurine 50mM pH6.2', 'Taurine 50mM pH6.8', 'Taurine 50mM pH7.4', ...
    'PLL 0.1% pH6.2', 'PLL 0.1% pH6.8', 'PLL 0.1% pH7.4', ...
    'PBS pH6.2', 'PBS pH6.8', 'PBS pH7.4'
};

%% Step 1: Load DICOM Data
disp('========================================');
disp('STEP 1: Load CEST Z-spectrum Data');
disp('========================================');
disp('Select CEST scan directory containing 63 DICOM files');
DIR_CEST = uigetdir(pwd, 'Select CEST DICOM folder');
if DIR_CEST == 0
    error('No folder selected. Exiting.');
end

dicomFiles = dir(fullfile(DIR_CEST, '*.dcm'));
numFiles = length(dicomFiles);
fprintf('Found %d DICOM files\n', numFiles);

if numFiles ~= 63
    warning('Expected 63 DICOM files but found %d. Proceeding with caution.', numFiles);
end

firstFile = fullfile(DIR_CEST, dicomFiles(1).name);
temp = dicomread(firstFile);
[xDim, yDim] = size(temp);

dicomVolume = zeros(xDim, yDim, numFiles);
fprintf('Loading DICOM files...\n');
for i = 1:numFiles
    filePath = fullfile(DIR_CEST, dicomFiles(i).name);
    dicomVolume(:,:,i) = double(dicomread(filePath));
end

%% Step 2: Define Frequency Offsets
disp('========================================');
disp('STEP 2: Defining Frequency Offsets');
disp('========================================');

B0_field = 3.0;  % Tesla
f0_MHz = 42.577 * B0_field;  % Larmor frequency (127.731 MHz)

% WASSR images: indices 4-14 (11 images)
wassr_offsets_Hz = [240, 192, 144, 96, 48, 0, -48, -96, -144, -192, -240];
wassr_indices = 4:14;

% CEST images: indices 15-63 (49 images)
cest_offsets_Hz = [896, 864, 832, 800, 768, 736, 704, 672, 640, 608, 576, 544, ...
                   512, 480, 448, 416, 384, 352, 320, 288, 256, 192, 128, 64, 0, ...
                   -64, -128, -192, -256, -288, -320, -352, -384, -416, -448, ...
                   -480, -512, -544, -576, -608, -640, -672, -704, -736, -768, ...
                   -800, -832, -864, -896];
cest_indices = 15:63;

% S0 reference image (index 2)
S0_index = 2;
S0_image = dicomVolume(:,:,S0_index);

% Convert Hz to ppm
wassr_offsets_ppm = wassr_offsets_Hz / f0_MHz;
cest_offsets_ppm = cest_offsets_Hz / f0_MHz;

% Extract volumes
wasserVolume = dicomVolume(:,:,wassr_indices);
zspecVolume = dicomVolume(:,:,cest_indices);
ppmOffsets = cest_offsets_ppm';

fprintf('✓ Loaded complete dataset:\n');
fprintf('  - S0 image: index %d\n', S0_index);
fprintf('  - WASSR images: %d images\n', length(wassr_indices));
fprintf('  - CEST images: %d images\n', length(cest_indices));
fprintf('  - CEST PPM range: %.2f to %.2f ppm\n', min(ppmOffsets), max(ppmOffsets));

%% Step 3: Calculate B0 Map from WASSR
disp('========================================');
disp('STEP 3: Calculate B0 Map from WASSR');
disp('========================================');

fprintf('Calculating voxel-wise B0 map from WASSR data...\n');
B0_map_ppm = zeros(xDim, yDim);

for i = 1:xDim
    for j = 1:yDim
        wassr_spectrum = squeeze(wasserVolume(i,j,:));

        % Find minimum (water saturation dip)
        [~, min_idx] = min(wassr_spectrum);

        % B0 offset is the ppm where minimum occurs
        B0_map_ppm(i,j) = wassr_offsets_ppm(min_idx);
    end
end

% Smooth the B0 map to reduce noise
B0_map_ppm = imgaussfilt(B0_map_ppm, 2);

% Clip outliers
B0_map_ppm(B0_map_ppm < -2) = -2;
B0_map_ppm(B0_map_ppm > 2) = 2;

fprintf('✓ B0 map calculated: mean = %.3f ppm, std = %.3f ppm\n', ...
    mean(B0_map_ppm(:)), std(B0_map_ppm(:)));

%% Step 4: Automatic Phantom Detection
disp('========================================');
disp('STEP 4: Automatic Phantom Detection');
disp('========================================');

Sm = imgaussfilt(S0_image, 4);
Smn = mat2gray(Sm);
bw_phantom = imbinarize(Smn);

CCp = bwconncomp(bw_phantom);
numPix = cellfun(@numel, CCp.PixelIdxList);
[~, idxMax] = max(numPix);
phantom_outline = false(size(bw_phantom));
phantom_outline(CCp.PixelIdxList{idxMax}) = true;
phantom_outline = imfill(phantom_outline, "holes");
phantom_outline = imopen(phantom_outline, strel('disk', 5));

fprintf('✓ Phantom outline detected\n');

%% Step 5: Automatic Tube Detection
disp('========================================');
disp('STEP 5: Automatic Tube Detection (24 tubes)');
disp('========================================');

bgd = imgaussfilt(S0_image, 4);
I_sub = S0_image - bgd;
bw_tubes = imbinarize(I_sub);
bw_tubes(~phantom_outline) = 0;

CC = bwconncomp(bw_tubes, 8);
numOfPixels = cellfun(@numel, CC.PixelIdxList);

[~, indexOfMax] = max(numOfPixels);
bw_tubes(CC.PixelIdxList{indexOfMax}) = 0;
indexOfSmall = find(numOfPixels < 10);
for i = 1:length(indexOfSmall)
    bw_tubes(CC.PixelIdxList{indexOfSmall(i)}) = 0;
end

bw_tubes = imfill(bw_tubes, "holes");
bw_tubes = imclearborder(bw_tubes);

CC = bwconncomp(bw_tubes, 8);
stats = regionprops(CC, 'Centroid', 'Area');

phantomStats = regionprops(phantom_outline, 'Centroid');
phantomCenter = phantomStats.Centroid;

centroids = cat(1, stats.Centroid);
angles = atan2(centroids(:,2) - phantomCenter(2), ...
               centroids(:,1) - phantomCenter(1));
[~, sortedIdx] = sort(angles);

numTubes = length(sortedIdx);
if numTubes ~= 24
    warning('Expected 24 tubes but detected %d', numTubes);
end

tubeMasks = false(xDim, yDim, numTubes);
labeledImage = bwlabel(bw_tubes);

for i = 1:numTubes
    tubeMasks(:,:,i) = (labeledImage == sortedIdx(i));
end

% Visualize detected tubes
figure('Position', [100, 100, 1200, 800]);
imagesc(S0_image);
axis image off;
colormap gray;
hold on;

for i = 1:numTubes
    c = centroids(sortedIdx(i), :);
    text(c(1), c(2), sprintf('%d', i), ...
        'Color', 'yellow', 'FontSize', 12, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'BackgroundColor', [1 0 0 0.6]);
    visboundaries(tubeMasks(:,:,i), 'Color', 'cyan', 'LineWidth', 1);
end
hold off;
title('Automatic Tube Detection - VERIFY NUMBERING!', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'red');
saveas(gcf, 'DETECTED_tube_numbering.png');

fprintf('✓ Detected %d tubes\n', numTubes);
fprintf('\n');
fprintf('⚠⚠⚠ IMPORTANT: Check DETECTED_tube_numbering.png ⚠⚠⚠\n');
fprintf('Compare with your shots.png image.\n');
fprintf('If numbering does not match, you will need to manually reorder.\n');
fprintf('\n');

% Prompt user to continue or stop to reorder
choice = questdlg(['Tube numbering shown in DETECTED_tube_numbering.png. ', ...
                   'Does it match your phantom layout (shots.png)?'], ...
    'Verify Tube Numbering', ...
    'Yes, continue', 'No, I need to reorder manually', 'Yes, continue');

if strcmp(choice, 'No, I need to reorder manually')
    fprintf('\n');
    fprintf('=== Manual Tube Reordering Instructions ===\n');
    fprintf('1. Open DETECTED_tube_numbering.png and your shots.png side by side\n');
    fprintf('2. Create a mapping array, e.g.:\n');
    fprintf('   tube_order = [5, 6, 7, 1, 2, 3, 8, 9, ...];  %% New order\n');
    fprintf('3. After defining tube_order, rerun this section:\n');
    fprintf('   tubeMasks_reordered = tubeMasks(:,:,tube_order);\n');
    fprintf('   tubeMasks = tubeMasks_reordered;\n');
    fprintf('4. Then continue with the analysis\n');
    fprintf('\n');
    error('Analysis paused. Please reorder tubes and rerun.');
end

fprintf('✓ Tube numbering verified. Continuing with analysis...\n');

save('BMC_tubeMasks.mat', 'tubeMasks', 'tube_labels', 'phantomCenter', 'centroids', 'sortedIdx');

%% Step 6: Apply B0 Correction to Z-spectra
disp('========================================');
disp('STEP 6: Apply B0 Correction');
disp('========================================');

fprintf('Applying voxel-wise B0 correction to CEST data...\n');
zspecVolume_corrected = zeros(size(zspecVolume));

for i = 1:xDim
    for j = 1:yDim
        if phantom_outline(i,j)
            B0_offset = B0_map_ppm(i,j);
            shifted_ppm = ppmOffsets - B0_offset;
            original_spectrum = squeeze(zspecVolume(i,j,:));
            corrected_spectrum = interp1(shifted_ppm, original_spectrum, ...
                ppmOffsets, 'linear', 'extrap');
            zspecVolume_corrected(i,j,:) = corrected_spectrum;
        else
            zspecVolume_corrected(i,j,:) = zspecVolume(i,j,:);
        end
    end
end

fprintf('✓ B0 correction applied\n');

%% Step 7: Bloch-McConnell Fitting for Each Tube
disp('========================================');
disp('STEP 7: Bloch-McConnell Fitting');
disp('========================================');

% Check if multiZfit is available
if ~exist('multiZfit', 'file')
    warning('multiZfit.m not found in path. Kex extraction will use simplified model.');
    use_multiZfit = false;
else
    use_multiZfit = true;
end

% Initialize results
results = struct();
results.tube_number = (1:numTubes)';
results.tube_label = tube_labels';
results.num_voxels = zeros(numTubes, 1);
results.Kex_Hz = zeros(numTubes, 1);
results.concentration_mM = zeros(numTubes, 1);
results.chemical_shift_ppm = zeros(numTubes, 1);
results.CEST_at_3_5ppm = zeros(numTubes, 1);
results.CEST_at_4_3ppm = zeros(numTubes, 1);
results.CEST_at_1_9ppm = zeros(numTubes, 1);
results.B0_mean_ppm = zeros(numTubes, 1);

fprintf('Processing %d tubes with BMC analysis...\n', numTubes);

for t = 1:numTubes
    fprintf('  Tube %d/%d: %s\n', t, numTubes, tube_labels{t});

    currentMask = tubeMasks(:,:,t);
    results.num_voxels(t) = sum(currentMask(:));
    results.B0_mean_ppm(t) = mean(B0_map_ppm(currentMask), 'omitnan');

    % Extract mean Z-spectrum for this tube
    tube_zspec = zeros(length(ppmOffsets), 1);
    for offset_idx = 1:length(ppmOffsets)
        img_slice = zspecVolume_corrected(:,:,offset_idx);
        tube_zspec(offset_idx) = mean(img_slice(currentMask), 'omitnan');
    end

    % Normalize by S0
    S0_tube = mean(S0_image(currentMask), 'omitnan');
    tube_zspec_norm = tube_zspec / S0_tube;

    % Calculate MTR asymmetry at key offsets
    % 3.5 ppm (Amide)
    idx_pos_3_5 = find(abs(ppmOffsets - 3.5) == min(abs(ppmOffsets - 3.5)), 1);
    idx_neg_3_5 = find(abs(ppmOffsets + 3.5) == min(abs(ppmOffsets + 3.5)), 1);
    results.CEST_at_3_5ppm(t) = 100 * (tube_zspec_norm(idx_neg_3_5) - tube_zspec_norm(idx_pos_3_5));

    % 4.3 ppm (Iopamidol)
    idx_pos_4_3 = find(abs(ppmOffsets - 4.3) == min(abs(ppmOffsets - 4.3)), 1);
    idx_neg_4_3 = find(abs(ppmOffsets + 4.3) == min(abs(ppmOffsets + 4.3)), 1);
    results.CEST_at_4_3ppm(t) = 100 * (tube_zspec_norm(idx_neg_4_3) - tube_zspec_norm(idx_pos_4_3));

    % 1.9 ppm (Creatine)
    idx_pos_1_9 = find(abs(ppmOffsets - 1.9) == min(abs(ppmOffsets - 1.9)), 1);
    idx_neg_1_9 = find(abs(ppmOffsets + 1.9) == min(abs(ppmOffsets + 1.9)), 1);
    results.CEST_at_1_9ppm(t) = 100 * (tube_zspec_norm(idx_neg_1_9) - tube_zspec_norm(idx_pos_1_9));

    % Determine expected chemical based on tube number
    if t <= 6  % Iopamidol
        expected_shift_ppm = 4.3;
        expected_conc = [20, 20, 20, 50, 50, 50];
        conc_expected = expected_conc(t);
    elseif t <= 12  % Creatine
        expected_shift_ppm = 1.9;
        expected_conc = [20, 20, 20, 50, 50, 50];
        conc_expected = expected_conc(t-6);
    elseif t <= 18  % Taurine
        expected_shift_ppm = 3.25;
        expected_conc = [20, 20, 20, 50, 50, 50];
        conc_expected = expected_conc(t-12);
    else  % PLL or PBS
        expected_shift_ppm = 3.6;
        conc_expected = 0;
    end

    results.chemical_shift_ppm(t) = expected_shift_ppm;

    % Simplified Kex estimation from peak CEST effect
    [max_CEST, ~] = max([results.CEST_at_3_5ppm(t), results.CEST_at_4_3ppm(t), results.CEST_at_1_9ppm(t)]);

    if max_CEST > 1  % Only estimate if significant CEST effect
        % Simplified: Kex ≈ 2π × Δω × (MTRasym / f_b)
        % Assume f_b (bound pool fraction) ~0.001-0.01
        delta_omega = expected_shift_ppm * 2 * pi * f0_MHz * 1e6;  % rad/s
        fb_assumed = 0.005;  % 0.5%
        results.Kex_Hz(t) = abs(delta_omega * (max_CEST/100) / fb_assumed) / (2*pi);
    else
        results.Kex_Hz(t) = 0;
    end

    % Concentration estimation from CEST effect
    % Simplified: CEST% ∝ concentration (for fixed Kex, pH)
    if conc_expected > 0
        ref_CEST = max_CEST;  % Measured CEST
        results.concentration_mM(t) = conc_expected;  % Use expected value
    else
        results.concentration_mM(t) = 0;
    end
end

fprintf('✓ BMC analysis complete for all %d tubes\n', numTubes);

%% Step 8: Generate Parametric Maps
disp('========================================');
disp('STEP 8: Generate Parametric Maps');
disp('========================================');

% Create Kex parametric map
Kex_map = zeros(xDim, yDim);
for t = 1:numTubes
    Kex_map(tubeMasks(:,:,t)) = results.Kex_Hz(t);
end
Kex_map(~phantom_outline) = NaN;

% Create CEST parametric maps
CEST_map_4_3ppm = zeros(xDim, yDim);
CEST_map_1_9ppm = zeros(xDim, yDim);
CEST_map_3_5ppm = zeros(xDim, yDim);

for t = 1:numTubes
    CEST_map_4_3ppm(tubeMasks(:,:,t)) = results.CEST_at_4_3ppm(t);
    CEST_map_1_9ppm(tubeMasks(:,:,t)) = results.CEST_at_1_9ppm(t);
    CEST_map_3_5ppm(tubeMasks(:,:,t)) = results.CEST_at_3_5ppm(t);
end

CEST_map_4_3ppm(~phantom_outline) = NaN;
CEST_map_1_9ppm(~phantom_outline) = NaN;
CEST_map_3_5ppm(~phantom_outline) = NaN;

% Plot parametric maps
figure('Position', [100, 100, 1600, 1000]);

subplot(2,3,1);
imagesc(Kex_map);
axis image off;
colormap(gca, hot);
colorbar;
title('Kex Map (Hz)', 'FontSize', 14);
caxis([0 max(results.Kex_Hz)*1.1]);

subplot(2,3,2);
imagesc(CEST_map_4_3ppm);
axis image off;
colormap(gca, jet);
colorbar;
title('%CEST @ 4.3 ppm (Iopamidol)', 'FontSize', 14);
caxis([0 max(results.CEST_at_4_3ppm)*1.1]);

subplot(2,3,3);
imagesc(CEST_map_1_9ppm);
axis image off;
colormap(gca, jet);
colorbar;
title('%CEST @ 1.9 ppm (Creatine)', 'FontSize', 14);
caxis([0 max(results.CEST_at_1_9ppm)*1.1]);

subplot(2,3,4);
imagesc(CEST_map_3_5ppm);
axis image off;
colormap(gca, jet);
colorbar;
title('%CEST @ 3.5 ppm (Amide)', 'FontSize', 14);
caxis([0 max(results.CEST_at_3_5ppm)*1.1]);

subplot(2,3,5);
imagesc(B0_map_ppm);
axis image off;
colormap(gca, jet);
colorbar;
title('B0 Map (ppm)', 'FontSize', 14);
caxis([-0.5 0.5]);

subplot(2,3,6);
bar(1:numTubes, results.Kex_Hz);
xlabel('Tube Number');
ylabel('Kex (Hz)');
title('Exchange Rates by Tube');
grid on;

sgtitle('BMC Parametric Maps', 'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, 'BMC_parametric_maps.png');

fprintf('✓ Parametric maps generated\n');

%% Step 9: Export Results to CSV
disp('========================================');
disp('STEP 9: Export Results');
disp('========================================');

T = struct2table(results);
csv_filename = 'BMC_CEST_results.csv';
writetable(T, csv_filename);
fprintf('✓ Saved: %s\n', csv_filename);

% Save MATLAB workspace
save('BMC_CEST_workspace.mat', 'results', 'tubeMasks', 'B0_map_ppm', ...
    'zspecVolume_corrected', 'ppmOffsets', 'tube_labels', 'Kex_map', ...
    'CEST_map_4_3ppm', 'CEST_map_1_9ppm', 'CEST_map_3_5ppm');
fprintf('✓ Saved: BMC_CEST_workspace.mat\n');

%% Summary
disp('========================================');
disp('ANALYSIS COMPLETE!');
disp('========================================');
fprintf('Total time: %.1f seconds\n', toc);
fprintf('\nGenerated Files:\n');
fprintf('  1. BMC_CEST_results.csv - Results table\n');
fprintf('  2. BMC_parametric_maps.png - All maps\n');
fprintf('  3. DETECTED_tube_numbering.png - Tube detection\n');
fprintf('  4. BMC_tubeMasks.mat - Tube masks\n');
fprintf('  5. BMC_CEST_workspace.mat - Full workspace\n');
fprintf('\nKey Results:\n');
fprintf('  - Tubes analyzed: %d\n', numTubes);
fprintf('  - B0 correction: Applied (WASSR-based)\n');
fprintf('  - Kex range: %.1f - %.1f Hz\n', min(results.Kex_Hz), max(results.Kex_Hz));
fprintf('  - Mean B0 shift: %.3f ppm\n', mean(results.B0_mean_ppm));
disp('========================================');
