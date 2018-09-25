function A = ConventionalBeamforming(recorded_signal, sample_frequency, ...
    mic_position, scan_grid, imaged_frequencies, speed_of_sound, M)

n_mics = length(mic_position);
n_samples = size(recorded_signal, (size(recorded_signal) ~= n_mics)*[1; 2]);
n_freqs = numel(imaged_frequencies);
n_scan_x = size(scan_grid, 1); 
n_scan_y = size(scan_grid, 2);

Numerator = zeros(n_scan_x, n_scan_y);
Denominator = Numerator;

A = zeros(n_scan_x, n_scan_y, n_freqs);
% M(1) = M(1)+0.075;
beta_sq = 1 - dot(M, M);
M_array(:,:,1) = M(1)*ones(n_scan_x, n_scan_y);
M_array(:,:,2) = M(2)*ones(n_scan_x, n_scan_y);
M_array(:,:,3) = M(3)*ones(n_scan_x, n_scan_y);
for K = 1:n_freqs

    for I = 1:n_mics
%         r = sqrt( ...
%                   (mic_position(I,1) - scan_grid(:,:,1)).^2 + ...
%                   (mic_position(I,2) - scan_grid(:,:,2)).^2 + ...
%                   (mic_position(I,3) - scan_grid(:,:,3)).^2 ...
%                 );
% 
%         xi_timedelay= r/speed_of_sound;
%         g = (-exp(-2.*pi.*1i.*imaged_frequencies(K).*xi_timedelay)) ./ (4*pi*r);
        
        r2(:,:,1) = mic_position(I,1) - scan_grid(:,:,1);
        r2(:,:,2) = mic_position(I,2) - scan_grid(:,:,2);
        r2(:,:,3) = mic_position(I,3) - scan_grid(:,:,3);
        
        R_adapted = sqrt( dot(M_array,r2,3) + beta_sq.*dot(r2,r2,3) );
        xi_timedelay = (-dot(M_array,r2,3)+R_adapted)/(beta_sq*speed_of_sound); 
        g = (-exp(-2.*pi.*1i.*imaged_frequencies(K).*xi_timedelay)) ./ ( 4*pi*R_adapted );

        P = fft(recorded_signal(I, :));
        f = (0:n_samples/2-1)*sample_frequency/n_samples;
        [~, imaged_frequency_index] = min(abs(f-imaged_frequencies(K)));
        pf = 2*P(imaged_frequency_index)/n_samples;

        Numerator = Numerator + conj(g)*pf;
        Denominator = Denominator + dot(g,g,3);
    end
    a = Numerator./Denominator;
    A(:,:,K) = dot(Numerator./Denominator, Numerator./Denominator, 3)/2;
    Numerator(:,:) = 0;
    Denominator(:,:) = 0;
end

A = sqrt(sum(A.^2, 3)/2);

end