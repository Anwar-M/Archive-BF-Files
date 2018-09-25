function A = ConventionalBeamforming(x, CSM, mic_positions, frequencies, speed_of_sound, M)

    n_mics = length(mic_positions);
    n_freqs = numel(frequencies);

    r2 = mic_positions - x'*ones(1, n_mics);
    beta_sq = 1 - dot(M, M);
    M_array = M'*ones(1, n_mics);

    A = zeros(n_freqs, 1);

    for I = 1:n_freqs

        R_adapted = sqrt( dot(M_array,r2,1).^2 + beta_sq.*dot(r2,r2,1) );
        xi_timedelay = (-dot(M_array,r2,1) + R_adapted)/(beta_sq*speed_of_sound); 
        g = (-exp(-2.*pi.*1i.*frequencies(I).*xi_timedelay)) ./ ( 4*pi*R_adapted );

        A(I) = conj(g)*(squeeze(CSM(I,:,:)).')*g.'/(dot(g,g)^2);

    end

end