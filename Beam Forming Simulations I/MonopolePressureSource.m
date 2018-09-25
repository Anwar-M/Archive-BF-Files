% Equation (19) of Acoustic beamforming for the ranking of aircraft noise
% by Sijtsma NLR. Source is sent at t = 0. Signal is expected as an vector
% as a function of time (not efficient).

function [t, p] = MonopolePressureSource(sigma, fs, x, xi, M, c)
    beta_sq = 1 - dot(M, M);
    D = sqrt( dot(M, x - xi)^2 + beta_sq*dot(x-xi, x-xi) );
    dt_e = ( -dot(M, x - xi) + D ) / (c * beta_sq);
    
    disp(dt_e);
    
    t_shift_sample = round(dt_e * fs);
    t = dt_e:(1/fs):dt_e + (length(sigma)-1)/fs;
    
%     if t_shift_sample > length(sigma) 
%         error('Shift in time larger than signal.'); 
%     else
%         sigma = sigma(t_shift_sample:end);
%     end
    
    p = + sigma / (4 * pi * D);
    
    beta_sq = 1 - dot(M_wind, M_wind);
    ti = -4e-3:1/fs:4e-3-1/fs;
    p = zeros(n_mics, 8e-3*fs);
    ti = ti + 4e-3;
    yi = 2*sqrt(2)*cos(2*pi*f_source*ti);
    
    for I = 1:size(xi,1)
        for J = 1:n_mics
            
            r = sqrt( dot(M_wind, x(J,:) - xi(I,:))^2 + beta_sq*dot(x(J,:)-xi(I,:), x(J,:)-xi(I,:)) );
            ph = ( -dot(M_wind, x(J,:) - xi(I,:)) + r ) / (c_sound * beta_sq);
            
%             r = sqrt(sum((x(J, :) - xi(I, :)).^2));
%             ph = r/c_sound;
            delayed_first_sample = round(ph*fs);
            p(J, :) = p(J, :) + 2*sqrt(2)*cos(2*pi*f_source*(ti-ph))/(4 * pi * r);
        end
    end    
end