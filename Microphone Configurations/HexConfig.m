function x = HexConfig(d)

    x = zeros(32,2);

    dphi = pi/3;

    k = 1;
    for j = 1:3
        ds = d*j;
        for i = 1:round(2*pi/dphi)
            x(k+1,1) = ds*cos(i*dphi);
            x(k+1,2) = ds*sin(i*dphi);
            k = k + 1;
        end
    end
    for i = 1:6
        x(k+1,1) = sqrt(3)*d*cos((i-1)*dphi+dphi/2);
        x(k+1,2) = sqrt(3)*d*sin((i-1)*dphi+dphi/2);
        k = k + 1;
    end
    x(26, :) = [x(2,1) x(14,2)];
    x(27,:) = [x(10,1) x(9,2)];
    x(28,:) = [(x(16,1)+x(10,1))/2 x(23,2)];
    x(29,:) = [x(5,1) x(17,2)];
    x(30,:) = [x(13,1) x(12,2)];
    x(31,:) = [(x(13,1)+x(19,1))/2 x(20,2)];
    x(32,1) = x(2,1); x(32,2) = x(18,2);
    
end