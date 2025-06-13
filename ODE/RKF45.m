function sol = RKF45(func, t0, tf, X0, dtInit, epsil)
    X = [X0];
    novars = length(X);
    t = [t0];
    i = 1;
    dt = [dtInit];
    while t(i)<tf;
        K1 = dt(i)*func(X(i,:), t(i));
        K2 = dt(i)*func(X(i,:) + 1/4*K1, t(i)+1/4*dt(i));
        K3 = dt(i)*func(X(i,:) + 3/32*K1+9/32*K2, t(i)+3/8*dt(i));
        K4 = dt(i)*func(X(i,:) + 1932/2197*K1 - 7200/2197*K2 + 7296/2197*K3, t(i) + 12/13*dt(i));
        K5 = dt(i)*func(X(i,:) + 439/216*K1 - 8*K2 + 3680/513*K3 - 845/4104*K4, t(i)+dt(i));
        K6 = dt(i)*func(X(i,:) - 8/27*K1 + 2*K2 - 3544/2565*K3 + 1859/4104*K4 - 11/40*K5, t(i)+dt(i));
        X_X1  = X(i,:) + 25/216 * K1 + 1408/2565*K3 + 2197/4104*K4 - 1/5*K5;
        X = [X; X_X1];
        TE = max(abs(-1/360 * K1 + 128/4275 * K3 + 2197/75240 * K4 - 1/50 * K5 - 2/55 * K6));
        s = 0.9*(epsil/TE)^(1/5);
        if (s*dt(i)+t(i) < tf)
            dt = [dt; s*dt(i)];
        else
            dt = [dt; tf-t(i)];
        end
        t = [t; t(i)+dt(i)];
        i += 1;
    end
    sol = [t dt X];
endfunction