function animate(params, radii, angles, t, dt, figno, delay, descriptor)
    if ( (delay != 0) & (size(radii)(2) == size(angles)(2)))
        % Animation
        x = zeros(size(angles));
        y = zeros(size(angles));
        nopen = size(radii)(2);
        for i=1:nopen
            x(:,i) = radii(i) * cos(angles(:,i));
            y(:,i) = radii(i) * sin(angles(:,i));
        end
        % Set up figure
        figure(figno); clf;
        axis equal off;
        axis([-2.2 2.2 -2.2 2.2]);
        hold on;
        % Output GIF filename
        gifname = ["Figure ", num2str(figno), " ", num2str(nopen)];
        gifname = [gifname, " coupled pendulum(s) ", num2str(delay)];
        gifname = [gifname, " delay ", descriptor, ".gif"];

        frame(length(1:length(t))) = struct('cdata', [], 'colormap', []);
        k = 1;
        title=[num2str(nopen), " coupled pendulum(s), delay=", num2str(delay)];
        fields = fieldnames(params);
        for j = 1:numel(fields)
            fname = fields{j};
            val_str = num2str(params.(fname));
            str = [", ", fname, "=", val_str];
            title = [title, str];
        end
        for i = 1:length(t)
            cla;
            plot([0 x(i,1)], [0 y(i,1)], 'LineWidth', 2);      % Rod 1
            for j=2:nopen
                plot([x(i,j-1) x(i,j)], [y(i,j-1) y(i,j)], 'LineWidth', 2);
                plot(x(i,j-1), y(i,j-1), 'bo', 'MarkerSize', 8);
                plot(x(i,j), y(i,j), 'ro', 'MarkerSize', 8);
            end
            
            title([title, sprintf("t = %.2f s", t(i))]);

            drawnow;
            frame(k) = getframe(gcf);
            k = k+1;
        end


        for k=1:delay:length(frame)
            img = frame2im(frame(k));
            [A, map] = rgb2ind(img);

            if k == 1
            imwrite(A, map, gifname, "gif", "LoopCount", Inf, "DelayTime", dt(k));
            else
            imwrite(A, map, gifname, "gif", "WriteMode", "append", 
            "DelayTime", dt(k));
            end
        end

        disp(["GIF saved as ", gifname]);
    elseif (delay == 0)
        disp("Delay is set to zero!")
    else
        disp("Radii argument has more columns than angles argument.")
    end
endfunction