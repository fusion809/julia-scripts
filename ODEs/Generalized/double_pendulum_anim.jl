using CairoMakie, Dierckx;
function pendpos(r1, r2, theta1, theta2)
        if (length(r1) == 1)
                x1 = r1*cos.(theta1);
                y1 = r1*sin.(theta1);
        else
                x1 = r1.*cos.(theta1);
                y1 = r1.*sin.(theta1);
        end
        if (length(r2) == 1)
                x2 = x1 .+ r2*cos.(theta2);
                y2 = y1 .+ r2*sin.(theta2);
        else
                x2 = x1 .+ r2.*cos.(theta2);
                y2 = y1 .+ r2.*sin.(theta2);
        end
        return x1, y1, x2, y2;
end

function spline_double_pendulum(t, t_uni, x1, y1, x2, y2)
        splx1          = Spline1D(t, x1)
        splx2          = Spline1D(t, x2)
        sply1          = Spline1D(t, y1)
        sply2          = Spline1D(t, y2)
        x1_uni         = evaluate(splx1, t_uni)
        x2_uni         = evaluate(splx2, t_uni)
        y1_uni         = evaluate(sply1, t_uni)
        y2_uni         = evaluate(sply2, t_uni)
        return x1_uni, y1_uni, x2_uni, y2_uni
end

function double_pendulum_anim(t, r1, r2, theta1, theta2, N=10001, padrat = 0.1, filename="Figure 13 Double elastic pendulum animation k is 10 without dissipation.mp4")
        x1, y1, x2, y2 = pendpos(r1, r2, theta1, theta2);
        # Setup animation
        tf             = maximum(t)
        t0             = minimum(t)
        dt             = (tf-t0)/N;
        t_uni          = t0:dt:tf;
        if (length(r1) > 1 || length(r2) > 1)
                padding = padrat*maximum(r1 .+ r2);
        else
                padding = padrat*(r1+r2);
        end
        line1_data     = Observable([[0.0, x1[1]], [0.0, y1[1]]])
        line2_data     = Observable([[x1[1], x2[1]], [y1[1], y2[1]]])
        mass_data      = Observable((x=[x1[1], x2[1]], y=[y1[1], y2[1]]))
        time_text      = Observable("t = 0.00 s")
        f              = CairoMakie.Figure(resolution = (600, 600))
        ax             = Axis(f[1, 1], aspect = 1, 
        limits = (minimum([x1 x2])-padding, maximum([x1 x2])+padding, 
        minimum([y1 y2])-padding, maximum([y1 y2])+padding))
        pendulum_line1 = lines!(ax, line1_data[][1], line1_data[][2], color = :blue)
        pendulum_line2 = lines!(ax, line2_data[][1], line2_data[][2], color = :red)
        masses         = scatter!(ax, mass_data[].x, mass_data[].y, color = [:blue, :red], markersize = 10)
        text!(ax, time_text, position = Point2f(0, 2.1), align = (:left, :top), 
        fontsize = 18, color = :black)
        x1_uni, y1_uni, x2_uni, y2_uni = spline_double_pendulum(t, t_uni, x1, y1, x2, y2)

        # Create and save animation
        record(f, "graphics/$filename", 1:N; framerate = round(Int, 1/dt)) do i
        pendulum_line1[1] = [0, x1_uni[i]]
        pendulum_line1[2] = [0, y1_uni[i]]
        pendulum_line2[1] = [x1_uni[i], x2_uni[i]]
        pendulum_line2[2] = [y1_uni[i], y2_uni[i]]
        masses[1]         = [x1_uni[i], x2_uni[i]]
        masses[2]         = [y1_uni[i], y2_uni[i]]
        time_text[]       = "t = $(round(t_uni[i], digits = 2)) s"
        end
end