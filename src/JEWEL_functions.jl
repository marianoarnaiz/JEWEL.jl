"JEWEL.jl's FUNCTIONS FILE"

## Using
using ImageMorphology, OffsetArrays, Distances, CairoMakie, FFMPEG, Printf, GLMakie, Glob

######################################################################
#Clean up function
function cleanup()
    dir = "pics"

    if isdir(dir)
        png_files = glob("*.png", dir)
        for file in png_files
            rm(file, force=true)
        end
    else
        @warn "Directory 'pics' does not exist."
    end

    GC.gc()  # Trigger garbage collection
end


######################################################################

function sphere_mask(radius::Int; gaussian::Bool=false)
    d = 2 * radius + 1
    center = radius + 1

    if gaussian
        σ² = (radius / 2)^2
        sphere = zeros(Float32, d, d, d)

        for i in 1:d, j in 1:d, k in 1:d
            r² = (i - center)^2 + (j - center)^2 + (k - center)^2
            sphere[i, j, k] = exp(-r² / (2 * σ²))
        end

    else
        sphere = zeros(Bool, d, d, d)

        for i in 1:d, j in 1:d, k in 1:d
            if (i - center)^2 + (j - center)^2 + (k - center)^2 ≤ radius^2
                sphere[i, j, k] = true
            end
        end
    end

    return sphere
end

######################################################################

function print_report_I()
    # REPORT 1
    # Colors and headers
    printstyled("#################################################\n", color=:cyan, bold=true)
    printstyled("3D elastic FDTD wave propagation in isotropic medium\n", color=:green, bold=true)
    printstyled("Displacement formulation with Cerjan (1985) boundary conditions\n", color=:green)
    printstyled("#################################################\n", color=:cyan, bold=true)

    # Model info
    printstyled("Model 🪨:\n", bold=true, color=:green)
    println("\t$(nz) x $(ny) x $(nx)\tgrid (nz x ny x nx)")
    println("\t$(round(dz, sigdigits=2)) x $(round(dy, sigdigits=2)) x $(round(dx, sigdigits=2)) [m] (dz x dy x dx)")
    println("\t$(round(nx*dx, sigdigits=3)) x $(round(ny*dy, sigdigits=3)) x $(round(nz*dz, sigdigits=3)) [m] model size")
    println("\t Min: $(round(minimum(vp), sigdigits=3)) → Max: $(round(maximum(vp), sigdigits=3)) Vp [m/s] ")
    println("\t Min: $(round(minimum(vs), sigdigits=3)) → Max: $(round(maximum(vs), sigdigits=3)) Vs [m/s] ")
    println("\t Min: $(round(minimum(rho), sigdigits=3)) →Max: $(round(maximum(rho), sigdigits=3)) rho [kg/m³] ")

    # Time info
    printstyled("Time 🕖:\n", bold=true, color=:blue)
    println("\t$t_total [s] total")
    println("\t$dt [s] dt")
    println("\t$nt time steps")

    # Source info
    printstyled("Source 💥:\n", bold=true, color=:red)
    println("\t$f0 Hz dominant frequency")
    println("\tSource starts at $t0 s")

    # Other
    printstyled("Other:\n", bold=true, color=:yellow)
    println("\t$CFL CFL number")
    println("\t$min_wavelengh [m] shortest wavelength")
    println("\t$GPλ grid points-per-wavelength")



    printstyled("#################################################\n", color=:cyan, bold=true)

end

######################################################################

# Padding function: pad 1 zero in all directions (equivalent to padarray(A,[1 1 1]))
# zero-pad by 1 in every dimension
function dim_plus_2(A::Array{T,3}) where {T}
    pad = zeros(T, size(A) .+ 2)
    pad[2:end-1, 2:end-1, 2:end-1] .= A
    return pad
end

######################################################################
# first‐order centered derivatives (mirror MATLAB indexing)
function d_x(A::Array{T,3}) where {T}
    return co_dx .* (A[2:end-1, 2:end-1, 3:end] .- A[2:end-1, 2:end-1, 1:end-2])
end

function d_y(A::Array{T,3}) where {T}
    return co_dy .* (A[2:end-1, 3:end, 2:end-1] .- A[2:end-1, 1:end-2, 2:end-1])
end

function d_z(A::Array{T,3}) where {T}
    return co_dz .* (A[3:end, 2:end-1, 2:end-1] .- A[1:end-2, 2:end-1, 2:end-1])
end

######################################################################

# Second order derivatives
function d_xx(A::Array{T,3}) where {T}
    return d_x(dim_plus_2(d_x(A)))
end

function d_yy(A::Array{T,3}) where {T}
    return d_y(dim_plus_2(d_y(A)))
end

function d_zz(A::Array{T,3}) where {T}
    return d_z(dim_plus_2(d_z(A)))
end

# Mixed derivatives
function d_xz(A::Array{T,3}) where {T}
    return d_z(dim_plus_2(d_x(A)))
end

function d_xy(A::Array{T,3}) where {T}
    return d_y(dim_plus_2(d_x(A)))
end

function d_yz(A::Array{T,3}) where {T}
    return d_z(dim_plus_2(d_y(A)))
end


######################################################################
# FUNCTION TO SET BOUNDARY CONDITIOINS
function SetBoundaryConditions(nx::Int, ny::Int, nz::Int; top=:reflect, bottom=:reflect, right=:reflect, left=:reflect, front=:reflect, back=:reflect)

    # Set absorbing layer thickness and damping rate
    abs_thick = min(floor(Int, 0.15 * nx), floor(Int, 0.15 * nz))
    abs_rate = 0.3 / abs_thick

    # Padding to support absorbing layers (with +2 padding)
    weights = ones(nx + 2, ny + 2, nz + 2)

    # Define margins per direction
    lmargin = (abs_thick, abs_thick, abs_thick)
    rmargin = (abs_thick, abs_thick, abs_thick)

    # Loop through all padded indices
    for ix in 1:nx+2, iy in 1:ny+2, iz in 1:nz+2
        i = j = k = 0.0

        # X-direction (left-right)
        if ix < lmargin[1] + 1 && left == :absorb
            i = lmargin[1] + 1 - ix
        elseif ix > nx - rmargin[1] && right == :absorb
            i = ix - (nx - rmargin[1])
        end

        # Y-direction (back-front)
        if iy < lmargin[2] + 1 && back == :absorb
            j = lmargin[2] + 1 - iy
        elseif iy > ny - rmargin[2] && front == :absorb
            j = iy - (ny - rmargin[2])
        end

        # Z-direction (bottom-top)
        if iz < lmargin[3] + 1 && bottom == :absorb
            k = lmargin[3] + 1 - iz
        elseif iz > nz - rmargin[3] && top == :absorb
            k = iz - (nz - rmargin[3])
        end

        if i == 0 && j == 0 && k == 0
            continue
        end

        rr = abs_rate^2 * (i^2 + j^2 + k^2)
        weights[ix, iy, iz] = exp(-rr)
    end

    return weights
end


######################################################################
# FUNCTION TO RUN FULLY ELASTIC  SIMULATION
function JEWEL_sim!(
    rho, lam, mu,
    dt, dx, dy, dz,
    weights, nt,
    dist4pr, force_x, force_y, force_z,
    mt, receivers,
    IT_DISPLAY, t
)
    printstyled(" \n Starts FDTD of 3D Elastic Wave Equation \n", color=:green, bold=true)

    # PRE-ALLOCATION OF WAVEFIELD TENSORS
    # +2 for ghost points on all sides
    ux3 = zeros(nx + 2, ny + 2, nz + 2)  # Wavefield at t
    uy3 = zeros(nx + 2, ny + 2, nz + 2)
    uz3 = zeros(nx + 2, ny + 2, nz + 2)
    ux2 = zeros(nx + 2, ny + 2, nz + 2)  # Wavefield at t-1
    uy2 = zeros(nx + 2, ny + 2, nz + 2)
    uz2 = zeros(nx + 2, ny + 2, nz + 2)
    ux1 = zeros(nx + 2, ny + 2, nz + 2)  # Wavefield at t-2
    uy1 = zeros(nx + 2, ny + 2, nz + 2)
    uz1 = zeros(nx + 2, ny + 2, nz + 2)



    # Precomputed constants
    dt2rho = (dt^2) ./ rho
    lam_2mu = lam .+ 2 .* mu
    lam_mu = lam .+ mu


    c = 1
    for it in 1:nt
        printstyled("Working on time step $it of $nt \n", color=:blue, italic=true)

        # Clear current wavefield
        fill!(ux3, 0.0)
        fill!(uy3, 0.0)
        fill!(uz3, 0.0)

        # Second-order derivatives
        dux_dxx = d_xx(ux2)
        dux_dyy = d_yy(ux2)
        dux_dzz = d_zz(ux2)
        dux_dxz = d_xz(ux2)
        dux_dxy = d_xy(ux2)

        duy_dxx = d_xx(uy2)
        duy_dyy = d_yy(uy2)
        duy_dzz = d_zz(uy2)
        duy_dxy = d_xy(uy2)
        duy_dyz = d_yz(uy2)

        duz_dxx = d_xx(uz2)
        duz_dyy = d_yy(uz2)
        duz_dzz = d_zz(uz2)
        duz_dxz = d_xz(uz2)
        duz_dyz = d_yz(uz2)

        # Compute RHS
        sigmas_ux = lam_2mu .* dux_dxx .+ mu .* dux_dyy .+ mu .* dux_dzz .+ lam_mu .* duz_dxz .+ lam_mu .* duy_dxy
        sigmas_uy = mu .* duy_dxx .+ lam_2mu .* duy_dyy .+ mu .* duy_dzz .+ lam_mu .* duz_dyz .+ lam_mu .* dux_dxy
        sigmas_uz = mu .* duz_dxx .+ mu .* duz_dyy .+ lam_2mu .* duz_dzz .+ lam_mu .* duy_dyz .+ lam_mu .* dux_dxz

        # Update displacement fields
        ux3[2:end-1, 2:end-1, 2:end-1] .= 2.0 .* ux2[2:end-1, 2:end-1, 2:end-1] .- ux1[2:end-1, 2:end-1, 2:end-1] .+ sigmas_ux .* dt2rho
        uy3[2:end-1, 2:end-1, 2:end-1] .= 2.0 .* uy2[2:end-1, 2:end-1, 2:end-1] .- uy1[2:end-1, 2:end-1, 2:end-1] .+ sigmas_uy .* dt2rho
        uz3[2:end-1, 2:end-1, 2:end-1] .= 2.0 .* uz2[2:end-1, 2:end-1, 2:end-1] .- uz1[2:end-1, 2:end-1, 2:end-1] .+ sigmas_uz .* dt2rho

        # Source
        if mt.flag == 0
            ux3 .+= dist4pr .* force_x[it]
            uy3 .+= dist4pr .* force_y[it]
            uz3 .+= dist4pr .* force_z[it]
        else
            ux3[2:end-1, 2:end-1, 2:end-1] .+= force_x[it] .* (mt.xx .* d_x(dist4pr) .+ mt.xy .* d_y(dist4pr) .+ mt.xz .* d_z(dist4pr))
            uy3[2:end-1, 2:end-1, 2:end-1] .+= force_x[it] .* (mt.yy .* d_y(dist4pr) .+ mt.xy .* d_x(dist4pr) .+ mt.yz .* d_z(dist4pr))
            uz3[2:end-1, 2:end-1, 2:end-1] .+= force_x[it] .* (mt.zz .* d_z(dist4pr) .+ mt.yz .* d_y(dist4pr) .+ mt.xz .* d_x(dist4pr))
        end

        # Velocity & Acceleration
        vx = @. (ux3 - ux2) / dt
        vy = @. (uy3 - uy2) / dt
        vz = @. (uz3 - uz2) / dt
        ax = @. (ux3 - 2 * ux2 + ux1) / dt^2
        ay = @. (uy3 - 2 * uy2 + uy1) / dt^2
        az = @. (uz3 - 2 * uz2 + uz1) / dt^2

        # Absorbing boundary rotation
        ux1 .= ux2 .* weights
        ux2 .= ux3 .* weights
        uy1 .= uy2 .* weights
        uy2 .= uy3 .* weights
        uz1 .= uz2 .* weights
        uz2 .= uz3 .* weights

        # Receiver Recording
        for ri = 1:size(receivers, 1)
            i, j, k = receivers[ri].iloc
            receivers[ri].u[it, :] .= (ux2[i, j, k], uy2[i, j, k], uz2[i, j, k])
            receivers[ri].v[it, :] .= (vx[i, j, k], vy[i, j, k], vz[i, j, k])
            receivers[ri].a[it, :] .= (ax[i, j, k], ay[i, j, k], az[i, j, k])
        end

        # Optional Plotting

        if mod(it, IT_DISPLAY) == 0

            GLMakie.activate!() #activate GL plotting
            it_time = t[it]
            U = @. sqrt(ux2^2 + uy2^2 + uz2^2)
            fig = Figure(size=(1200, 1200))
            ax4 = Axis3(fig[1, 1], title=title = "Waveform at t = $it_time s \n |U|", xlabel="X [m]", ylabel="Y [m]", zlabel="Z [m]", aspect=:data)
            #hm4 = volume!(ax4, U, colormap=:viridis, transparency=true)
            hm4 = volume!(ax4, (0, W), (0, L), (0, H), U, colormap=:viridis, transparency=true)
            scatter!(ax4,rec_loc, marker=:dtriangle, color=:white)
            scatter!(ax4,source_loc, marker=:star5, color=:red)
            Colorbar(fig[2, 1], hm4, label="|U| [m]", vertical=false)

            # GLMakie.activate!() #activate GL plotting
            # it_time = t[it]
            # U = @. sqrt(ux2^2 + uy2^2 + uz2^2)
            # fig = Figure(size=(1200, 800))
            # ax1 = Axis3(fig[1, 1], title="Waveform at t = $it_time s \n vz", xlabel="X", ylabel="Z")
            # ccamp=maximum(abs.(vz))*0.5
            # hm1=volume!(ax1, vz, colormap=wave_cmap, algorithm=:absorption, absorption=5.0f0, colorrange=(-ccamp,ccamp))
            # Colorbar(fig[1, 1][1, 2], hm1, label="v [m/s]")
            # ax2 = Axis3(fig[1, 2][1, 1], title="vy", xlabel="X", ylabel="Z")
            # ccamp=maximum(abs.(vy))*0.5
            # hm2 = volume!(ax2, vy, colormap=wave_cmap, algorithm=:absorption, absorption=5.0f0, colorrange=(-ccamp,ccamp))
            # Colorbar(fig[1, 2][1, 2], hm2, label="v [m/s]")
            # ax3 = Axis3(fig[2, 1], title="vx", xlabel="X", ylabel="Z")
            # ccamp=maximum(abs.(vx))*0.5
            # hm3=volume!(ax3, vx, colormap=wave_cmap, algorithm=:absorption, absorption=5.0f0, colorrange=(-ccamp,ccamp))
            # Colorbar(fig[2, 1][1, 2], hm3, label="v [m/s]")
            # ax4 = Axis3(fig[2, 2][1, 1], title="|U|", xlabel="X", ylabel="Z")
            # hm4 = volume!(ax4, U, colormap=:viridis)
            # Colorbar(fig[2, 2][1, 2], hm4, label="|U| [m]")

            #save("pics/Wavefild_$it.png", fig, backend=CairoMakie)
            save("pics/$( @sprintf("%04d", c) ).png", fig)
            c = c + 1 #update counter
        end
    end
end


######################################################################
# FUNCTION TO RUN ELASTIC  SIMULATION with Anelastic Attenuation! see graves 1996
function JEWEL_wA_sim!(
    rho, lam, mu, Qs,
    dt, dx, dy, dz,
    weights, nt,
    dist4pr, force_x, force_y, force_z,
    mt, receivers,
    IT_DISPLAY, t, f0
)
    printstyled(" \n Starts FDTD of 3D Elastic Wave Equation \n", color=:green, bold=true)


    # PRE-ALLOCATION OF WAVEFIELD TENSORS
    # +2 for ghost points on all sides
    ux3 = zeros(nx + 2, ny + 2, nz + 2)  # Wavefield at t
    uy3 = zeros(nx + 2, ny + 2, nz + 2)
    uz3 = zeros(nx + 2, ny + 2, nz + 2)
    ux2 = zeros(nx + 2, ny + 2, nz + 2)  # Wavefield at t-1
    uy2 = zeros(nx + 2, ny + 2, nz + 2)
    uz2 = zeros(nx + 2, ny + 2, nz + 2)
    ux1 = zeros(nx + 2, ny + 2, nz + 2)  # Wavefield at t-2
    uy1 = zeros(nx + 2, ny + 2, nz + 2)
    uz1 = zeros(nx + 2, ny + 2, nz + 2)

    # Precomputed constants
    dt2rho = (dt^2) ./ rho
    lam_2mu = lam .+ 2 .* mu
    lam_mu = lam .+ mu

    #precompute Attenuations tensor
    #Attenution = exp.(-pi*f0*dt ./fill(Qs, size(lam_2mu)))
    Attenution = exp((-pi * f0 * dt) / Qs)
    c = 1
    for it in 1:nt
        printstyled("Working on time step $it of $nt \n", color=:blue, italic=true)

        # Clear current wavefield
        fill!(ux3, 0.0)
        fill!(uy3, 0.0)
        fill!(uz3, 0.0)

        # Second-order derivatives
        dux_dxx = d_xx(ux2)
        dux_dyy = d_yy(ux2)
        dux_dzz = d_zz(ux2)
        dux_dxz = d_xz(ux2)
        dux_dxy = d_xy(ux2)

        duy_dxx = d_xx(uy2)
        duy_dyy = d_yy(uy2)
        duy_dzz = d_zz(uy2)
        duy_dxy = d_xy(uy2)
        duy_dyz = d_yz(uy2)

        duz_dxx = d_xx(uz2)
        duz_dyy = d_yy(uz2)
        duz_dzz = d_zz(uz2)
        duz_dxz = d_xz(uz2)
        duz_dyz = d_yz(uz2)

        # Compute RHS
        sigmas_ux = (lam_2mu .* dux_dxx .+ mu .* dux_dyy .+ mu .* dux_dzz .+ lam_mu .* duz_dxz .+ lam_mu .* duy_dxy)# .* Attenution
        sigmas_uy = (mu .* duy_dxx .+ lam_2mu .* duy_dyy .+ mu .* duy_dzz .+ lam_mu .* duz_dyz .+ lam_mu .* dux_dxy)# .* Attenution
        sigmas_uz = (mu .* duz_dxx .+ mu .* duz_dyy .+ lam_2mu .* duz_dzz .+ lam_mu .* duy_dyz .+ lam_mu .* dux_dxz)# .* Attenution

        # Update displacement fields
        ux3[2:end-1, 2:end-1, 2:end-1] .= (2.0 .* ux2[2:end-1, 2:end-1, 2:end-1] .- ux1[2:end-1, 2:end-1, 2:end-1] .+ sigmas_ux .* dt2rho) .* Attenution
        uy3[2:end-1, 2:end-1, 2:end-1] .= (2.0 .* uy2[2:end-1, 2:end-1, 2:end-1] .- uy1[2:end-1, 2:end-1, 2:end-1] .+ sigmas_uy .* dt2rho) .* Attenution
        uz3[2:end-1, 2:end-1, 2:end-1] .= (2.0 .* uz2[2:end-1, 2:end-1, 2:end-1] .- uz1[2:end-1, 2:end-1, 2:end-1] .+ sigmas_uz .* dt2rho) .* Attenution

        # Source
        if mt.flag == 0
            ux3 .+= dist4pr .* force_x[it]
            uy3 .+= dist4pr .* force_y[it]
            uz3 .+= dist4pr .* force_z[it]
        else
            stfamp = mt.source_signal[it]*mt.source_physical_scale   #source time function amplitud
            ux3[2:end-1, 2:end-1, 2:end-1] .+= stfamp .* (mt.Mxx .* d_x(dist4pr) .+ mt.Mxy .* d_y(dist4pr) .+ mt.Mxz .* d_z(dist4pr))
            uy3[2:end-1, 2:end-1, 2:end-1] .+= stfamp .* (mt.Myy .* d_y(dist4pr) .+ mt.Mxy .* d_x(dist4pr) .+ mt.Myz .* d_z(dist4pr))
            uz3[2:end-1, 2:end-1, 2:end-1] .+= stfamp .* (mt.Mzz .* d_z(dist4pr) .+ mt.Myz .* d_y(dist4pr) .+ mt.Mxz .* d_x(dist4pr))
        end

        # Velocity & Acceleration
        vx = @. (ux3 - ux2) / dt
        vy = @. (uy3 - uy2) / dt
        vz = @. (uz3 - uz2) / dt
        ax = @. (ux3 - 2 * ux2 + ux1) / dt^2
        ay = @. (uy3 - 2 * uy2 + uy1) / dt^2
        az = @. (uz3 - 2 * uz2 + uz1) / dt^2

        # Absorbing boundary rotation
        ux1 .= ux2 .* weights
        ux2 .= ux3 .* weights
        uy1 .= uy2 .* weights
        uy2 .= uy3 .* weights
        uz1 .= uz2 .* weights
        uz2 .= uz3 .* weights

        # Receiver Recording
        for ri = 1:size(receivers, 1)
            i, j, k = receivers[ri].iloc
            receivers[ri].u[it, :] .= (ux2[i, j, k], uy2[i, j, k], uz2[i, j, k])
            receivers[ri].v[it, :] .= (vx[i, j, k], vy[i, j, k], vz[i, j, k])
            receivers[ri].a[it, :] .= (ax[i, j, k], ay[i, j, k], az[i, j, k])
        end

        # Optional Plotting

        if mod(it, IT_DISPLAY) == 0

            GLMakie.activate!() #activate GL plotting
            it_time = t[it]
            U = @. sqrt(ux2^2 + uy2^2 + uz2^2)
            fig = Figure(size=(1200, 1200))
            ax4 = Axis3(fig[1, 1], title=title = "Waveform at t = $it_time s \n |U|", xlabel="X [m]", ylabel="Y [m]", zlabel="Z [m]", aspect=:data)
            #hm4 = volume!(ax4, U, colormap=:viridis, transparency=true)
            hm4 = volume!(ax4, (0, W), (0, L), (0, H), U, colormap=:viridis, transparency=true)
            scatter!(ax4,rec_loc, marker=:dtriangle, color=:white)
            scatter!(ax4,source_loc, marker=:star5, color=:red)
            Colorbar(fig[2, 1], hm4, label="|U| [m]", vertical=false)

            # GLMakie.activate!() #activate GL plotting
            # it_time = t[it]
            # U = @. sqrt(ux2^2 + uy2^2 + uz2^2)
            # fig = Figure(size=(1200, 800))
            # ax1 = Axis3(fig[1, 1], title="Waveform at t = $it_time s \n vz", xlabel="X", ylabel="Z")
            # ccamp=maximum(abs.(vz))*0.5
            # hm1=volume!(ax1, vz, colormap=wave_cmap, algorithm=:absorption, absorption=5.0f0, colorrange=(-ccamp,ccamp), transparency=true)
            # Colorbar(fig[1, 1][1, 2], hm1, label="v [m/s]")
            # ax2 = Axis3(fig[1, 2][1, 1], title="vy", xlabel="X", ylabel="Z")
            # ccamp=maximum(abs.(vy))*0.5
            # hm2 = volume!(ax2, vy, colormap=wave_cmap, algorithm=:absorption, absorption=5.0f0, colorrange=(-ccamp,ccamp), transparency=true)
            # Colorbar(fig[1, 2][1, 2], hm2, label="v [m/s]")
            # ax3 = Axis3(fig[2, 1], title="vx", xlabel="X", ylabel="Z")
            # ccamp=maximum(abs.(vx))*0.5
            # hm3=volume!(ax3, vx, colormap=wave_cmap, algorithm=:absorption, absorption=5.0f0, colorrange=(-ccamp,ccamp), transparency=true)
            # Colorbar(fig[2, 1][1, 2], hm3, label="v [m/s]")
            # ax4 = Axis3(fig[2, 2][1, 1], title="|U|", xlabel="X", ylabel="Z")
            # hm4 = volume!(ax4, U, colormap=:viridis)
            # Colorbar(fig[2, 2][1, 2], hm4, label="|U| [m]")

            #save("pics/Wavefild_$it.png", fig, backend=CairoMakie)
            save("pics/$( @sprintf("%04d", c) ).png", fig)
            c = c + 1 #update counter
        end
    end
end





######################################################################
# A function that makes a video
function pisc2mp4(imagesdirectory::String="pics", framerate::Int=2, outputfile::String="Wavefield.mp4")
    ffmpeg_cmd = `-framerate $(framerate) -f image2 -i $(imagesdirectory)/%04d.png -vf scale=1200:1200 -c:v libx264 -pix_fmt yuv420p -y $(outputfile)`
    FFMPEG.ffmpeg_exe(ffmpeg_cmd)
end


function plot_wavefield_components(receivers, t)
jldsave("receivers_20GPWL.jld2",receivers=receivers)
    for ri in 1:size(receivers, 1)

        right_size = size(receivers[ri].u[:, 1], 1)
        fig = Figure(size=(800, 1200))

        # Displacement
        ax1 = Axis(fig[1, 1], ylabel="u [m]", title="RECORDS OF: Receiver $(ri) \n Displacement Records")
        lines!(ax1, t[1:right_size], receivers[ri].u[:, 1], color=:red, label="Ux")
        lines!(ax1, t[1:right_size], receivers[ri].u[:, 2], color=:green, label="Uy")
        lines!(ax1, t[1:right_size], receivers[ri].u[:, 3], color=:blue, label="Uz")
        axislegend(ax1, position=:lb)

        # Velocity
        ax2 = Axis(fig[2, 1], ylabel="v [m/s]", title="Velocity Records")
        lines!(ax2, t[1:right_size], receivers[ri].v[:, 1], color=:red, label="vx")
        lines!(ax2, t[1:right_size], receivers[ri].v[:, 2], color=:green, label="vy")
        lines!(ax2, t[1:right_size], receivers[ri].v[:, 3], color=:blue, label="vz")
        axislegend(ax2, position=:lb)

        # Acceleration
        ax3 = Axis(fig[3, 1], xlabel="t [s]", ylabel="a [m/s²]", title="Acceleration Records")
        lines!(ax3, t[1:right_size], receivers[ri].a[:, 1], color=:red, label="ax")
        lines!(ax3, t[1:right_size], receivers[ri].a[:, 2], color=:green, label="ay")
        lines!(ax3, t[1:right_size], receivers[ri].a[:, 3], color=:blue, label="az")
        axislegend(ax3, position=:lb)


        save("receiver_$(ri).png", fig)
    end
end







######################################################################
# FOCAL TO MOMENT TENSOR
function focal_to_moment_tensor(strike::Float64, dip::Float64, rake::Float64, Mw::Float64)
    # Convert degrees to radians
    strike = deg2rad(strike)
    dip = deg2rad(dip)
    rake = deg2rad(rake)

    # Compute scalar seismic moment in dyne·cm, convert to N·m (1e7 dyne·cm = 1 N·m)
    M0 = 10.0^(1.5 * Mw + 9.1) / 1e7

    # Trig functions
    cosλ, sinλ = cos(rake), sin(rake)
    cosφ, sinφ = cos(strike), sin(strike)
    cosδ, sinδ = cos(dip), sin(dip)

    # Moment tensor components in NED coordinate system
    Mrr = -M0 * (sinδ * cosλ * sin(2 * strike) + sin(2 * dip) * sinλ * sin(strike)^2)
    Mtt = M0 * (sinδ * cosλ * sin(2 * strike) - sin(2 * dip) * sinλ * cos(strike)^2)
    Mpp = -Mrr - Mtt
    Mrt = -M0 * (sinδ * cosλ * cos(2 * strike) + 0.5 * sin(2 * dip) * sinλ * sin(2 * strike))
    Mrp = -M0 * (cosδ * cosλ * cos(strike) + cos(2 * dip) * sinλ * sin(strike))
    Mtp = -M0 * (cosδ * cosλ * sin(strike) - cos(2 * dip) * sinλ * cos(strike))

    # Return Cartesian x-y-z components (assumes z=down, y=north, x=east)
    return (
        Mxx=Mtt,
        Myy=Mpp,
        Mzz=Mrr,
        Mxy=Mtp,
        Mxz=Mrp,
        Myz=Mrt
    )
end
