"Can I make a 3D equation model?"

#GLMakie.activate!() #activate GL plotting
# INCLUDE
include("src/JEWEL_functions.jl")
include("src/JEWEL_structures.jl")

######################################################################
######################        INPUTS         #########################
######################################################################
cleanup() # garbage collector and clean memory
######################################################################
## 0. Visual: Plot event IT_DISPLAY time steps
IT_DISPLAY = 25

######################################################################
## 1. Receivers locations
rec_loc = [(13.1, 12, 1.9 + 7); (18.8, 12, 1.9 + 7)]

######################################################################
## 2. SOURCE PARAMETERS:
source_loc = (13.5, 15.5, 0.0 + 7) # location of the source in the domain (x,y,z). All in m.
f0 = 250#1190#1250#830.0                # dominant frequency of the wavelet [Hz]
t0 = 1.20 / f0           # excitation time [s]
factor = -1            # amplitude coefficient of source time function
angle_force = 90.0       # spatial orientation in degrees (90° = x-axis)
#focal SOURCE
strike = 210.0
dip = 89.0
rake = 148.0
Mw = -1.41
Mfault = focal_to_moment_tensor(strike, dip, rake, Mw)

######################################################################
## 3. MODEL PARAMETERS
# Model dimensions [m]
W = 7 + 7 + 7#7.0    # X Dimesion, width of the model in  [m]
L = 35#35.0    # Y Dimesion, lengh of the model in  [m]
H = 20#20.0     # Z Dimesion, height of the model in [m]
# Vp,Vs,rho
α = 5400.0    # Vp in [m/s]
β = α / 1.9   # Vs in [m/s]
ρ = 2800.0   # Density in [kg/m^3]
Qs = 21 # Attenuation parameter

######################################################################
## 4. SIMULATION STABILITY PARAMETERS
CFL = 0.5 # Courant number, should be < 1
GPλ = 11  # Grids per wavelengh (λ) 20
t_total = 0.01  # total recording time [s] 0.02

######################################################################
######################   PREPARE PARAMETRE   #########################
######################################################################

######################################################################
## 5. COMPUTE APPROPIATE dx,dy,dz and number of grid points
cfmax = 1.5   # fmax constant, use 2 for Ricker and 1st Der of Gaussian and 1.5 for Gaussian.
fmax = f0 * cfmax # max frequency based on the wavelet
dz = min((α / (fmax * GPλ)), (β / (fmax * GPλ))) # based on vp and vs velocity, max frequency of the source and GPλ compute dz
dx = dz # same as dz
dy = dz # same as dz
#number of grid points
nx, ny, nz = ceil(Int, W / dx), ceil(Int, L / dy), ceil(Int, H / dz)

######################################################################
# Elastic parameters
vp = fill(α, nx, ny, nz) # compressional wave velocity matrix [m/s]
vs = fill(β, nx, ny, nz)   # shear wave velocity matrix [m/s]
rho = fill(ρ, nx, ny, nz)   # density matrix [m/s]

######################################################################
# Lamé parameters
lam = rho .* (vp .^ 2 .- 2 .* vs .^ 2)    # first Lamé parameter
mu = rho .* vs .^ 2                     # shear modulus [N/m^2]

######################################################################
## TIME STEPPING
dt = CFL * min(dx, dz) / sqrt(maximum(vp)^2 + 2 * maximum(vs)^2) #based on CFL value
nt = round(Int, t_total / dt)        # number of time steps
t = range(0.0, stop=t_total, step=dt)  # time vector


######################################################################
# Source location
isrc, jsrc, ksrc = round.(Int, source_loc ./ (dx, dy, dz))
src0 = (isrc, jsrc, ksrc) #source location vector


######################################################################
## RECEIVERS LOCATIONS AND INITIATE
irec_loc = Vector{Tuple{Int64,Int64,Int64}}(undef, size(rec_loc, 1)) #begin the vector of index location
for ri = 1:size(rec_loc, 1)
    irec_loc[ri] = round.(Int, rec_loc[ri] ./ (dx, dy, dz))
end

# Create the vector of Receivers
receivers = [Receiver(rec_loc[a], irec_loc[a], zeros(nt, 3), zeros(nt, 3), zeros(nt, 3)) for a in 1:size(rec_loc, 1)]


######################################################################
a = π^2 * f0^2  # We don't know
dt2rho_src = dt^2 / rho[ksrc, jsrc, isrc] #Density at the source
min_wavelengh = minimum(vs[findall(vs .> 0.1)]) / f0     # shortest wavelength bounded by velocity in the air

######################################################################
# Source time function

source_signal = factor .* exp.(-a .* (t .- t0) .^ 2);                                # Gaussian wavelet
#source_signal =  -factor .* 2.0 .* a .*(t .- t0) .* exp.( -a .* ( t.- t0) .^2);                    # First derivative of a Gaussian wavelet
#source_signal = -factor .* (1 .- 2 .* a .* (t .- t0) .^ 2) .* exp.(-a .* (t .- t0) .^ 2) # Ricker wavelet: second derivative of Gaussian

######################################################################
# Force components (normalized by volume and scaled by direction)
vol = dx * dy * dz
force_x = sind(angle_force) .* source_signal .* dt2rho_src / vol
force_y = cosd(angle_force) .* source_signal .* dt2rho_src / vol
force_z = sind(angle_force) .* source_signal .* dt2rho_src / vol
# Set forces to zero. Example: force_z if using only horizontal force
force_x .= 0.0
force_y .= 0.0
#force_z .= 0.0

######################################################################
# Moment tensor setup
# mt = MomentTensor(
#     0,           # flag: 0 = force source, 1 = moment tensor
#     1.0,         # Mxx
#     1.0,         # Myy
#     1.0,         # Mzz
#     0.0,         # Mxy
#     0.0,         # Mxz
#     0.0,         # Myz
#     factor,      # scaling factor
#     source_signal,  # source time function
#     (dt2rho_src / vol) # source_physical_scale, source density
# )

mt = MomentTensor(
    1,           # flag: 0 = force source, 1 = moment tensor
    Mfault.Mxx,         # Mxx
    Mfault.Myy,         # Myy
    Mfault.Mzz,         # Mzz
    Mfault.Mxy,         # Mxy
    Mfault.Mxz,         # Mxz
    Mfault.Myz,         # Myz
    factor,      # scaling factor
    source_signal,  # source time function
    (dt2rho_src / vol) # source_physical_scale, source density
)
######################################################################
# DISTRIBUTED SOURCE OVER A SPHERE
# Radius of the spherical source [grid nodes]
szb = ceil(Int, min(nx, nz) ÷ (20 * 2)) # TO CHECK
# Create a 3D array with ones forming a sphere
sphere_b = sphere_mask(szb, gaussian=false)

#make the souce model as large as the medium
sphere_e = zeros(nx, ny, nz)
# Indices where the sphere will be inserted
xrange = (isrc-szb):(isrc+szb)
zrange = (ksrc-szb):(ksrc+szb)
yrange = (jsrc-szb):(jsrc+szb)
#insert source in the model
sphere_e[xrange, yrange, zrange] .= sphere_b

# Distance from source location to each point in the source volume
dist = zeros(size(sphere_e))

@inbounds for i in 1:nx, j in 1:ny, k in 1:nz
    dist[i, j, k] = euclidean((i, j, k), src0)
end

# Exponential source amplitude decay
dist4pr = exp.(-(dist .^ 2) ./ 2)
dist4pr = dim_plus_2(dist4pr) # we padd the array to have the same size as the derivatives

# Attenution = exp.( (-pi*f0*.dist) ./ (vs.*fill(Qs, size(lam))) )
# Attenution = exp.( (-pi*f0*.dt) ./ (fill(Qs, size(lam))) )

######################################################################
# SET BOUNDARY CONDITIONS
weights = SetBoundaryConditions(nx, ny, nz; top=:absorb, bottom=:absorb, right=:reflect, left=:reflect, front=:reflect, back=:absorb)

######################################################################
# Stop here for report 1
print_report_I()
#continue

######################################################################
#Constants for the derivatives
const co_dx = 1 / (2 * dx)
const co_dy = 1 / (2 * dy)
const co_dz = 1 / (2 * dz)

const co_dxx = 1 / dx^2
const co_dyy = 1 / dy^2
const co_dzz = 1 / dz^2
const co_dxy = 1 / (4.0 * dx * dy)
const co_dxz = 1 / (4.0 * dx * dz)
const co_dyz = 1 / (4.0 * dy * dz)

######################################################################
# RUN ELASTIC SIMULATION
######################################################################

# JEWEL_sim!(
#     rho, lam, mu,
#     dt, dx, dy, dz,
#     weights, nt,
#     dist4pr, force_x, force_y, force_z,
#     mt, receivers,
#     IT_DISPLAY, t
# )


######################################################################
# RUN ANAELASTIC SIMULATION
######################################################################

JEWEL_wA_sim!(
    rho, lam, mu, Qs,
    dt, dx, dy, dz,
    weights, nt,
    dist4pr, force_x, force_y, force_z,
    mt, receivers,
    IT_DISPLAY, t, f0
)

######################################################################
# Stop here to make a video
pisc2mp4()
# Plot the wavfileds
plot_wavefield_components(receivers, t)
###########################
# Final Clean Update
#cleanup() # garbage collector and clean memory
