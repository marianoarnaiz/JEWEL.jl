"JEWEL.jl's Structures and colormaps File"


using GLMakie, ColorSchemes, Colors

######################################################################
# STRUCTURES

# MomentTensor: Structure for the moment tensor of a source
struct MomentTensor
    flag::Int # Flag: tensor source = 1, use force source = 0
    Mxx::Float64 # xx Moment tensor components
    Myy::Float64 # yy Moment tensor components
    Mzz::Float64 # zz Moment tensor components
    Mxy::Float64 # xy Moment tensor components
    Mxz::Float64 # xz Moment tensor components
    Myz::Float64 # yz Moment tensor components
    factor::Float64  # amplitude coefficient
    source_signal::Vector{Float64} # source time function vector
    source_physical_scale::Float64
end



# Model: Structure for the model parameters! All of them!
struct Model
    #Model Dimesions
    W::Float64
    H::Float64
    L::Float64
    #Model steps
    dx::Float64
    dy::Float64
    dz::Float64
    #Model numbers of grid points for each dimension
    nx::Int
    ny::Int
    nz::Int
    # Geophysical properties arrays
    vp::Array{Float64,3}
    vs::Array{Float64,3}
    rho::Array{Float64,3}
    #Mechanical properties arrays
    lam::Array{Float64,3}
    mu::Array{Float64,3}
    #Bounrady conditions
    weights::Array{Float64,3}
    dist4pr::Array{Float64,3}
end


# MomentTensor: Structure for the moment tensor of a source
struct Receiver
    loc::Tuple{Float64,Float64,Float64} # xx Moment tensor components
    iloc::Tuple{Int64,Int64,Int64} # xx Moment tensor components
    u::Matrix{Float64} # 3 column array for displacement (u) record
    v::Matrix{Float64} # 3 column array for velocity (v) record
    a::Matrix{Float64} # 3 column array for acceleration (v) record
end


######################################################################
# COLORMAPS
# Extract the 11 colors of RdBu_11
cmap = ColorSchemes.RdBu_11
# Copy the colors to a mutable array
colors = copy(cmap.colors)
# Replace the middle color (6th) with pure white
colors[6] = RGBf(1, 1, 1)  # or just colorant"white"
# Create a new continuous colormap with these modified colors
wave_cmap = cgrad(colors, categorical=true)
