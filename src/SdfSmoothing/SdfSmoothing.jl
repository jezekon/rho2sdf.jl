module SdfSmoothing

export RBFs_smoothing, export_rbf_results, rbf_sample, extended_raw_sdf

using JLD2
using SparseArrays
using IterativeSolvers
using FastGaussQuadrature
using Printf
using LinearAlgebra
using StaticArrays
using Base.Threads: @threads

using Rho2sdf.MeshGrid
using Rho2sdf.SignedDistances: fill_far_distances, evalDistances, Sign_Detection,
    find_local_coordinates
using Rho2sdf.DataExport
using Rho2sdf

# Signed distance of the design domain, for the clip of the RBF field (DomainDistance):
include("DomainDistance.jl")

# Raw SDF of the iso-contour alone, extended beyond the design domain (extended_raw_sdf):
include("ExtendedRawSdf.jl")

# Volume of the body defined by the RBF field (RBFBandVolume):
include("CalcVolumeFromSDF.jl")

# RBF smoothing of the SDF (interpolation or approximation) and its export:
include("RBFs4Smoothing.jl")


end


