"""
XTC2RDF

A small tool to calculate Radial Distribution Function from .xtc file
for single kind of atoms.

Author: Maksim Posysoev
Date: 24 June 2024, 04 July 2024
License: GPL-3.0

# Usage

julia xtc2rdf.jl -t <trajectory.xtc> -o <output.dat> -n 200 -b 0.05

Distance range: 200 * 0.05 = 10 Å
"""

using Chemfiles
using Printf
using Dates
using ArgParse

const DEFAULT_NBINS = 200
const DEFAULT_BIN_WIDTH = 0.05

"""
    calculate_distance_matrix(frame::Frame) -> Matrix{Float64}

Build a distance matrix for all atoms in the given frame.
"""
function calculate_distance_matrix(frame::Frame)
    num_atoms = length(frame)
    distance_matrix = zeros(Float64, num_atoms, num_atoms)
    @inbounds for i in 1:num_atoms
        for j in 1:(i - 1)
            dist = distance(frame, i - 1, j - 1)
            distance_matrix[i, j] = dist
            distance_matrix[j, i] = dist
        end
    end
    return distance_matrix
end

"""
    accumulate_histogram!(distance_matrix::Matrix{Float64}, histogram::Vector{Float64},
                                        num_atoms::Int, bin_width::Float64, num_bins::Int)

Accumulate histogram of atom-atom distances.
"""
function accumulate_histogram!(
        distance_matrix::Matrix{Float64}, histogram::Vector{Float64},
        num_atoms::Int, bin_width::Float64, num_bins::Int)
    @inbounds for i in 1:num_atoms
        for j in 1:(i - 1)
            bin_index = floor(Int, distance_matrix[i, j] / bin_width) + 1
            if bin_index <= num_bins
                histogram[bin_index] += 2  # Count each pair twice for symmetry
            end
        end
    end
    return histogram
end

"""
    normalize_histogram!(histogram::Vector{Float64}, num_atoms::Int, bin_width::Float64,
                                            num_bins::Int, box_dimensions::Vector{Float64})

Normalize the histogram to produce the RDF.
"""
function normalize_histogram!(
        histogram::Vector{Float64}, num_atoms::Int, bin_width::Float64,
        num_bins::Int, box_dimensions::Vector{Float64})
    box_volume = prod(box_dimensions)
    num_pairs = num_atoms * (num_atoms - 1)
    bins = collect(1:num_bins) .* bin_width
    shell_volumes = 4 * π .* bins .^ 2 .* bin_width
    normalization_factors = box_volume / num_pairs ./ shell_volumes
    histogram .*= normalization_factors
    return histogram
end

"""
    compute_rdf(traj_path::String, num_bins::Int, bin_width::Float64) -> Vector{Float64}

Compute the Radial Distribution Function from a trajectory file.
"""
function compute_rdf(traj_path::String, num_bins::Int, bin_width::Float64)
    trajectory = Trajectory(traj_path)
    num_frames = size(trajectory)
    println("Processing trajectory: $traj_path")
    println("Total frames: $num_frames")

    frame = read_step(trajectory, 0)
    num_atoms = length(frame)
    println("Number of atoms: $num_atoms")

    rdf = zeros(Float64, num_bins)
    histogram = zeros(Float64, num_bins)

    for frame_index in 0:(num_frames - 1)
        frame = read_step(trajectory, frame_index)
        box_dimensions = lengths(UnitCell(frame))
        distance_matrix = calculate_distance_matrix(frame)
        fill!(histogram, 0.0)
        accumulate_histogram!(distance_matrix, histogram, num_atoms, bin_width, num_bins)
        normalize_histogram!(histogram, num_atoms, bin_width, num_bins, box_dimensions)
        rdf .+= histogram
    end

    rdf ./= num_frames
    return rdf
end

"""
    save_rdf(output_path::String, bins::Vector{Float64}, rdf::Vector{Float64})

Write the computed RDF to a file.
"""
function save_rdf(output_path::String, bins::Vector{Float64}, rdf::Vector{Float64})
    @assert length(bins)==length(rdf) "Bins and RDF must have the same length"
    open(output_path, "w") do io
        println(io, "# Radial Distribution Function Data")
        println(io, "# r (Å)    RDF")
        for (r, g) in zip(bins, rdf)
            @printf(io, "%8.4f  %12.6f\n", r, g)
        end
    end
end

"""
    parse_command_line() -> Dict{String, Any}

Parse command-line arguments.
"""
function parse_command_line()
    s = ArgParseSettings(description = "Calculate Radial Distribution Function from XTC trajectory.")

    @add_arg_table! s begin
        "--trajectory", "-t"
        help = "Path to the trajectory file (XTC format)"
        arg_type = String
        required = true
        "--nbins", "-n"
        help = "Number of bins for RDF calculation"
        arg_type = Int
        default = DEFAULT_NBINS
        "--binwidth", "-w"
        help = "Width of each bin for RDF calculation (Å)"
        arg_type = Float64
        default = DEFAULT_BIN_WIDTH
        "--output", "-o"
        help = "Output file name for the RDF data"
        arg_type = String
        default = "rdf_output.dat"
    end

    return parse_args(s)
end

function main()
    start_time = now()
    args = parse_command_line()

    trajectory_path = args["trajectory"]
    n_bins = args["nbins"]
    bin_width = args["binwidth"]
    output_path = args["output"]

    println("XTC -> RDF")
    println("Distance range: from 0 to $(bin_width*n_bins) Å")

    if !isfile(trajectory_path)
        error("Trajectory file not found: $trajectory_path")
    end

    if n_bins <= 0 || bin_width <= 0
        error("n_bins and bin_width should be positive numbers")
    end

    rdf_data = compute_rdf(trajectory_path, n_bins, bin_width)
    bins = collect(1:n_bins) .* bin_width
    save_rdf(output_path, bins, rdf_data)
    println("RDF data written to $output_path")

    end_time = now()
    elapsed_time = canonicalize(end_time - start_time)
    println("Total execution time: $elapsed_time")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
