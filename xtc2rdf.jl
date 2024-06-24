"""
XTC2RDF

A small tool to calculate Radial Distribution Function from .xtc file 
for single kind of atoms.

Author: Maksim Posysoev
Date: 24 June 2024
License: GPL-3.0
"""

using Chemfiles
using Printf
using Dates
using ArgParse

# Function to build distance matrix
function buildDistanceMatrix(frame)
    N = length(frame)
    distanceMatrix = zeros(Float64, N, N)
    @inbounds for i in 0:N-1
        @inbounds for j in 0:N-1
            distanceMatrix[i+1, j+1] = distance(frame, i, j)
        end
    end
    return distanceMatrix
end

# Function to accumulate histogram
function hist!(distanceMatrix, hist, N, binWidth, Nbins)
    for i in 1:N
        @fastmath for j in 1:i-1
            histIndex = floor(Int, 1 + distanceMatrix[i, j] / binWidth)
            if histIndex <= Nbins
                hist[histIndex] += 1
            end
        end
    end
    return hist
end

# Function to normalize histogram
function normalizehist!(hist, N, binWidth, Nbins, box)
    boxVolume = box[1] * box[2] * box[3]
    Npairs::Int = N * (N - 1) / 2
    bins = [bin * binWidth for bin = 1:Nbins]
    shellVolumes = [4 * π * binWidth * bins[i]^2 for i in eachindex(bins)]
    rdfNorm = ones(Float64, Nbins)
    for i in eachindex(rdfNorm)
        rdfNorm[i] = boxVolume / Npairs / shellVolumes[i]
    end
    hist .*= rdfNorm
    return hist
end

# Function to compute RDF
function computeRDF(trajname, Nbins, binWidth)
    # Start the timer
    startTime = Dates.now()
    println("Starting at: ", startTime)
    
    # Read trajectory
    traj = Trajectory(trajname)
    nframes = Int(size(traj))
    println("Reading $(trajname)\n")
    println("The trajectory contains $(nframes) frames")
    frame = read_step(traj, 0)
    N = length(frame)
    println("Number of atoms: $(N)")

    # Initialize the RDF array
    RDF = zeros(Float64, Nbins)

    # Accumulate histogram
    for i in 0:nframes-1
        frame = read_step(traj, i)
        box = lengths(UnitCell(frame))
        distanceMatrix = buildDistanceMatrix(frame)
        hist = zeros(Float64, Nbins)
        hist!(distanceMatrix, hist, N, binWidth, Nbins)
        normalizehist!(hist, N, binWidth, Nbins, box)
        # Accumulate normalized histogram
        RDF .+= hist
    end

    # Normalize RDF by the number of frames
    RDF ./= nframes

    # Stop the timer
    stopTime = Dates.now()
    wallTime = Dates.canonicalize(stopTime - startTime)
    println("Stopping at: ", stopTime, "\n")
    println("Walltime: ", wallTime)

    return RDF
end

# Function to write RDF to a file
function writeRDF(outname, bins, rdf, slicing=1)
    @assert length(bins) == length(rdf)
    open(outname, "w") do io
        println(io, "# RDF data")
        println(io, "# r, Å;   \t\tRDF")
        for i in eachindex(bins[1:slicing:end])
            println(io, @sprintf("%6.7f\t%12.7f", bins[1:slicing:end][i], rdf[1:slicing:end][i]))
        end
    end
end

# Main function to parse arguments and run the computation
function main()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--trajname"
        help = "Path to the trajectory file (XTC format)"
        arg_type = String

        "--Nbins"
        help = "Number of bins for RDF calculation"
        arg_type = Int
        default = 200

        "--binWidth"
        help = "Width of each bin for RDF calculation"
        arg_type = Float64
        default = 0.05

        "--outname"
        help = "Output file name for the RDF data"
        arg_type = String
        default = "RDF_output.rdf"
    end

    args = parse_args(s)

    trajname = args["trajname"]
    Nbins = args["Nbins"]
    binWidth = args["binWidth"]
    outname = args["outname"]

    rdf_data = computeRDF(trajname, Nbins, binWidth)
    bins = [bin * binWidth for bin in 1:Nbins]
    writeRDF(outname, bins, rdf_data)
end

# Run the main function if the script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
