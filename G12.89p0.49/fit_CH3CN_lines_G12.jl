using LsqFit
using Optim
using LinearAlgebra
using Random
using Distributions
using MCMCDiagnosticTools
using KernelDensity
using ThreadTools
using Serialization
using DelimitedFiles
using FITSIO
using Statistics
using StatsBase
using DataStructures
#using PairPlots
using CairoMakie
using LaTeXStrings
#using BlackBoxOptim

# example of command to run the script:
# julia -t 14 --heap-size-hint=80% --check-bounds=no fit_CH3CN_lines_G12.jl &> log.txt 

Random.seed!(11051989)

const DATA_DIR = "G12"

const FILE_NAME_WITH_MODELS_GRID = joinpath(DATA_DIR, "grid_results_CH3CN_G12.txt")

const model_ids = open(joinpath(DATA_DIR, "modelids_G12.jls"), "r") do file
    deserialize(file)
end

const MOLECULAR_MASS = 41.05 # g/mol

const THRESHOLD_IN_SIGMAS = 5
const MIN_LINE_FWHM_IN_KMS = 2.5
const MAX_LINE_FWHM_IN_KMS = 5.0
const VLSR_IN_KMS = 33.8
const MIN_ALLOWED_LINE_VEL_WRT_VLSR_KMS = -8
const MAX_ALLOWED_LINE_VEL_WRT_VLSR_KMS = 8
#const LINE_VERTICAL_LABELS = OrderedDict("12_{11}->11_{11}" => 220.2350315, "12_{10}->11_{10}" => 220.3236311, "12_9->11_9" => 220.4039006, "12_8->11_8" => 220.4758078, "12_7->11_7" => 220.539324, "12_6->11_6" => 220.5944237, "12_5->11_5" => 220.6410844, "12_4->11_4" => 220.6792874, "12_3->11_3" => 220.709017, "12_2->11_2" => 220.7302611, "12_1->11_1" => 220.7430111, "12_0->11_0" => 220.7472617)
const LINE_VERTICAL_LABELS = OrderedDict("12_8->11_8" => 220.4758078, "12_7->11_7" => 220.539324, "12_6->11_6" => 220.5944237, "12_5->11_5" => 220.6410844, "12_4->11_4" => 220.6792874, "12_3->11_3" => 220.709017, "12_2->11_2" => 220.7302611, "12_1->11_1" => 220.7430111, "12_0->11_0" => 220.7472617)
const INDEXES_OF_LINES_THAT_WILL_BE_SELECTED_BASED_ON_VLSR_BUT_NOT_ON_INTENSITY = [1, 3, 4, 8, 9]
const MODEL_FREQ_INDEXES_TO_NOT_INCLUDE_IN_ABSENT_LINE_LIST = [8, 9]
const MODEL_FREQ_INDEXES_TO_NOT_INCLUDE_IN_PHYSICAL_COND_ESTIMATION = [8, 9]
const MODEL_FREQ_INDEX_FOR_OPTICAL_DEPTH_OUTPUT = 5

const FITS_FILE = "line27_targetImageForSelfcal.image.pbcor_cut_channel.fits"
const FITS_FOR_NOISE_ESTIMATES = "line27_targetImageForSelfcal.image.pbcor_cut_all_channels.fits"
const CONT_CHANNELS_FOR_NOISE_ESTIMATES = "50~57;98~99;189~190;231~236;316~326;399~405;524~526;641~643;710~714;769~773;788~796;803~811;842~843;883~886;925~936;961~981;1010~1012;1047~1067;1175~1191;1209~1212;1231~1233;1269~1271;1343~1350;1409~1415;1423~1431;1493~1497;1626~1629;1666~1667;1695~1697;1735~1736;1756~1766;1776~1780;1813~1822;1853~1861;1881~1890;1899~1910"

const HDENSITIES0 = collect(3.0:0.02:10.0) # [log cm^-3]
const GAS_TEMPERATURES = collect(10:1:950) # [K]
const SPECIFIC_COLUMN_DENSITIES = collect(8.5:0.02:14.08) # [log cm^-3 s]
const MOLECULAR_ABUNDANCE = 1.e-6 # wrt H2
const MINIMUM_CLOUD_SIZE_PROJECTED_ON_SKY = 100.0 * 1.496e+13 # [cm]
const MAXIMUM_CLOUD_SIZE_PROJECTED_ON_SKY = 100000. * 1.496e+13 # [cm]
const MAX_SPECIFIC_COLUMN_DENSITY = 13.0 # [log cm^-3 s]

const Au = [0.0005116329002001736, 0.0006082282208606637, 0.0006918777912929745, 0.0007627413847878795, 0.0008209859267298996, 0.0008662156493747154, 0.0008983886789123312, 0.0009178310034746356, 0.0009243321150775603]
const gu = [50., 50., 100., 50., 50., 100., 50., 50., 50.]
const Eu_cm = [365.2936, 290.965, 226.5128, 171.9517, 127.2938, 92.5493, 67.7262, 52.8301, 47.8643]

function get_partition_function(T)
    return 0.48597 * T^1.5929
end

##############################################################################################

if !isdir(DATA_DIR)
    mkdir(DATA_DIR)
end

if !isdir(joinpath(DATA_DIR, "figures/"))
    mkdir(joinpath(DATA_DIR, "figures/"))
end

if !isdir(joinpath(DATA_DIR, "figures/spectra/"))
    mkdir(joinpath(DATA_DIR, "figures/spectra/"))
end

if !isdir(joinpath(DATA_DIR, "figures/fitspectra/"))
    mkdir(joinpath(DATA_DIR, "figures/fitspectra/"))
end

if !isdir(joinpath(DATA_DIR, "figures/rot_diagram/"))
    mkdir(joinpath(DATA_DIR, "figures/rot_diagram/"))
end

const c_cgs = 2.99792458e10  # cm/s
const h_cgs = 6.62607015e-27  # erg·s
const k_cgs = 1.380649e-16  # erg/K

const Eu_K = Eu_cm .* (h_cgs * c_cgs / k_cgs)

const CLOUD_SIZE_FACTOR = 1.e5 / MOLECULAR_ABUNDANCE

const plotting_lock = ReentrantLock()

function parse_number_ranges(s::String)
    ranges = split(s, ";")
    numbers = Int[]
    for r in ranges
        start, stop = map(x -> parse(Int, x), split(r, "~"))
        push!(numbers, start:stop...)
    end
    return numbers
end

const CONT_CHANNELS_FOR_NOISE_ESTIMATES_VECTOR = parse_number_ranges(CONT_CHANNELS_FOR_NOISE_ESTIMATES)

function sample(llhood::Function, numwalkers::Integer, x0::Matrix{Float64}, numsamples_perwalker::Integer, thinning::Integer=numwalkers, a::Number=2.; rng::Random.AbstractRNG=Random.GLOBAL_RNG)
	if numsamples_perwalker < 2
		numsamples_perwalker = 2
	end
	x = copy(x0)
	chain = Array{Float64}(undef, div(numsamples_perwalker, thinning), numwalkers, size(x0, 2))
    lastllhoodvals = zeros(Float64, numwalkers)
    for i = 1:numwalkers
        lastllhoodvals[i] = llhood(@view(x[i, :]))
    end
	chain[1, :, :] = x0
	batch1 = 1:div(numwalkers, 2)
	batch2 = div(numwalkers, 2) + 1:numwalkers
	divisions = [(batch1, batch2), (batch2, batch1)]
    zs = zeros(Float64, numwalkers)
    newllhoods = zeros(Float64, numwalkers)
    proposals = zeros(Float64, (numwalkers,size(x0, 2)))
	rand_inactive = 0.0
    inv_a = 1.0 / a
	for i = 1:numsamples_perwalker
		for ensembles in divisions
			active, inactive = ensembles
            newllhoods = rand(rng, length(active))
            for ia in 1:length(active)
                @views zs[ia] = ((a - 1) * newllhoods[ia] + 1)^2 * inv_a
				rand_inactive = rand(rng, inactive)
                for j in 1:size(proposals, 2)
                    @views proposals[ia, j] = zs[ia] * x[active[ia], j] + (1 - zs[ia]) * x[rand_inactive, j]
                end
                @views newllhoods[ia] = llhood(proposals[ia, :])
            end
			for (j, walkernum) in enumerate(active)
				@views logratio = (size(x, 2) - 1) * log(zs[j]) + newllhoods[j] - lastllhoodvals[walkernum]
				if log(rand(rng)) < logratio
					lastllhoodvals[walkernum] = newllhoods[j]
					x[walkernum, :] = proposals[j, :]
				end
				if i % thinning == 0
					@views chain[div(i, thinning), walkernum, :] = x[walkernum, :]
				end
			end
		end
	end
	return chain
end

function flattenmcmcarray(chain::Array{Float64})
	numsteps, numwalkers, numdims = size(chain)
	newchain = Array{Float64}(undef, numdims, numwalkers * numsteps)
	@views for j = 1:numsteps
		for i = 1:numwalkers
			newchain[:, i + (j - 1) * numwalkers] .= chain[j, i, :]
		end
	end
	return newchain
end

const HDENSITIES = reverse(sort(HDENSITIES0))

const DlogN_dV = SPECIFIC_COLUMN_DENSITIES[2] - SPECIFIC_COLUMN_DENSITIES[1]
const DlognH = HDENSITIES[2] - HDENSITIES[1]
const DTg = GAS_TEMPERATURES[2] - GAS_TEMPERATURES[1]

const len_nH = length(HDENSITIES)
const len_Tg = length(GAS_TEMPERATURES)

function get_models_grid_frequency_values(filename)
    open(filename) do f
        line = readline(f)
        values_str = split(line, ",")[2]
        values = parse.(Float64, split(values_str))
        return values
    end
end
const MODEL_LINE_FREQS = get_models_grid_frequency_values(FILE_NAME_WITH_MODELS_GRID)
#const MODEL_LINE_FREQS = [220.4758078, 220.539324, 220.59442370000002, 220.64108439999998, 220.6792874, 220.709017, 220.7302611, 220.7430111, 220.7472617]

const modelTbs = open(joinpath(DATA_DIR, "modelTbs_G12.jls"), "r") do file
    deserialize(file)
end
const modelTaus = open(joinpath(DATA_DIR, "modelTaus_G12.jls"), "r") do file
    deserialize(file)
end

#const modelTbs = zeros(2, 2)
#const modelTaus = zeros(2, 2)

function one_to_two_dim(index)
    # Calculate row and column indices based on 1-based indexing
    row = ceil(Int, index / 3)
    col = mod(index - 1, 3) + 1
    return [row, col]
end

function two_to_one_dim(row, col)
    # Calculate the 1-based index from the given row and column
    index = (row - 1) * 3 + col
    return index
end

mutable struct LineProfileModel
    freq::Float64
    obs_X::Vector{Float64}
    obs_Y::Vector{Float64}
    p::Vector{Float64}
    err_p::Vector{Float64}
    function LineProfileModel()
        new(0.0, Float64[], Float64[], Float64[], Float64[])
    end
    function LineProfileModel(p, obs_X, obs_Y)
        new(0.0, obs_X, obs_Y, p, zeros(length(p)))
    end
end

mutable struct Blend
    obs_X::Vector{Float64}
    obs_Y::Vector{Float64}
    modelids::Vector{Int}
    function Blend()
        new(Float64[], Float64[], Int[])
    end
end

mutable struct Cell
    id::Int64
    logN_dV::Float64
    logN::Float64
    lognH::Float64
    Tg::Float64
    Nrot::Float64
    Nrot_err::Float64
    Trot::Float64
    Trot_err::Float64
    fill_fac::Float64
    fill_fac_err::Float64
    fill_fac_rot::Float64
    fill_fac_rot_err::Float64
    tau::Float64
    min_tau::Float64
    cloud_size::Float64
    err_logN_dV::Float64
    err_lognH::Float64
    err_Tg::Float64
    fwhm::Float64
    max_vlsr_diff::Float64
    x::Int64
    y::Int64
    rms::Float64
    vlsr::Float64
    err_vlsr::Float64
    freqs::Vector{Float64}
    obs_spec::Vector{Float64}
    obs_lines_freq
    obs_lines
    rhat_phys::Vector{Float64}
    ess_phys::Vector{Float64}
    gaussians::Vector{LineProfileModel}
    all_gaussians::Vector{LineProfileModel}
    num_of_lines_to_fit::Int

    function Cell(id::Int64, freqs::Vector{Float64}, obs_spec_in::Vector{Float64}, rms::Float64, x::Int64, y::Int64)
        first_non_nan_index = findfirst(!isnan, obs_spec_in)
        freqs = freqs[first_non_nan_index:end]
        obs_spec = obs_spec_in[first_non_nan_index:end]
        new(id, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, x, y, rms, VLSR_IN_KMS, NaN, freqs, obs_spec, Any, Any, Float64[], Float64[], LineProfileModel[], LineProfileModel[], 0)
    end
end

function separate_non_nan_ranges(x, arr)
    # Find where the NaN values are
    mask = isnan.(arr)
    # Find where the condition changes
    idx = findall(diff(mask) .!= 0)
    # Add the start and end indices
    idx = vcat(1, idx .+ 1, length(arr) + 1)
    # Separate the array into non-NaN ranges
    result_x = [x[idx[i]:idx[i+1]-1] for i in 1:length(idx)-1 if !mask[idx[i]]]
    result_arr = [arr[idx[i]:idx[i+1]-1] for i in 1:length(idx)-1 if !mask[idx[i]]]
    return result_x, result_arr
end

function gauss_profile(x::Vector{Float64}, p)
    p32 = - 1.0 / p[3]^2
    return p[1] .* exp.((x .- p[2]).^2 .* p32)
end

function gauss_profile(x::Float64, p1::Float64, p2::Float64, inv_p3::Float64)
    return p1 * exp(- ((x - p2) * inv_p3)^2)
end

function gauss_profile(x::Float64, p)
    return p[1] * exp(- ((x - p[2]) / p[3])^2)
end

function minimize_gauss_profile(x::Vector{Float64}, obs::Vector{Float64}, p0::Vector{Float64}, l::Vector{Float64}, u::Vector{Float64})
    function err_gauss(p)
        return sum((gauss_profile(x, p) .- obs) .^2)
    end
    function err_gauss_lim(p)
        if !all(l .<= p .<= u)
            return Inf
        end
        return sum((gauss_profile(x, p) .- obs).^2)
    end
    J = Array{Float64}(undef, length(x), length(p0))
    function err_gauss_profile_grad_lsqfit(x, p)
        inv_p1 = 1.0 / p[1]
        inv_p3 = 1.0 / p[3]
        inv_p32 = inv_p3^2
        model = 0.0
        diff_x = x[1] - p[2]
        for i in eachindex(x)
            model = gauss_profile(x[i], p[1], p[2], inv_p3)
            J[i, 1] = model * inv_p1
            J[i, 2] = 2 * model * diff_x * inv_p32
            J[i, 3] = J[i, 2] * diff_x * inv_p3
        end
        return J
    end
    function Avv!(dir_deriv, p, v)
        Hes = 0.0
        modelHes = 0.0
        diff_x = 0.0
        p32 = 0.0
        diff_x2 = 0.0
        for ix in eachindex(x)
            dir_deriv[ix] = 0.0
            modelHes = 2 * gauss_profile(x[ix], p)
            diff_x = x[ix] - p[2]
            diff_x2 = diff_x^2
            p32 = p[3]^2
            for i in 1:3
                for j in 1:3
                    Hes = 0.0
                    if (i == 1 && j == 2) || (i == 2 && j == 1) 
                        Hes = modelHes * diff_x / (p[1] * p32)
                    elseif (i == 1 && j == 3) || (i == 3 && j == 1)
                        Hes = modelHes * diff_x2 / (p[1] * p32 * p[3])
                    elseif i == 2 && j == 2
                        Hes = modelHes * ( 2 * diff_x2 - p32 ) / p32^2
                    elseif (i == 2 && j == 3) || (i == 3 && j == 2)
                        Hes = 2 * diff_x * modelHes * ( diff_x2 - p32 ) / (p32 * p[3]^3)
                    elseif i == 3 && j == 3
                        Hes = diff_x2 * modelHes * ( 2 * diff_x2 - 3 * p32 ) / (p32 * p[3]^4)
                    end
                    dir_deriv[ix] += Hes * v[i] * v[j]
                end
            end
        end
    end
    #p0 = Optim.minimizer(optimize(err_gauss, p0, SimulatedAnnealing(), Optim.Options(iterations = 1000000)))
    p0 = Optim.minimizer(optimize(err_gauss, p0, ParticleSwarm(lower = l, upper = u, n_particles = 10 * length(p0)), Optim.Options(iterations = 10000)))
    #p0 = Optim.minimizer(optimize(err_gauss_lim, p0, NelderMead(), Optim.Options(iterations = 10000, successive_f_tol=10000)))
    fit = curve_fit(gauss_profile, err_gauss_profile_grad_lsqfit, x, obs, p0, maxIter=1000000000, lower=l, upper=u, avv! = Avv!, min_step_quality=0)
    return LsqFit.coef(fit)
end
         
function minimize_gauss_blend(x::Vector{Float64}, obs::Vector{Float64}, num_of_peaks::Int, p0::Vector{Float64}, l::Vector{Float64}, u::Vector{Float64}, rms::Float64)
    model = zeros(length(x))
    function err_gauss_blend(p)
        peaks_are_ordered = true
        for i in 0:num_of_peaks-2
            @views peaks_are_ordered *= p[i*3+2] < p[(i+1)*3+2]
        end
        if !peaks_are_ordered
            return Inf
        end
        model .= 0.0
        for i in 1:num_of_peaks
            model .+= gauss_profile(x, @view(p[(i-1)*3+1:(i-1)*3+3]))
        end
        return sum((model .- obs).^2)
    end
    function err_gauss_blend_lim(p)
        peaks_are_ordered = true
        for i in 0:num_of_peaks-2
            @views peaks_are_ordered *= p[i*3+2] < p[(i+1)*3+2]
        end
        if !all(l .<= p .<= u) || !peaks_are_ordered
            return Inf
        end
        model .= 0.0
        for i in 1:num_of_peaks
            model .+= gauss_profile(x, @view(p[(i-1)*3+1:(i-1)*3+3]))
        end
        return sum((model .- obs).^2)
    end
    function model_lsqfit(x, p)
        model .= 0.0
        for i in 1:num_of_peaks
            model .+= gauss_profile(x, @view(p[(i-1)*3+1:(i-1)*3+3]))
        end
        return model
    end
    Jacobian = Array{Float64}(undef, length(x), length(p0))
    function model_lsqfit_grad(x, p_i)
        modelJ = 0.0
        inv_p1 = 1.0
        inv_p3 = 1.0
        diff_x = 0.0
        for i in 1:num_of_peaks
            @views temp_p = p_i[(i-1)*3+1:(i-1)*3+3]
            inv_p1 = 1.0 / temp_p[1]
            inv_p3 = 1.0 / temp_p[3]
            inv_p32 = inv_p3^2
            for ix in eachindex(x)
                modelJ = gauss_profile(x[ix], temp_p[1], temp_p[2], inv_p3)
                Jacobian[ix, (i - 1) * 3 + 1] = modelJ * inv_p1
                diff_x = x[ix] - temp_p[2]
                Jacobian[ix, (i - 1) * 3 + 2] = 2 * modelJ * diff_x * inv_p32
                Jacobian[ix, (i - 1) * 3 + 3] = Jacobian[ix, (i - 1) * 3 + 2] * diff_x * inv_p3
            end
        end
        return Jacobian
    end
    function Avv!(dir_deriv, p, v)
        Hes = 0.0
        modelHes = 0.0
        diff_x = 0.0
        p32 = 0.0
        diff_x2 = 0.0
        for ix in eachindex(x)
            dir_deriv[ix] = 0.0
            for ip in 1:num_of_peaks
                @views temp_p = p[(ip-1)*3+1:(ip-1)*3+3]
                modelHes = 2 * gauss_profile(x[ix], temp_p)
                diff_x = x[ix] - temp_p[2]
                diff_x2 = diff_x^2
                p32 = temp_p[3]^2
                for i in 1:3
                    for j in 1:3
                        Hes = 0.0
                        if (i == 1 && j == 2) || (i == 2 && j == 1) 
                            Hes = modelHes * diff_x / (temp_p[1] * p32)
                        elseif (i == 1 && j == 3) || (i == 3 && j == 1)
                            Hes = modelHes * diff_x2 / (temp_p[1] * p32 * temp_p[3])
                        elseif i == 2 && j == 2
                            Hes = modelHes * ( 2 * diff_x2 - p32 ) / p32^2
                        elseif (i == 2 && j == 3) || (i == 3 && j == 2)
                            Hes = 2 * diff_x * modelHes * ( diff_x2 - p32 ) / (p32 * temp_p[3]^3)
                        elseif i == 3 && j == 3
                            Hes = diff_x2 * modelHes * ( 2 * diff_x2 - 3 * p32 ) / (p32 * temp_p[3]^4)
                        end
                        dir_deriv[ix] += Hes * v[(ip-1)*3+i] * v[(ip-1)*3+j]
                    end
                end
            end
        end
    end
    function model_nlopt(p::Vector{Float64}, grad::Vector{Float64})
        model .= 0.0
        for i in 1:num_of_peaks
            model .+= gauss_profile(x, @view(p[(i-1)*3+1:(i-1)*3+3]))
        end
        return sum((model .- obs).^2)
    end
    function model_nlopt_constraints(p::Vector{Float64}, grad::Vector{Float64})
        peaks_are_ordered = true
        for i in 0:num_of_peaks-2
            peaks_are_ordered *= p[i*3+2] < p[(i+1)*3+2]
        end
        if !peaks_are_ordered || !all(l .<= p .<= u) !! any(isinf, p)
            return Inf
        end
        return -1.0
    end
    #p0 = Optim.minimizer(optimize(err_gauss_blend_lim, p0, SimulatedAnnealing(), Optim.Options(iterations = 1000000, successive_f_tol=10000)))
    p0 = Optim.minimizer(optimize(err_gauss_blend, p0, ParticleSwarm(lower = l, upper = u, n_particles = 10 * length(p0)), Optim.Options(iterations = 10000)))
    #p0 = Optim.minimizer(optimize(err_gauss_blend_lim, p0, SimulatedAnnealing(), Optim.Options(iterations = 1000000, successive_f_tol=10000)))
    #p0 = Optim.minimizer(optimize(err_gauss_blend_lim, p0, NelderMead(), Optim.Options(iterations = 10000, successive_f_tol=10000)))
    #res = compare_optimizers(err_gauss_blend; SearchRange = collect(zip(l, u)), NumDimensions = length(p0), MaxFuncEvals = 10000000)
    #res = bboptimize(err_gauss_blend, p0; SearchRange = collect(zip(l, u)), NumDimensions = length(p0), MaxFuncEvals = 1000000, TraceMode = :silent, PopulationSize = 100 * length(p0), Method=:adaptive_de_rand_1_bin_radiuslimited)
    #p0 = best_candidate(res)
    fit = curve_fit(model_lsqfit, model_lsqfit_grad, x, obs, p0, lower=l, upper=u, maxIter=1000000000, avv! = Avv!, min_step_quality=0)
    return LsqFit.coef(fit)
    #=opt = NLopt.Opt(:GN_ORIG_DIRECT, length(p0))
    NLopt.lower_bounds!(opt, l)
    NLopt.upper_bounds!(opt, u)
    #NLopt.xtol_rel!(opt, 1.e-6)
    NLopt.min_objective!(opt, model_nlopt)
    NLopt.inequality_constraint!(opt, model_nlopt_constraints, 1.e-8)
    min_f, min_x, ret = NLopt.optimize(opt, p0)
    println("ret: ", ret)
    id_p0 = 1
    for i in 1:num_of_peaks
        for j in 1:3
            if !(l[id_p0] <= min_x[id_p0] <= u[id_p0])
                println("#warning_new initial guess is out of bounds for $j : $(l[id_p0]) <= $(min_x[id_p0]) <= $(u[id_p0])")
                println(iblend.obs_X[1], " to ", iblend.obs_X[end])
                println(obs_X0[1], " to ", obs_X0[end])
                #println(obs_Y0)
                #println(rms)
            end
            id_p0 += 1
        end
    end
    opt = NLopt.Opt(:LN_COBYLA, length(p0))
    NLopt.lower_bounds!(opt, l)
    NLopt.upper_bounds!(opt, u)
    #NLopt.xtol_rel!(opt, 1.e-14)
    NLopt.min_objective!(opt, model_nlopt)
    NLopt.inequality_constraint!(opt, model_nlopt_constraints, 1.e-8)
    min_f, min_x, ret = NLopt.optimize(opt, min_x)
    fit = curve_fit(model_lsqfit, model_lsqfit_grad, x, obs, min_x, maxIter=1000000, lower=l, upper=u, x_tol=1.e-10, g_tol=1.e-12)
    #return min_x, LsqFit.stderror(fit)
    return LsqFit.coef(fit), LsqFit.stderror(fit)=#
end

function get_intensity_errors_using_mcmc(x::Vector{Float64}, obs::Vector{Float64}, num_of_peaks::Int, p0::Vector{Float64}, p0_err::Vector{Float64}, rms::Float64, l::Vector{Float64}, u::Vector{Float64})
    model = zeros(length(x))
    inv_rms2 = -0.5 / (rms * rms)
    #const_fac = - 0.5 * length(obs) * log(2.0 * pi * rms * rms)
    #inv_errs = 1.0 ./ p0_err
    function log_prob_blend(p)
        peaks_are_ordered = true
        for i in 0:num_of_peaks-2
            @views peaks_are_ordered *= p[i*3+2] < p[(i+1)*3+2]
        end
        if !all(l .<= p .<= u) || !peaks_are_ordered
            return -Inf
        end
        model .= 0.0
        for i in 1:num_of_peaks
            model .+= gauss_profile(x, @view(p[(i-1)*3+1:(i-1)*3+3]))
        end
        return sum((model .- obs).^2) * inv_rms2 #- 0.5 * sum(((p .- p0).* inv_errs).^2) #+ const_fac
    end
    local_rng = Random.TaskLocalRNG()
    Random.seed!(local_rng, 11051989)
    numwalkers = 2 * length(p0) + 2
    burnin = 1000
    numsamples_perwalker = 100000
    p00 = zeros(numwalkers, length(p0))
    p00[1, :] .= copy(p0)
    p_cand = zeros(length(p0))
    normal_dist = Normal(0.0, 1.0)
    res_modified = copy(p0)
    for iw in 2:numwalkers
        while true
            p_cand .= res_modified .* (1.0 .+ 1.e-7 .* rand(local_rng, normal_dist, length(p_cand)))
            if isfinite(log_prob_blend(@view(p_cand[:])))
                p00[iw, :] .= p_cand
                break
            end
        end
    end
    chain = sample(log_prob_blend, numwalkers, p00, burnin, 1, rng=local_rng)
    chain = sample(log_prob_blend, numwalkers, chain[end, :, :], numsamples_perwalker, 1, rng=local_rng)
    flatchain = flattenmcmcarray(chain)
    res_modified_l = copy(p0)
    res_modified_u = copy(p0)
    for ipar in eachindex(p0)
        res_modified[ipar] = median(@view(flatchain[ipar, :]))
        res_modified_l[ipar] = percentile(@view(flatchain[ipar, :]), 16)
        res_modified_u[ipar] = percentile(@view(flatchain[ipar, :]), 84)
    end
    #println("gauss rhat: ", rhat(chain))
    #fig_gauss = pairplot(flatchain')
    #save(joinpath(DATA_DIR, "figures/corner_plot.png"), fig_gauss)
    return res_modified, (res_modified_u .- res_modified_l) .* 0.5
end

function get_K_from_line_label(line_label)
    res = split(line_label, "->")[1]
    res = split(res, "_")[2]
    res = replace(res, "{" => "")
    res = replace(res, "}" => "")
    return res
end

function check_point_is_gaussian(left_X, right_X, left_Y, right_Y, peak_Y, peak_X, min_width, max_width, init_width, rms)
    temp_model_obs_X = vcat(left_X, right_X)
    if length(temp_model_obs_X) < 4
        return true
    end
    deltaFreq = 0.5 * abs(temp_model_obs_X[2] - temp_model_obs_X[1])
    temp_model_obs_Y = vcat(left_Y, right_Y)
    sorting_ids = sortperm(temp_model_obs_X)
    temp_model_obs_X = temp_model_obs_X[sorting_ids]
    temp_model_obs_Y = temp_model_obs_Y[sorting_ids]
    p0 = [peak_Y, peak_X, init_width]
    p = minimize_gauss_profile(temp_model_obs_X, temp_model_obs_Y, p0, [peak_Y - rms, peak_X - deltaFreq, min_width], [peak_Y + rms, peak_X + deltaFreq, max_width])
    difference = temp_model_obs_Y .- gauss_profile(temp_model_obs_X, p)
    max_diff = maximum(abs.(difference))
    if std(difference) < rms && max_diff < rms
        return true
    else
        return false
    end
end

# CLEAN like algorithm to fit a spectrum with Gaussians
function fit_spectrum_with_gaussians(obs_X0, obs_Y0, rms, threshold, fig, fig_row, fig_col)
    obs_X = copy(obs_X0)
    obs_Y = copy(obs_Y0)
    deltaFreq = 0.5 * abs(obs_X[end] - obs_X[1])
    min_doppler_b = MIN_LINE_FWHM_IN_KMS / (2. * sqrt(log(2.))) / 3.e5 * abs(obs_X[end] + obs_X[1]) * 0.5
    max_doppler_b = MAX_LINE_FWHM_IN_KMS / (2. * sqrt(log(2.))) / 3.e5 * abs(obs_X[end] + obs_X[1]) * 0.5
    models = Vector{LineProfileModel}()
    index_of_maximum_flux = argmax(abs.(obs_Y))
    need_to_search_lines = true
    need_to_change_sign = false
    if obs_Y[index_of_maximum_flux] < 0
        obs_Y = -copy(obs_Y)
        need_to_change_sign = true
    end
    if obs_Y[index_of_maximum_flux] < threshold
        need_to_search_lines = false
    end
    # Fit the spectrum peak with Gaussian, exclude the fitted Gaussian from spectrum and continue while remaining peak flux is above threshold
    while need_to_search_lines
        temp_model_obs_X = Vector{Float64}()
        temp_model_obs_Y = Vector{Float64}()
        bad_peak = false
        if index_of_maximum_flux > length(obs_X) - 2 || index_of_maximum_flux < 3
            bad_peak = true
        else
            append!(temp_model_obs_X, obs_X[index_of_maximum_flux])
            append!(temp_model_obs_Y, obs_Y[index_of_maximum_flux])
            left_pointsX = Vector{Float64}()
            left_pointsY = Vector{Float64}()
            right_pointsX = Vector{Float64}()
            right_pointsY = Vector{Float64}()
            left_point_id = index_of_maximum_flux - 1
            append!(left_pointsX, obs_X[left_point_id])
            append!(left_pointsY, obs_Y[left_point_id])
            left_point_id -= 1
            right_point_id = index_of_maximum_flux + 1
            append!(right_pointsX, obs_X[right_point_id])
            append!(right_pointsY, obs_Y[right_point_id])
            right_point_id += 1
            num_of_added_points = 3
            points_added = true
            while points_added
                add_left_point = false
                if obs_Y[left_point_id] <= left_pointsY[end] && left_point_id > 1 && obs_Y[left_point_id] >= 0.0 && check_point_is_gaussian(vcat(left_pointsX, obs_X[left_point_id]), right_pointsX, vcat(left_pointsY, obs_Y[left_point_id]), right_pointsY, obs_Y[index_of_maximum_flux], obs_X[index_of_maximum_flux], min_doppler_b, max_doppler_b, min_doppler_b, rms)
                    add_left_point = true
                end
                add_right_point = false
                if obs_Y[right_point_id] <= right_pointsY[end] && right_point_id < length(obs_X) && obs_Y[right_point_id] >= 0.0 && check_point_is_gaussian(left_pointsX, vcat(right_pointsX, obs_X[right_point_id]), left_pointsY, vcat(right_pointsX, obs_X[right_point_id]), obs_Y[index_of_maximum_flux], obs_X[index_of_maximum_flux], min_doppler_b, max_doppler_b, min_doppler_b, rms)
                    add_right_point = true
                end
                points_added = false
                if add_left_point
                    append!(left_pointsX, obs_X[left_point_id])
                    append!(left_pointsY, obs_Y[left_point_id])
                    left_point_id -= 1
                    num_of_added_points += 1
                    points_added = true
                end
                if add_right_point
                    append!(right_pointsX, obs_X[right_point_id])
                    append!(right_pointsY, obs_Y[right_point_id])
                    right_point_id += 1
                    num_of_added_points += 1
                    points_added = true
                end
            end
            if num_of_added_points < 4
                bad_peak = true
            else
                temp_model_obs_X = vcat(temp_model_obs_X, left_pointsX)
                temp_model_obs_Y = vcat(temp_model_obs_Y, left_pointsY)
                temp_model_obs_X = vcat(temp_model_obs_X, right_pointsX)
                temp_model_obs_Y = vcat(temp_model_obs_Y, right_pointsY)
            end
        end
        if bad_peak
            obs_Y = deleteat!(obs_Y, index_of_maximum_flux)
            obs_X = deleteat!(obs_X, index_of_maximum_flux)
            if length(obs_Y) == 0
                need_to_search_lines = false
            else
                index_of_maximum_flux = argmax(obs_Y)
                if obs_Y[index_of_maximum_flux] < threshold
                    need_to_search_lines = false
                end
            end
            continue
        end
        sorting_ids = sortperm(temp_model_obs_X)
        temp_model_obs_X = temp_model_obs_X[sorting_ids]
        temp_model_obs_Y = temp_model_obs_Y[sorting_ids]
        p0 = [obs_Y[index_of_maximum_flux], obs_X[index_of_maximum_flux], min_doppler_b]
        p = minimize_gauss_profile(temp_model_obs_X, temp_model_obs_Y, p0, [obs_Y[index_of_maximum_flux] - rms, obs_X[index_of_maximum_flux] - deltaFreq, min_doppler_b], [obs_Y[index_of_maximum_flux] + rms, obs_X[index_of_maximum_flux] + deltaFreq, max_doppler_b])
        if p[1] >= threshold
            push!(models, LineProfileModel(p, temp_model_obs_X, temp_model_obs_Y))
            obs_Y .= obs_Y .- gauss_profile(obs_X, p)
            index_of_maximum_flux = argmax(obs_Y)
        end
        if obs_Y[index_of_maximum_flux] < threshold || p[1] < threshold
            need_to_search_lines = false
        end
    end
    #display(Plots.plot!(collect(range(obs_X[1], stop = obs_X[end], length = 10000)), gauss_profile(collect(range(obs_X[1], stop = obs_X[end], length = 10000)), p), label = "model"))

    blends = Vector{Blend}()
    if length(models) > 0
        models = sort!(models, by = x -> x.p[2])
        if length(models) > 1
            # Filter out the Gaussians that are blended too much
            models_to_remove = []
            for i in 1:length(models)-1
                if models[i].p[2] + 0.1 * models[i].p[3] > models[i+1].p[2] - 0.1 * models[i+1].p[3]
                    if models[i].p[1] >= models[i+1].p[1]
                        push!(models_to_remove, i+1)
                    else
                        push!(models_to_remove, i)
                    end
                end
            end
            if length(models_to_remove) > 0
                models_to_remove = sort!(unique!(models_to_remove))
                models = deleteat!(models, models_to_remove)
            end
        end
        push!(blends, Blend())
        push!(blends[1].modelids, 1)
    end
    i_blend = 1

    # Unify somewhat the line widths
    mean_line_width = mean([x.p[3] for x in models])
    mean_line_width_err = std([x.p[3] for x in models])
    min_doppler_b_old = min_doppler_b
    max_doppler_b_old = max_doppler_b
    if mean_line_width_err > 1.e-14
        min_doppler_b = max(min_doppler_b, mean_line_width - mean_line_width_err)
        max_doppler_b = min(max_doppler_b, mean_line_width + mean_line_width_err)
    end
    
    # Combine Gaussians into blends. The Gaussian parameters within the blend will be fitted simultaneously. The Gaussian fit for each blend will be performed independently
    for i in 2:length(models)
        if models[i-1].p[2] + 3.0 / sqrt(2.) * models[i-1].p[3] > models[i].p[2] - 3.0 / sqrt(2.) * models[i].p[3]
            push!(blends[i_blend].modelids, i)
        else
            i_blend += 1
            push!(blends, Blend())
            push!(blends[i_blend].modelids, i)
        end
    end

    # Fit each blend
    all_models = []
    for iblend in blends
        # Get the observed spectrum points that correspond to a given blend
        found_enough_points = false
        start_X = 1
        end_X = length(obs_X0)
        left_model_id = 1
        right_model_id = length(iblend.modelids)
        for i in 2:length(iblend.modelids)-1
            if models[iblend.modelids[i]].p[2] - models[iblend.modelids[i]].p[3] < models[iblend.modelids[left_model_id]].p[2] - models[iblend.modelids[left_model_id]].p[3]
                left_model_id = i
            end
            if models[iblend.modelids[i]].p[2] + models[iblend.modelids[i]].p[3] > models[iblend.modelids[right_model_id]].p[2] + models[iblend.modelids[right_model_id]].p[3]
                right_model_id = i
            end
        end
        for widths in 2.0:0.1:5.0
            for i in 2:length(obs_X0)
                if models[iblend.modelids[left_model_id]].p[2] - widths * models[iblend.modelids[left_model_id]].p[3] / sqrt(2.) <= obs_X0[i]
                    start_X = i - 1
                    break
                end
            end
            for i in length(obs_X0)-1:-1:1
                if models[iblend.modelids[right_model_id]].p[2] + widths * models[iblend.modelids[right_model_id]].p[3] / sqrt(2.) >= obs_X0[i]
                    end_X = i + 1
                    break
                end
            end
            start_X = max(1, start_X)
            end_X = min(end_X, length(obs_X0))
            if end_X - start_X >= 3 * length(iblend.modelids)
                found_enough_points = true
                break
            else
                start_X = max(1, start_X - 1)
                end_X = min(end_X + 1, length(obs_X0))
                if end_X - start_X >= 3 * length(iblend.modelids)
                    found_enough_points = true
                    break
                end
            end
        end
        max_num_of_models = div(end_X - start_X + 1, 4)
        if max_num_of_models == 0
            println("#warning cannot find enough points for freq range $(obs_X0[1]) to $(obs_X0[end])")
            max_num_of_models = 1
        end
        iblend.obs_X = copy(obs_X0[start_X:end_X])
        iblend.obs_Y = copy(obs_Y0[start_X:end_X])
        if need_to_change_sign
            iblend.obs_Y = - iblend.obs_Y
        end
        res_models = sort!(models[iblend.modelids], by = x -> abs(x.p[1]), rev=true)
        max_num_of_models = min(max_num_of_models, length(res_models))
        res_models = res_models[1:max_num_of_models]
        res_models = sort!(res_models, by = x -> x.p[2])
        res_models_copy = deepcopy(res_models)
        max_obs_Y = maximum(iblend.obs_Y)
        # Fit of the observed spectrum points with Gaussians within a given blend. To minimize the number of Gaussians and avoid degenerate solutions, one removes Gaussians from the fit that are too much blended or too weak
        models_were_removed = true
        while models_were_removed
            p0 = zeros(length(res_models) * 3)
            id_p0 = 1
            l = zeros(length(res_models) * 3)
            u = zeros(length(res_models) * 3)
            for i in eachindex(res_models)
                for j in 1:3
                    p0[id_p0] = res_models_copy[i].p[j]
                    if j == 1
                        l[id_p0] = max(2 * rms, res_models[i].p[j] - rms)
                        u[id_p0] = min(max_obs_Y + rms, res_models[i].p[j] + rms)
                    end
                    if j == 2
                        l[id_p0] = res_models_copy[i].p[j] - 3 * res_models_copy[i].p[3]
                        u[id_p0] = res_models_copy[i].p[j] + 3 * res_models_copy[i].p[3]
                    end
                    if j == 3
                        p0[id_p0] = mean_line_width
                        l[id_p0] = min_doppler_b
                        u[id_p0] = max_doppler_b
                    end
                    if p0[id_p0] <= l[id_p0]
                        p0[id_p0] = l[id_p0] * 1.0001
                    end
                    if p0[id_p0] >= u[id_p0]
                        p0[id_p0] = u[id_p0] * 0.9999
                    end
                    id_p0 += 1
                end
            end
            result = minimize_gauss_blend(iblend.obs_X, iblend.obs_Y, length(res_models), p0, l, u, rms)
            for i in eachindex(res_models)
                for j in 1:3
                    res_models[i].p[j] = result[3 * (i-1) + j]
                end
            end
            models_were_removed = false
            models_to_remove = []
            res_models = sort!(res_models, by = x -> x.p[2])
            for i in 1:length(res_models)-1
                if res_models[i].p[2] + 0.5 * res_models[i].p[3] >= res_models[i+1].p[2] - 0.5 * res_models[i+1].p[3]
                    if abs(res_models[i].p[1]) >= abs(res_models[i+1].p[1])
                        push!(models_to_remove, i+1)
                    else
                        push!(models_to_remove, i)
                    end
                end
                if abs(res_models[i].p[1]) < 3 * rms
                    push!(models_to_remove, i)
                end
                if res_models[i].p[2] + 2 * res_models[i].p[3] < iblend.obs_X[1] || res_models[i].p[2] - 2 * res_models[i].p[3] > iblend.obs_X[end]
                    push!(models_to_remove, i)
                end
            end
            if abs(res_models[end].p[1]) < 3 * rms || res_models[end].p[2] + 2 * res_models[end].p[3] < iblend.obs_X[1] || res_models[end].p[2] - 2 * res_models[end].p[3] > iblend.obs_X[end]
                push!(models_to_remove, length(res_models))
            end
            if length(models_to_remove) > 0
                models_to_remove = sort!(unique!(models_to_remove))
                res_models = deleteat!(res_models, models_to_remove)
                res_models_copy = deleteat!(res_models_copy, models_to_remove)
                if length(res_models) > 0 
                    models_were_removed = true
                end
            end
        end
        # Get the fit errors using Monte-Carlo Markov Chain method
        result = zeros(3 * length(res_models))
        result_err = zeros(3 * length(res_models))
        l = zeros(3 * length(res_models))
        u = zeros(3 * length(res_models))
        for i in eachindex(res_models)
            for j in 1:3
                result[3 * (i-1) + j] = res_models[i].p[j]
                if j == 1
                    result_err[3 * (i-1) + j] = rms
                    l[3 * (i-1) + j] = 1.8 * rms
                    u[3 * (i-1) + j] =  max_obs_Y + 5 * rms
                end
                if j == 2
                    result_err[3 * (i-1) + j] = 0.1 * min_doppler_b
                    l[3 * (i-1) + j] = max(iblend.obs_X[1] - 5 * max_doppler_b, res_models[i].p[j] - 5 * res_models[i].p[3])
                    u[3 * (i-1) + j] = min(iblend.obs_X[end] + 5 * max_doppler_b, res_models[i].p[j] + 5 * res_models[i].p[3])
                end
                if j == 3
                    result_err[3 * (i-1) + j] = 0.1 * min_doppler_b
                    l[3 * (i-1) + j] = 0.1 * min_doppler_b_old
                    u[3 * (i-1) + j] = 10.0 * max_doppler_b_old
                end
            end
        end
        _, result_err = get_intensity_errors_using_mcmc(iblend.obs_X, iblend.obs_Y, length(res_models), result, result_err, rms, l, u)
        for i in eachindex(res_models)
            for j in 1:3
                res_models[i].err_p[j] = result_err[3 * (i-1) + j]
            end
        end
        for i in eachindex(res_models)
            if need_to_change_sign
                res_models[i].p[1] = - res_models[i].p[1]
            end
        end
        all_models = vcat(all_models, res_models)
    end

    axs = Axis(fig[fig_row, fig_col], xgridvisible = false, ygridvisible = false)

    if need_to_change_sign
        CairoMakie.hlines!(axs, [-threshold], color = :black, linestyle = :dash)
    else
        CairoMakie.hlines!(axs, [threshold], color = :black, linestyle = :dash)
    end

    CairoMakie.scatter!(axs, obs_X0, obs_Y0, color=:black, label="observations", markersize = 14)

    if length(all_models) > 0
        models = sort!(all_models, by = x -> x.p[2])
        xrange = collect(range(obs_X0[1], stop = obs_X0[end], length = 10000))
        difference = copy(obs_Y0)
        total_model = zeros(length(xrange))
        for model in models
            total_model .+= gauss_profile(xrange, model.p)
            difference .-= gauss_profile(obs_X0, model.p)
        end

        CairoMakie.lines!(axs, xrange, total_model, color = :green, label = "total model")
        CairoMakie.lines!(axs, obs_X0, difference, color = :blue, label = "difference")

        for (i, model) in enumerate(models)
            if i == 1
                CairoMakie.lines!(axs, xrange, gauss_profile(xrange, model.p), color = :red, label = (i == 1 ? "single model" : ""))
            else
                CairoMakie.lines!(axs, xrange, gauss_profile(xrange, model.p), color = :red)
            end
            #CairoMakie.errorbars!(axs, [model.p[2]], [model.p[1]], [model.err_p[1]]; color = :red)
        end
    end
    #ymax = axs.finallimits[].origin[2] + axs.finallimits[].widths[2]
    ymax = maximum(obs_Y0)
    for key in keys(LINE_VERTICAL_LABELS)
        if obs_X0[1] <= LINE_VERTICAL_LABELS[key] * (1 - VLSR_IN_KMS / 3.e5) <= obs_X0[end]
            CairoMakie.vlines!(axs, [LINE_VERTICAL_LABELS[key] * (1 - VLSR_IN_KMS / 3.e5)], color = :black, linestyle = :dot)
            line_K = get_K_from_line_label(key)
            CairoMakie.text!(axs, LINE_VERTICAL_LABELS[key] * (1 - VLSR_IN_KMS / 3.e5), ymax * 1.1, text = latexstring("\$K = $(line_K)\$"), fontsize = 14, rotation = 0, space = :data)
        end
    end
    return models, axs
end

# From all Gaussian components, find those that have a maximimum flux and are closest to the line frequency taking into account the Doppler shift
function get_needed_gaussians_and_vlsr(models)
    needed_gaussians = []
    individual_vlsr = []
    individual_vlsr_err = []
    line_label_id = 0
    for (line_label, freq0) in LINE_VERTICAL_LABELS
        line_label_id += 1
        if line_label_id in INDEXES_OF_LINES_THAT_WILL_BE_SELECTED_BASED_ON_VLSR_BUT_NOT_ON_INTENSITY
            continue
        end
        closest_models_ids = []
        for model_i in eachindex(models)
            min_freq_diff = Inf
            closest_line_label = "none"
            for (line_label1, freq1) in LINE_VERTICAL_LABELS
                freq_shfted1 = freq1 * (1 - VLSR_IN_KMS / 3.e5)
                freq_diff = abs(models[model_i].p[2] - freq_shfted1)
                if freq_diff < min_freq_diff
                    min_freq_diff = freq_diff
                    closest_line_label = line_label1
                end
            end
            if closest_line_label == line_label
                push!(closest_models_ids, model_i)
            end
        end
        if length(closest_models_ids) == 0
            continue
        end
        model_id_with_max_peak = argmax([abs(x.p[1]) for x in models[closest_models_ids]])
        models[closest_models_ids[model_id_with_max_peak]].freq = freq0
        line_vlsr = (freq0 / models[closest_models_ids[model_id_with_max_peak]].p[2] - 1) * 3e5
        if MIN_ALLOWED_LINE_VEL_WRT_VLSR_KMS < line_vlsr - VLSR_IN_KMS < MAX_ALLOWED_LINE_VEL_WRT_VLSR_KMS
            needed_gaussians = vcat(needed_gaussians, models[closest_models_ids[model_id_with_max_peak]])
            push!(individual_vlsr, line_vlsr)
            push!(individual_vlsr_err, models[closest_models_ids[model_id_with_max_peak]].err_p[2] / freq0 * 3e5)
        end
    end
    mean_vlsr = VLSR_IN_KMS
    vlsr_err = NaN
    vlsr_max_diff = NaN
    if length(individual_vlsr) > 0
        mean_vlsr = mean(individual_vlsr)
        vlsr_max_diff = abs(maximum(individual_vlsr) - minimum(individual_vlsr))
        vlsr_err = sqrt(sum(individual_vlsr_err .^ 2)) / length(individual_vlsr)
    end
    
    individual_vlsr_01 = []
    individual_vlsr_err_01 = []
    line_label_id = 0
    for (line_label, freq0) in LINE_VERTICAL_LABELS
        line_label_id += 1
        if !(line_label_id in INDEXES_OF_LINES_THAT_WILL_BE_SELECTED_BASED_ON_VLSR_BUT_NOT_ON_INTENSITY)
            continue
        end
        closest_models_ids = []
        individual_vlsr = []
        for model_i in eachindex(models)
            min_freq_diff = Inf
            closest_line_label = "none"
            for (line_label1, freq1) in LINE_VERTICAL_LABELS
                freq_shfted1 = freq1 * (1 - mean_vlsr / 3.e5)
                freq_diff = abs(models[model_i].p[2] - freq_shfted1)
                if freq_diff < min_freq_diff
                    min_freq_diff = freq_diff
                    closest_line_label = line_label1
                end
            end
            if closest_line_label == line_label
                push!(closest_models_ids, model_i)
            end
        end
        if length(closest_models_ids) == 0
            continue
        end
        model_id_with_closest_peak = argmin([abs(x.p[2] - freq0 * (1 - mean_vlsr / 3.e5)) for x in models[closest_models_ids]])
        models[closest_models_ids[model_id_with_closest_peak]].freq = freq0
        line_vlsr = (freq0 / models[closest_models_ids[model_id_with_closest_peak]].p[2] - 1) * 3e5
        if MIN_ALLOWED_LINE_VEL_WRT_VLSR_KMS < line_vlsr - VLSR_IN_KMS < MAX_ALLOWED_LINE_VEL_WRT_VLSR_KMS
            needed_gaussians = vcat(needed_gaussians, models[closest_models_ids[model_id_with_closest_peak]])
            push!(individual_vlsr_01, line_vlsr)
            push!(individual_vlsr_err_01, models[closest_models_ids[model_id_with_closest_peak]].err_p[2] / freq0 * 3e5)
        end
    end
    if isnan(vlsr_err) && length(individual_vlsr_01) > 0
        mean_vlsr = mean(individual_vlsr_01)
        vlsr_max_diff = abs(maximum(individual_vlsr_01) - minimum(individual_vlsr_01))
        vlsr_err = sqrt(sum(individual_vlsr_err_01 .^ 2)) / length(individual_vlsr_01)
    end
    if isnan(vlsr_err)
        mean_vlsr = NaN
    end
    if length(needed_gaussians) > 0
        needed_gaussians = sort!(needed_gaussians, by = x -> x.p[2])
    end

    return needed_gaussians, mean_vlsr, vlsr_err, vlsr_max_diff
end

function get_gaussians(obs_lines_freq, obs_lines, rms, cell_id)
    figs = CairoMakie.Figure(size = (1500, 1200), fontsize=18)
    fig = figs[1, 1] = CairoMakie.GridLayout()
    gaussians = []
    all_gaussians = []
    plot_id = 1
    max_row = 1
    axes_arr = []
    for (xaxis, line) in zip(obs_lines_freq, obs_lines)
        fig_row, fig_col = one_to_two_dim(plot_id)
        line_models, axs = fit_spectrum_with_gaussians(xaxis, line, rms, THRESHOLD_IN_SIGMAS * rms, fig, fig_row, fig_col)
        if length(line_models) != 0
            push!(all_gaussians, line_models)
        end
        push!(axes_arr, axs)
        max_row = max(max_row, fig_row)
        plot_id += 1
    end
    all_gaussians = vcat(all_gaussians...)
    if length(all_gaussians) > 0
        all_gaussians =sort!(all_gaussians, by = x -> x.p[2])
    end
    gaussians, mean_vlsr, err_vlsr, max_vlsr_diff = Vector{LineProfileModel}(), NaN, NaN, NaN
    if length(all_gaussians) > 0
        gaussians, mean_vlsr, err_vlsr, max_vlsr_diff = get_needed_gaussians_and_vlsr(all_gaussians)
    else
        println("#warning no lines found in cell $cell_id")
    end
    fig_row, fig_col = one_to_two_dim(length(axes_arr)+1)
    lock(plotting_lock) do
        axes_arr[two_to_one_dim(max_row, 1)].xlabel = "Frequency, GHz"
        axes_arr[two_to_one_dim(max_row, 1)].ylabel = "Intensity, K"
        #CairoMakie.Legend(fig[fig_row, fig_col], axes_arr[end])
        axes_arr[1].title = "Cell id = $cell_id"
    end
    CairoMakie.save(joinpath(DATA_DIR, "figures/spectra/$cell_id.png"), figs)
    return gaussians, all_gaussians, mean_vlsr, err_vlsr, max_vlsr_diff
end

# Fit rotational diagram
function perform_rotational_analysis(freq::Vector{Float64}, Eu::Vector{Float64}, Aul::Vector{Float64}, g::Vector{Float64}, W::Vector{Float64}, W_err::Vector{Float64}, deltav::Vector{Float64}, cell_id::Int64, logN0::Float64, logN0_err::Float64, temp0::Float64, temp0_err::Float64, omega_in::Float64 = 1.0, omega_in_err::Float64 = Inf)

    function line_flux_optthick(omega, ntot, q, temp, nu, gup, tup, aup, deltav)
        Nup = gup * ntot / q * exp(-tup / temp)
        tau = aup * (c_cgs / nu)^3. / (8 * pi * deltav) * Nup * expm1(h_cgs * nu / (k_cgs * temp))
        lineflux_dV = omega * Nup * aup * (c_cgs / nu)^2. * (h_cgs * c_cgs / (8 * pi * k_cgs))
        if tau > 0.0
            lineflux_dV *= (-expm1(-tau) / tau)
        end
        return lineflux_dV * 1.e-5 # cm/s -> km/s
    end

    inv_sigma_obs = 1 ./ W_err
    ln2pisigma = log.((2 * pi) .* W_err.^2)
    logNtotlimit1 = SPECIFIC_COLUMN_DENSITIES[1] + log10(MIN_LINE_FWHM_IN_KMS * 1.e5)
    logNtotlimit2 = SPECIFIC_COLUMN_DENSITIES[end] + log10(MAX_LINE_FWHM_IN_KMS * 1.e5)
    omegalimit1 = 1.e-8
    omegalimit2 = 1.0

    function rotloglikelihood(theta)
        logntot, temp, omega = theta
        @views if !(logNtotlimit1 <= logntot <= logNtotlimit2 && GAS_TEMPERATURES[1] <= temp <= GAS_TEMPERATURES[end] && omegalimit1 <= omega <= omegalimit2)
            return -Inf
        end
        ntot = 10.0 ^ logntot
        q = get_partition_function(temp)
        likelihoodsum = 0.0
        for iline in eachindex(freq)
            @views likelihoodsum += ((line_flux_optthick(omega, ntot, q, temp, freq[iline], g[iline], Eu[iline], Aul[iline], deltav[iline]) - W[iline]) * inv_sigma_obs[iline])^2 + ln2pisigma[iline]
        end
        likelihoodsum += ((logntot - logN0) / logN0_err)^2 + ((temp - temp0) / temp0_err)^2 + ((omega - omega_in) / omega_in_err)^2
        return -0.5 * likelihoodsum
    end

    function get_model_and_obs_lnNu_gu(logntot, temp, omega)
        ntot = 10.0^logntot
        q = get_partition_function(temp)
        Nup_g = zeros(length(freq))
        Nup_obs = zeros(length(freq))
        N_err = zeros(length(freq))
        for i in eachindex(freq)
            Nup_g[i] = ntot / q * exp(-Eu[i] / temp)
            tau = Aul[i] * (c_cgs / freq[i])^3. / (8 * pi * deltav[i]) * Nup_g[i] * g[i] * expm1(h_cgs * freq[i] / (k_cgs * temp))
            W_cgs = W[i] * 1e5
            Nup_obs[i] = - ((8 * pi * k_cgs * freq[i]^2 * W_cgs) / (h_cgs * c_cgs^3 * Aul[i])) * tau / expm1(-tau) / g[i] / omega
            N_err[i] = W_err[i] / W[i]
        end
        return log.(Nup_g), log.(Nup_obs), N_err
    end
    
    local_rng = Random.TaskLocalRNG()
    Random.seed!(local_rng, 11051989)
    numwalkers = 16
    logntotini = rand(local_rng, Uniform(logNtotlimit1, logNtotlimit2), numwalkers)
    tempini = rand(local_rng, Uniform(GAS_TEMPERATURES[1], GAS_TEMPERATURES[end]), numwalkers)
    omegaini = rand(local_rng, Uniform(omegalimit1, omegalimit2), numwalkers)
    p00 = hcat(logntotini, tempini, omegaini)
    chain = sample(rotloglikelihood, numwalkers, p00, 100, 1, rng=local_rng)
    chain = sample(rotloglikelihood, numwalkers, chain[end, :, :], 10000, 1, rng=local_rng)
    flatchain = flattenmcmcarray(chain)
    res_modified = zeros(size(p00, 2))
    res_modified_l = zeros(size(p00, 2))
    res_modified_u = zeros(size(p00, 2))
    for ipar in eachindex(res_modified)
        res_modified[ipar] = median(@view(flatchain[ipar, :]))
        res_modified_l[ipar] = percentile(@view(flatchain[ipar, :]), 16)
        res_modified_u[ipar] = percentile(@view(flatchain[ipar, :]), 84)
    end

    lgNu_mod, ln_Nu_gu, ln_Nu_gu_err = get_model_and_obs_lnNu_gu(res_modified[1], res_modified[2], res_modified[3])

    fig = CairoMakie.Figure(size = (1500, 1200), fontsize=18)
    ax = CairoMakie.Axis(fig[1, 1])
    CairoMakie.lines!(ax, Eu, lgNu_mod; color = :green)
    CairoMakie.errorbars!(ax, Eu, ln_Nu_gu, ln_Nu_gu_err; color = :red)
    ax.xlabel = "Upper level energy, K"
    ax.ylabel = "ln(Nu/g)"
    ax.title = "Cell id = $(cell_id)"
    CairoMakie.save(joinpath(DATA_DIR, "figures/rot_diagram/$(cell_id).png"), fig)
    return res_modified[1], (res_modified_u[1] - res_modified_l[1]) * 0.5, res_modified[2], (res_modified_u[2] - res_modified_l[2]) * 0.5, res_modified[3], (res_modified_u[3] - res_modified_l[3]) * 0.5
end

function get_phys_model(logN_dV0, lognH0, Tg0, line_list)
    logN_dV = logN_dV0 - SPECIFIC_COLUMN_DENSITIES[1]
    index_logN_dV = round(Int, logN_dV / DlogN_dV)
    lognH = lognH0 - HDENSITIES[1]
    index_lognH = round(Int, lognH / DlognH)
    Tg = Tg0 - GAS_TEMPERATURES[1]
    index_Tg = round(Int, Tg / DTg)
    model_index = 1 + index_logN_dV * len_nH * len_Tg + index_lognH * len_Tg + index_Tg
    if model_ids[model_index] == 0
        return NaN, NaN
    end
    return modelTbs[model_ids[model_index], line_list], modelTaus[model_ids[model_index], MODEL_FREQ_INDEX_FOR_OPTICAL_DEPTH_OUTPUT], minimum(modelTaus[model_ids[model_index], line_list])
end

#=const threshold_density_factor = 3.5 * 4.0 / 3.0 * (0.1e-4*0.5)^3 *1.e2 / (pi*(0.1e-4*0.5)^2) / 2.8 / 1.67262192e-24 * sqrt(1.5e15*4680)
function threshold_density(Tg0) # see https://academic.oup.com/mnras/article/482/4/4673/5173064 https://doi.org/10.1093/mnras/sty3038
    Tg = 1.0 / Tg0
    return log10(threshold_density_factor*sqrt(Tg)*exp(-4680*Tg))
end
=#

# Get the radiative transfer model that is closest to the observed spectrum and satisfies the constraints on the flux and physical parameters
function get_closest_phys_model(obs, inv_sigma_obs, line_list, absent_lines_list, threshold, max_Tgas, fwhm, fill_fac)
    function get_multidim_indices(model_index_given)
        S = model_index_given - 1
        index_logN_dV = div(S, (len_nH * len_Tg)) + 1
        S %= len_nH * len_Tg
        index_lognH = div(S, len_Tg) + 1
        S %= len_Tg
        index_Tg = S + 1
        return (
            index_logN_dV,
            index_lognH,
            index_Tg
        )
    end
    max_flux = maximum(abs.(obs) .+ 1.0 ./ inv_sigma_obs)
    min_diff = Inf
    min_diff_id = 1
    diff = 0.0
    i = 0
    new_modelTbs = zeros(Float64, length(modelTbs[1, :]))
    for model_id in model_ids
        i += 1
        if model_id == 0
            continue
        end
        new_modelTbs .= fill_fac .* modelTbs[model_id, :]
        @views if any(abs.(new_modelTbs[absent_lines_list]) .> threshold[absent_lines_list]) || any(x -> abs(x) > max_flux, new_modelTbs[line_list])
            continue
        end
        index_logN_dV, index_lognH, index_Tg = get_multidim_indices(i)
        cloud_size = 10^(SPECIFIC_COLUMN_DENSITIES[index_logN_dV] - HDENSITIES[index_lognH]) * fwhm * CLOUD_SIZE_FACTOR
        if GAS_TEMPERATURES[index_Tg] > max_Tgas || !(MINIMUM_CLOUD_SIZE_PROJECTED_ON_SKY <= cloud_size <= MAXIMUM_CLOUD_SIZE_PROJECTED_ON_SKY) || SPECIFIC_COLUMN_DENSITIES[index_logN_dV] > MAX_SPECIFIC_COLUMN_DENSITY
            continue
        end
        diff = 0.0
        for j in eachindex(line_list)
            diff += ((new_modelTbs[line_list[j]] - obs[j]) * inv_sigma_obs[j])^2
        end
        if diff < min_diff
            min_diff = diff
            min_diff_id = i
        end
    end
    index_logN_dV, index_lognH, index_Tg = get_multidim_indices(min_diff_id)
    return [SPECIFIC_COLUMN_DENSITIES[index_logN_dV], HDENSITIES[index_lognH], GAS_TEMPERATURES[index_Tg], fill_fac]
end

# Estimate physical conditions using MCMC
function emcee_phys(obs::Vector{Float64}, inv_sigma_obs::Vector{Float64}, ln2pisigma::Vector{Float64}, threshold::Vector{Float64}, line_list::Vector{Int64}, absent_lines_list::Vector{Int64}, max_Tgas::Float64, fwhm::Float64, fill_fac_in::Float64 = 1.0, fill_fac_err::Float64 = Inf)
    max_flux = maximum(abs.(obs) .+ 1.0 ./ inv_sigma_obs)
    new_modelTbs = zeros(Float64, length(modelTbs[1, :]))
    function log_prob_phys(p)
        logN_dV0, lognH0, Tg0, fill_fac = p[1], p[2], p[3], p[4]
        cloud_size = 10^(logN_dV0 - lognH0) * fwhm * CLOUD_SIZE_FACTOR
        if logN_dV0 < SPECIFIC_COLUMN_DENSITIES[1] || logN_dV0 > SPECIFIC_COLUMN_DENSITIES[end] || logN_dV0 > MAX_SPECIFIC_COLUMN_DENSITY
            return -Inf
        end
        if lognH0 < HDENSITIES[end] || lognH0 > HDENSITIES[1]
            return -Inf
        end
        if Tg0 < GAS_TEMPERATURES[1] || Tg0 > GAS_TEMPERATURES[end] || Tg0 > max_Tgas
            return -Inf
        end
        if !(MINIMUM_CLOUD_SIZE_PROJECTED_ON_SKY <= cloud_size <= MAXIMUM_CLOUD_SIZE_PROJECTED_ON_SKY)
            return -Inf
        end
        if !(1.e-8 <= fill_fac <= 1.0)
            return -Inf
        end
        logN_dV = logN_dV0 - SPECIFIC_COLUMN_DENSITIES[1]
        index_logN_dV = round(Int, logN_dV / DlogN_dV)
        lognH = lognH0 - HDENSITIES[1]
        index_lognH = round(Int, lognH / DlognH)
        Tg = Tg0 - GAS_TEMPERATURES[1]
        index_Tg = round(Int, Tg / DTg)
        model_index = 1 + index_logN_dV * len_nH * len_Tg + index_lognH * len_Tg + index_Tg
        if model_ids[model_index] == 0
            return -Inf
        end
        new_modelTbs .= fill_fac .* modelTbs[model_ids[model_index], :]
        if @views any(abs.(new_modelTbs[absent_lines_list]) .> threshold[absent_lines_list]) || any(x -> abs(x) > max_flux, new_modelTbs[line_list])
            return -Inf
        end
        diff = 0.0
        for i in eachindex(line_list)
            diff += ((new_modelTbs[line_list[i]] - obs[i]) * inv_sigma_obs[i])^2 #+ ln2pisigma[i]
        end
        diff += ((fill_fac - fill_fac_in) / fill_fac_err)^2
        return - 0.5 * diff
    end
    res = get_closest_phys_model(obs, inv_sigma_obs, line_list, absent_lines_list, threshold, max_Tgas, fwhm, fill_fac_in)
    local_rng = Random.TaskLocalRNG()
    Random.seed!(local_rng, 11051989)
    numwalkers = 2 * length(res) + 2
    burnin = 10000
    numsamples_perwalker = 1000000
    p0 = zeros(numwalkers, length(res))
    p0[1, :] .= copy(res)
    normal_dist = Normal(0.0, 1.0)
    p_cand = zeros(length(res))
    res_modified = copy(res)
    for iw in 2:numwalkers
        while true
            p_cand .= res_modified .* (1.0 .+ 1.e-6 .* rand(local_rng, normal_dist, length(p_cand)))
            if isfinite(log_prob_phys(@view(p_cand[:])))
                p0[iw, :] .= p_cand
                break
            end
        end
    end
    chain = sample(log_prob_phys, numwalkers, p0, burnin, 1, rng=local_rng)
    chain = sample(log_prob_phys, numwalkers, chain[end, :, :], numsamples_perwalker, 1, rng=local_rng)
    #return flattenmcmcarray(chain), rhat(chain), res
    return flattenmcmcarray(chain), [], res
end

function get_physical_conditions(cell)
    lines_list = Vector{Int64}()
    absent_lines_list = Vector{Int64}()
    yobs = Vector{Float64}()
    Tobs_dV = Vector{Float64}()
    Tobs_dV_err = Vector{Float64}()
    linewidth_obs = Vector{Float64}()
    sigmaobs = Vector{Float64}()
    for gauss_model in cell.gaussians
        min_freq_diff = 1.e60
        min_freq_id = -1
        for line_id in eachindex(MODEL_LINE_FREQS)
            freq_diff = abs(gauss_model.freq - MODEL_LINE_FREQS[line_id])
            if freq_diff < min_freq_diff
                min_freq_diff = freq_diff
                min_freq_id = line_id
            end
        end
        if !(min_freq_id in MODEL_FREQ_INDEXES_TO_NOT_INCLUDE_IN_PHYSICAL_COND_ESTIMATION)
            push!(lines_list, min_freq_id)
            push!(yobs, gauss_model.p[1])
            push!(linewidth_obs, gauss_model.p[3] / gauss_model.p[2])
            push!(Tobs_dV, gauss_model.p[1] * linewidth_obs[end] * 3.e5 * 1.665)
            push!(Tobs_dV_err, Tobs_dV[end] * sqrt((gauss_model.err_p[1] / gauss_model.p[1])^2. + (gauss_model.err_p[3] / gauss_model.p[3])^2.))
            push!(sigmaobs, gauss_model.err_p[1])
        end
    end
    if length(lines_list) < 5
        println("#warning not enough lines for cell $(cell.id)")
        return nothing
    end

    max_Tgas = (maximum(linewidth_obs) * 3.e10)^2.0 * (MOLECULAR_MASS / (1.380649e-16 * 6.02214076e23) * 0.5)
    cell.fwhm = mean(linewidth_obs) * 3.e5 * 1.665

    absent_lines_list = setdiff(1:length(MODEL_LINE_FREQS), lines_list)
    for freq_to_not_include in MODEL_FREQ_INDEXES_TO_NOT_INCLUDE_IN_ABSENT_LINE_LIST
        absent_lines_list = filter(x -> x != freq_to_not_include, absent_lines_list)
    end

    flux_limits = ones(length(MODEL_LINE_FREQS)) * cell.rms * THRESHOLD_IN_SIGMAS
    for line_id in eachindex(MODEL_LINE_FREQS)
        freq_shifted = MODEL_LINE_FREQS[line_id] * (1 - cell.vlsr / 3.e5)
        for (obs_X, obs_Y) in zip(cell.obs_lines_freq, cell.obs_lines)
            width_to_to_find_flux_limits_in_channels = ceil(Int, (2 * cell.fwhm / 3.e5 * freq_shifted) / (obs_X[2] - obs_X[1]))
            if obs_X[1] <= freq_shifted <= obs_X[end]
                center_id = 1
                for i in 1:length(obs_X) - 1
                    if obs_X[i] <= freq_shifted <= obs_X[i + 1]
                        center_id = i
                        break
                    end
                end
                left_point = max(center_id - width_to_to_find_flux_limits_in_channels, 1)
                right_point = min(center_id + width_to_to_find_flux_limits_in_channels, length(obs_X))
                flux_limits[line_id] = maximum(abs.(obs_Y[left_point:right_point]))
                break
            end
        end
    end

    cell.num_of_lines_to_fit = length(lines_list)

    flat_samples, cell.rhat_phys, min_pars = emcee_phys(yobs, 1.0 ./ sigmaobs, log.((2 * pi) .* (sigmaobs .* sigmaobs)), flux_limits, lines_list, absent_lines_list, max_Tgas, cell.fwhm)
    fill_fac_kde = kde(flat_samples[4, :])
    imax = findmax(fill_fac_kde.density)[2]
    cell.fill_fac = fill_fac_kde.x[imax]
    flat_samples, cell.rhat_phys, min_pars = emcee_phys(yobs, 1.0 ./ sigmaobs, log.((2 * pi) .* (sigmaobs .* sigmaobs)), flux_limits, lines_list, absent_lines_list, max_Tgas, cell.fwhm, cell.fill_fac, (percentile(@view(flat_samples[4, :]), 84) - percentile(@view(flat_samples[4, :]), 16)) * 0.5)

    cell.logN_dV = median(@view(flat_samples[1, :]))
    cell.lognH = median(@view(flat_samples[2, :]))
    cell.Tg = median(@view(flat_samples[3, :]))
    cell.fill_fac = median(@view(flat_samples[4, :]))
    cell.cloud_size = 10^(cell.logN_dV - cell.lognH) * cell.fwhm * CLOUD_SIZE_FACTOR / 1.496e+13
    cell.logN = cell.logN_dV + log10(cell.fwhm * 1.e5)
    #println("phys rhat $(cell.id): ", cell.rhat_phys)
    println("phys min khi2 pars $(cell.id): ", min_pars, " cloud size in AU: ", cell.cloud_size)
    println("phys pars $(cell.id): $(cell.logN_dV), $(cell.lognH), $(cell.Tg), $(cell.fill_fac), cloud size in AU: ", cell.cloud_size)

    theorTb, cell.tau, cell.min_tau = get_phys_model(cell.logN_dV, cell.lognH, cell.Tg, lines_list)
    if isa(theorTb, Vector)
        theorTb = theorTb .* cell.fill_fac
        fig = CairoMakie.Figure(size = (1500, 1200), fontsize=18)
        ax = CairoMakie.Axis(fig[1, 1])
        CairoMakie.scatter!(ax, MODEL_LINE_FREQS[lines_list], theorTb; color=:black, markersize=15)
        CairoMakie.errorbars!(ax, MODEL_LINE_FREQS[lines_list], yobs, sigmaobs; color = :red)
        ax.xlabel = "Frequency, GHz"
        ax.ylabel = "Intensity, K"
        ax.title = "Cell id = $(cell.id)"
        CairoMakie.save(joinpath(DATA_DIR, "figures/fitspectra/$(cell.id).png"), fig)
    elseif isnan(theorTb)
        println("#warning bad model for $(cell.id): ", theorTb)
        cell.logN_dV = NaN ; cell.lognH = NaN ; cell.Tg = NaN ; cell.fill_fac = NaN
    end

    function get_phys_err(arr)
        perc1 = percentile(arr, 16)
        perc2 = percentile(arr, 84)
        return (perc2 - perc1) * 0.5
    end

    cell.err_logN_dV = get_phys_err(@view(flat_samples[1, :]))
    cell.err_lognH = get_phys_err(@view(flat_samples[2, :]))
    cell.err_Tg = get_phys_err(@view(flat_samples[3, :]))
    cell.fill_fac_err = get_phys_err(@view(flat_samples[4, :]))

    cell.Nrot, cell.Nrot_err, cell.Trot, cell.Trot_err, cell.fill_fac_rot, cell.fill_fac_rot_err = perform_rotational_analysis(MODEL_LINE_FREQS[lines_list] .* 1.e9, Eu_K[lines_list], Au[lines_list], gu[lines_list], abs.(Tobs_dV), Tobs_dV_err, linewidth_obs .* (3.e10 * 1.665), cell.id, cell.logN, cell.err_logN_dV + log10(cell.fwhm * 1.e5), cell.Tg, cell.err_Tg, cell.fill_fac, cell.fill_fac_err)
    #println("phys rot pars $(cell.id): $(cell.Nrot) $(cell.Trot) $(cell.fill_fac_rot)")
    #fig_phys = pairplot(flat_samples')
    #save(joinpath(DATA_DIR, "figures/corner_plot_$(cell.id).png"), fig_phys)
end

function process_cell(cell::Cell)
    cell.obs_lines_freq, cell.obs_lines = separate_non_nan_ranges(cell.freqs, cell.obs_spec)
    #=local_rng = Random.TaskLocalRNG()  # Test line profile
    Random.seed!(local_rng, 11051989)
    normal_dist = Normal(0, cell.rms)
    width1 = 4.3 / (2. * sqrt(log(2.))) / 3.e5 * 220.57
    width2 = 3.3 / (2. * sqrt(log(2.))) / 3.e5 * 220.57
    for (xval, yval) in zip(cell.obs_lines_freq, cell.obs_lines)
        if (xval[1] <= 220.57 <= xval[end]) ||  (xval[1] >= 220.57 >= xval[end])
            yval .= zeros(length(yval))
            yval .+= gauss_profile(xval, [9.7, 220.57, width1])
            yval .+= gauss_profile(xval, [5.5, 220.57*(1-6.5/3.e5), width2])
            yval .+= rand(local_rng, normal_dist, length(yval))
            cell.rms = std(yval[1:14], corrected=true)
        else
            yval .= rand(local_rng, normal_dist, length(yval))
        end
    end=#
    cell.gaussians, cell.all_gaussians, cell.vlsr, cell.err_vlsr, cell.max_vlsr_diff = get_gaussians(cell.obs_lines_freq, cell.obs_lines, cell.rms, cell.id)
    #=for gauss in cell.all_gaussians
        println(gauss.p)
        println(gauss.err_p)
        println(cell.rms)
    end=#
    get_physical_conditions(cell)
end

function main()
    file_path = joinpath(DATA_DIR, FITS_FOR_NOISE_ESTIMATES)
    hdul = FITSIO.FITS(file_path, "r")
    header = FITSIO.read_header(hdul[1])
    image_data_cont = FITSIO.read(hdul[1])
    FITSIO.close(hdul)

    file_path = joinpath(DATA_DIR, FITS_FILE)
    hdul = FITSIO.FITS(file_path, "r")
    header = FITSIO.read_header(hdul[1])
    header_to_write = header
    image_data = FITSIO.read(hdul[1])
    channels = 1:header["NAXIS3"]
    frequencies = collect(header["CRVAL3"] .+ header["CDELT3"] .* (channels .- header["CRPIX3"]))
    frequencies = frequencies ./ 1e9
    # Conversion factor
    from_Jy_beam_to_K = 1.222e3 * 1.0e3 / (header["RESTFRQ"]^2 / 1.0e18) / (header["BMAJ"] * header["BMIN"] * 3600 * 3600)
    FITSIO.close(hdul)

    # Loop through the image data dimensions
    cell_id = 1
    cells = Vector{Cell}()
    for i in 1:size(image_data, 1)
        for j in 1:size(image_data, 2)
            spectrum = image_data[i, j, :] * from_Jy_beam_to_K
            if isnan(spectrum[1])
                continue
            end
            spectrum_cont = image_data_cont[i, j, CONT_CHANNELS_FOR_NOISE_ESTIMATES_VECTOR] * from_Jy_beam_to_K
            rms = std(spectrum_cont, corrected=true)
            spectrum[1:130] .= NaN
            spec = spectrum[419:439]
            max_line_flux_12_5_11_5 = maximum(spec)
            spec = spectrum[461:475]
            max_line_flux_12_4_11_4 = maximum(spec)
            spec = spectrum[490:505]
            max_line_flux_12_3_11_3 = maximum(spec)
            spec = spectrum[517:527]
            max_line_flux_12_2_11_2 = maximum(spec)
            spec = spectrum[530:547]
            max_line_flux_12_0_11_0 = maximum(spec)
            if max_line_flux_12_5_11_5 > THRESHOLD_IN_SIGMAS * rms && max_line_flux_12_4_11_4 > THRESHOLD_IN_SIGMAS * rms && max_line_flux_12_3_11_3 > THRESHOLD_IN_SIGMAS * rms && max_line_flux_12_2_11_2 > THRESHOLD_IN_SIGMAS * rms && max_line_flux_12_0_11_0 > THRESHOLD_IN_SIGMAS * rms
                cells = push!(cells, Cell(cell_id, frequencies, spectrum, rms, i, j))
                cell_id += 1
            end
        end
    end
    GC.gc()

    Threads.@threads :greedy for cell in cells
        process_cell(cell)
    end
    GC.gc()

    open(joinpath(DATA_DIR, "cells_G12.jls"), "w") do file
        serialize(file, cells)
    end

    image_NdV = copy(image_data[:, :, 1]) ; image_NdV .= NaN
    image_N_rot = copy(image_data[:, :, 1]) ; image_N_rot .= NaN
    image_T_rot = copy(image_data[:, :, 1]) ; image_T_rot .= NaN
    image_N = copy(image_data[:, :, 1]) ; image_N .= NaN
    image_nH2 = copy(image_data[:, :, 1]) ; image_nH2 .= NaN
    image_Tgas = copy(image_data[:, :, 1]) ; image_Tgas .= NaN
    image_tau = copy(image_data[:, :, 1]) ; image_tau .= NaN
    image_mintau = copy(image_data[:, :, 1]) ; image_mintau .= NaN
    image_vlsr = copy(image_data[:, :, 1]) ; image_vlsr .= NaN
    image_cellid = copy(image_data[:, :, 1]) ; image_cellid .= NaN
    image_num_lines = copy(image_data[:, :, 1]) ; image_num_lines .= NaN
    image_fwhm = copy(image_data[:, :, 1]) ; image_fwhm .= NaN
    image_cloud_size = copy(image_data[:, :, 1]) ; image_cloud_size .= NaN
    image_max_vlsr_diff = copy(image_data[:, :, 1]) ; image_max_vlsr_diff .= NaN
    image_rmsT = copy(image_data[:, :, 1]) ; image_rmsT .= NaN
    image_fillfac = copy(image_data[:, :, 1]) ; image_fillfac .= NaN

    image_errNdV = copy(image_data[:, :, 1]) ; image_errNdV .= NaN
    image_errNrot = copy(image_data[:, :, 1]) ; image_errNrot .= NaN
    image_errnH2 = copy(image_data[:, :, 1]) ; image_errnH2 .= NaN
    image_errTgas = copy(image_data[:, :, 1]) ; image_errTgas .= NaN
    image_errTrot = copy(image_data[:, :, 1]) ; image_errTrot .= NaN
    image_errvlsr = copy(image_data[:, :, 1]) ; image_errvlsr .= NaN
    image_errfillfac = copy(image_data[:, :, 1]) ; image_errfillfac .= NaN

    for one_cell in cells
        i = one_cell.x
        j = one_cell.y
        image_NdV[i, j] = one_cell.logN_dV ; image_errNdV[i, j] = one_cell.err_logN_dV
        image_N[i, j] = one_cell.logN ; #image_errNdV[i, j] = one_cell.err_logN_dV
        image_nH2[i, j] = one_cell.lognH ; image_errnH2[i, j] = one_cell.err_lognH
        image_Tgas[i, j] = one_cell.Tg ; image_errTgas[i, j] = one_cell.err_Tg
        image_vlsr[i, j] = one_cell.vlsr - VLSR_IN_KMS ; image_errvlsr[i, j] = one_cell.err_vlsr
        image_tau[i, j] = one_cell.tau
        image_mintau[i, j] = one_cell.min_tau
        image_cellid[i, j] = one_cell.id
        image_num_lines[i, j] = one_cell.num_of_lines_to_fit
        image_fwhm[i, j] = one_cell.fwhm
        image_cloud_size[i, j] = one_cell.cloud_size
        image_max_vlsr_diff[i, j] = one_cell.max_vlsr_diff
        image_rmsT[i, j] = one_cell.rms
        image_N_rot[i, j] = one_cell.Nrot ; image_errNrot[i, j] = one_cell.Nrot_err
        image_T_rot[i, j] = one_cell.Trot ; image_errTrot[i, j] = one_cell.Trot_err
        image_fillfac[i, j] = one_cell.fill_fac ; image_errfillfac[i, j] = one_cell.fill_fac_err
    end

    function write_into_fits_file(arr, name, input_header)
        f = FITSIO.FITS(joinpath(DATA_DIR, name), "w")
        write(f, arr, header=input_header)
        FITSIO.close(f)
    end

    header_to_write["NAXIS3"] = 1

    header_to_write["BTYPE"] = "log(N/dV)"
    header_to_write["BUNIT"] = "[cm^-3 s]"
    write_into_fits_file(image_NdV, "image_NdV_G12.fits", header_to_write)
    header_to_write["BTYPE"] = "Nrot"
    header_to_write["BUNIT"] = "cm^-3 s"
    write_into_fits_file(image_N_rot, "image_Nrot_G12.fits", header_to_write)
    header_to_write["BTYPE"] = "log(N)"
    header_to_write["BUNIT"] = "[cm^-2]"
    write_into_fits_file(image_N, "image_N_G12.fits", header_to_write)
    header_to_write["BTYPE"] = "log(nH2)"
    header_to_write["BUNIT"] = "[cm^-3]"
    write_into_fits_file(image_nH2, "image_nH2_G12.fits", header_to_write)
    header_to_write["BTYPE"] = "Tgas"
    header_to_write["BUNIT"] = "K"
    write_into_fits_file(image_Tgas, "image_Tgas_G12.fits", header_to_write)
    header_to_write["BTYPE"] = "Trot"
    header_to_write["BUNIT"] = "K"
    write_into_fits_file(image_T_rot, "image_Trot_G12.fits", header_to_write)
    header_to_write["BTYPE"] = "tau"
    header_to_write["BUNIT"] = "tau"
    write_into_fits_file(image_tau, "image_tau_G12.fits", header_to_write)
    header_to_write["BTYPE"] = "tau"
    header_to_write["BUNIT"] = "tau"
    write_into_fits_file(image_mintau, "image_mintau_G12.fits", header_to_write)
    header_to_write["BTYPE"] = "VLSR"
    header_to_write["BUNIT"] = "km/s"
    write_into_fits_file(image_vlsr, "image_vlsr_G12.fits", header_to_write)
    header_to_write["BTYPE"] = "cell id"
    header_to_write["BUNIT"] = "cell id"
    write_into_fits_file(image_cellid, "image_cellid_G12.fits", header_to_write)
    header_to_write["BTYPE"] = "number of lines"
    header_to_write["BUNIT"] = "number of lines"
    write_into_fits_file(image_num_lines, "image_num_lines_G12.fits", header_to_write)
    header_to_write["BTYPE"] = "FWHM"
    header_to_write["BUNIT"] = "km/s"
    write_into_fits_file(image_fwhm, "image_fwhm_G12.fits", header_to_write)
    header_to_write["BTYPE"] = "cloud size"
    header_to_write["BUNIT"] = "AU"
    write_into_fits_file(image_cloud_size, "image_cloud_size_G12.fits", header_to_write)
    header_to_write["BTYPE"] = "max VLSR diff"
    header_to_write["BUNIT"] = "km/s"
    write_into_fits_file(image_max_vlsr_diff, "image_max_vlsr_diff_G12.fits", header_to_write)
    header_to_write["BTYPE"] = "Intensity error"
    header_to_write["BUNIT"] = "K"
    write_into_fits_file(image_rmsT, "image_rmsT_G12.fits", header_to_write)
    header_to_write["BTYPE"] = "filling factor"
    header_to_write["BUNIT"] = "omega"
    write_into_fits_file(image_fillfac, "image_FillFacT_G12.fits", header_to_write)

    header_to_write["BTYPE"] = "error log(N/dV)"
    header_to_write["BUNIT"] = "[cm^-3 s]"
    write_into_fits_file(image_errNdV, "image_errNdV_G12.fits", header_to_write)
    header_to_write["BTYPE"] = "error Nrot"
    header_to_write["BUNIT"] = "cm^-3 s"
    write_into_fits_file(image_errNrot, "image_errNrot_G12.fits", header_to_write)
    header_to_write["BTYPE"] = "error log(nH2)"
    header_to_write["BUNIT"] = "[cm^-3]"
    write_into_fits_file(image_errnH2, "image_errnH2_G12.fits", header_to_write)
    header_to_write["BTYPE"] = "error Tgas"
    header_to_write["BUNIT"] = "K"
    write_into_fits_file(image_errTgas, "image_errTgas_G12.fits", header_to_write)
    header_to_write["BTYPE"] = "error Trot"
    header_to_write["BUNIT"] = "K"
    write_into_fits_file(image_errTrot, "image_errTrot_G12.fits", header_to_write)
    header_to_write["BTYPE"] = "error VLSR"
    header_to_write["BUNIT"] = "km/s"
    write_into_fits_file(image_errvlsr, "image_errvlsr_G12.fits", header_to_write)
    header_to_write["BTYPE"] = "filling factor error"
    header_to_write["BUNIT"] = "omega"
    write_into_fits_file(image_errfillfac, "image_errFillFacT_G12.fits", header_to_write)
end

main()
