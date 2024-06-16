using Distributions
using Distributed
using Plots
using Statistics
using StatsBase
using Random
using ArgParse


#=
Simplest possible simulation - just calculating the possibility (for each R0) that there is only one detected case
This is approximately the same as those simulations accpeted at 90 days after the first detected case
=#

function write_file(filename, A, B)
    open(filename, "w") do file
        # Iterate over the elements of A and B
        for i in 1:length(A)
            # Write A[i] and B[i] to the file, formatted as "A[i] B[i]\n"
            write(file, string(A[i], " ", B[i], "\n"))
        end
    end
    
    
end

Random.seed!(1)

function parse_commandline()
    
    # Parse command line options and parameters
    
    s = ArgParseSettings()
    
    @add_arg_table! s begin 
    "-o"
    default="julia_output"
    help="base filename"
    
    "-n"
    arg_type=Int64
    default=100000
    help="number of samples per R0"
end
return parse_args(s)
end


function simulate()

    args = parse_commandline()


    T_sim = 10000
    a1 = 7.402580682098853
    b1 = 1.7375081458523993
    a2 = 1.03138989929643
    b2 = 1.0025194828961315
    time_to_detection = 18

    p_detect = 0.1
    max_detected = 1
    max_cases = 2000

    rng_symptom() = floor(Int64, rand(Weibull(a1, b1))+0.5)
    rng_infect() = floor(Int64,rand(Weibull(a2, b2))+0.5)

    test_symptom = [rng_symptom() for i in 1:100000]
    @show(mean(test_symptom), std(test_symptom))

    test_infect = [rng_infect() for i in 1:100000]
    @show(mean(test_infect), std(test_infect))


    function sample(R0=3.0, N=1000)
        result = []
        for i in 1:N
            n_detected = 0
            n_cases = 0
            first_detection_time = nothing

            population = [(0, rng_symptom(), false)]
            old_population = []
            p = popfirst!(population)
            time = p[1]
            push!(old_population, p)
            n_cases += 1

            fail = false
            while time <= T_sim && n_cases <= max_cases
                n_infect = rand(Poisson(R0))
                for j = 1:n_infect
                    n_cases += 1
                    t_infect = p[2] + rng_infect()
                    t_symptom = t_infect + rng_symptom()
                    detected = rand() < p_detect
                    push!(population, (t_infect, t_symptom, detected))
                    if detected
                        t_detect = t_symptom + time_to_detection
                        n_detected += 1
                        if isnothing(first_detection_time) || t_detect < first_detection_time
                            first_detection_time = t_detect
                        end
                        if n_detected > max_detected
                            fail = true
                            break
                        end
                    end
                end

                if isempty(population) || fail
                    break
                end
                sort!(population, by=x->x[1])
                p = popfirst!(population)
                push!(old_population, p)
                time = p[1]
            end
            push!(result, (fail, first_detection_time, time, n_cases, n_detected))
        end
        result
    end

    R0_vals = LinRange(0.1, 4.0, 40)
    @show(R0_vals)
    results = fill([], 40)
    Threads.@threads for i in 1:40
        R0 = R0_vals[i]
	results[i] = sample(R0, args["n"])
    end

    p1 = plot(R0_vals, [sum([!s[1] && s[5] == 1 for s in r]) for r in results], title="Detection Summary")

    acc = [sum([(!s[1] && s[5]==1) for s in r]) for r in results]
    
    write_file(args["o"] * "/julia_R0_acceptance_rates.dat", R0_vals, acc)

    
    results_flat = reduce(vcat, results)
    
    #=
    p_cases = [r[4] for r in results_flat if !r[1] && r[5] == 1]
    if !isempty(p_cases)
        p2 = histogram(p_cases, bins=0:50, title="Case Distribution")
    end
    =#

    p_days = [r[2] for r in results_flat if !r[1] && r[5] == 1]
    c = counts(p_days, 0:50)
    c = c/sum(c)

    write_file(args["o"] * "/julia_output_time.dat", 0:50, c)

#=
    @show(c)
    if !isempty(p_days)
        p3 = histogram(p_days, bins=0:50, title="Day Distribution")
    end

    plot(p1, p2, p3, layout=(1, 3), size=(1200,400))
    savefig(args["o"] * "/output.png")

=#

end

simulate()
