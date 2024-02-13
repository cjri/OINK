using Distributions
using Distributed
using Plots
using Statistics
using StatsBase

function simulate()
    T_sim = 10000
    a1 = 7.4026
    b1 = 1.7375
    a2 = 1.0314
    b2 = 1.0025
    time_to_detection = 18

    p_detect = 0.1
    max_detected = 1
    max_cases = 2000

    rng_symptom() = floor(Int64, rand(Weibull(a1, b1))+0.5)
    rng_infect() = floor(Int64,rand(Weibull(a2, b2))+0.5)


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
    #R0_vals = LinRange(0.8, 0.9, 40)
    @show(R0_vals)
    results = fill([], 40)
    Threads.@threads for i in 1:40
        R0 = R0_vals[i]
	results[i] = sample(R0, 10000)
    end

    p1 = plot(R0_vals, [sum([!s[1] && s[5] == 1 for s in r]) for r in results], title="Detection Summary")

    results_flat = reduce(vcat, results)

    p_cases = [r[4] for r in results_flat if !r[1] && r[5] == 1]
    if !isempty(p_cases)
        p2 = histogram(p_cases, bins=0:50, title="Case Distribution")
    end

    p_days = [r[2] for r in results_flat if !r[1] && r[5] == 1]
    @show(counts(p_days, 0:50))
    if !isempty(p_days)
        p3 = histogram(p_days, bins=0:50, title="Day Distribution")
    end

    plot(p1, p2, p3, layout=(1, 3), size=(1200,400))
    savefig("output.png")

    results_flat
end

simulate()
