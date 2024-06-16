using ArgParse
using Distributions
using Distributed
using Plots
using Statistics
using StatsBase
using Random
using ProfileView
using DataStructures

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

function simulate()
    
    args = parse_commandline()
    
    T_sim = 10000
    a1 = 7.402580682098853
    b1 = 1.7375081458523993
    a2 = 1.03138989929643
    b2 = 1.0025194828961315
    time_to_detection = 18
    
    p_detect = 0.04
    max_detected = 1
    
    rng_symptom() = floor(Int64, rand(Weibull(a1, b1))+0.5)
    rng_infect() = floor(Int64,rand(Weibull(a2, b2))+0.5)
    
    min_delay = time_to_detection 

    function sample(R0=3.0, N=1000, t_measure=1)
        result = []
        for i in 1:N
	    Random.seed!(i)
            n_detected = 0
            n_cases = 0
            first_detection_time = nothing
            
            detections = []
	    population = BinaryHeap{Tuple{Int64,Int64,Bool}}(Base.By(x->x[2]), [])
            push!(population, (0, rng_symptom(), false))
            old_population = []
            p = pop!(population)
            time = p[2]
            push!(old_population, p)
            n_cases += 1
            
            fail = false
            bad_detection = false
            while time <= T_sim 
                if (!isnothing(first_detection_time))
                    if (time > (first_detection_time - min_delay - t_measure)) && bad_detection # probably could weaken this (see below)		   
                        fail = bad_detection
                        break
                    end
                    if (time > (first_detection_time - min_delay + t_measure)) && !bad_detection
                        fail = bad_detection
                        break
                    end
                end
                
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
                        push!(detections, t_detect)
                        if isnothing(first_detection_time) 
                            first_detection_time = t_detect
                            bad_detection = false
                        else
                            if (t_detect <= first_detection_time)
                                if first_detection_time - t_detect <= t_measure
                                    bad_detection = true
                                else                              
                                    bad_detection = false
                                end
                                first_detection_time = t_detect
                            else
                                if t_detect - first_detection_time  <= t_measure
                                    bad_detection=true
                                end
                            end
                        end
                    end
                end	
                if isempty(population)
                    if isnothing(first_detection_time) || bad_detection
                        fail=true
                        #@show("Empty", fail, n_detected)
                    end
                    #@show("B3")
                    
                    break
                end
#                sort!(population, by=x->x[2])
                p = pop!(population)
                push!(old_population, p)
                time = p[2]
            end
	    #=
	    sort!(detections)
            if !fail
                @assert((length(detections)==1) || (detections[2]-detections[1]>t_measure)) ## Check detections are correct
                @assert(time>(first_detection_time+t_measure-time_to_detection) || isempty(population)) ## Check we have gone past correct time, or empty pop
            end
            if fail
                #@show(first_detection_time, detections, time)
                #if !isnothing(first_detection_time)
                #		    @show(first_detection_time-time_to_detection-t_measure)
                #end
                @assert(isnothing(first_detection_time) || ((length(detections)>1) && (detections[2] - detections[1] <= t_measure) && (isempty(population) || (time >=first_detection_time-time_to_detection-t_measure )))) # doublecheck this time. Either terminate or go beyond the critical time where not possible to recover a matching simulation
            end
	    =#
            push!(result, (fail, first_detection_time, time, n_cases, n_detected))
        end
        result
    end
    
    R0_vals = LinRange(0.1, 4.0, 40)
    @show(R0_vals)
    results = fill([], 40)
    Threads.@threads for i in 1:40
        R0 = R0_vals[i]
        results[i] = sample(R0, Int64(args["n"]), Int64(args["e"]))
    end


    acc = [sum([!s[1] for s in r]) for r in results]
    
    write_file(args["o"] * "_R0.dat", R0_vals, acc)
    

    results_flat = reduce(vcat, results)

    p_days = [r[2] for r in results_flat if (!r[1])]

    if !isempty(p_days)
        c = counts(p_days, 0:50)
        c = c/length(p_days)
    end

    write_file(args["o"] * "_time.dat", 0:50, c)

end

function parse_commandline()
    
    # Parse command line options and parameters
    
    s = ArgParseSettings()
    
    @add_arg_table! s begin 
    "-o"
    default="julia_output"
    help="base filename"
    
    "-e"
    arg_type=Int64
    default=1
    help="end time"
    "-n"
    arg_type=Int64
    default=10000
    help="number of samples per R0"
end
return parse_args(s)
end


simulate()
