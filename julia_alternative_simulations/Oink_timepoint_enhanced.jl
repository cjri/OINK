using ArgParse
using Distributions
using Distributed
using Plots
using Statistics
using StatsBase
using Random

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
    p_detect_enhanced = 0.2
    detection_increase_delay = 18

    rng_symptom() = floor(Int64, rand(Weibull(a1, b1))+0.5)
    rng_infect() = floor(Int64,rand(Weibull(a2, b2))+0.5)
    
    
    function sample(R0=3.0, N=1000, t_measure=1)
        result = []
        for i in 1:N
            n_detected = 0
            n_cases = 0
            first_detection_time = nothing
            
            detections = []
	    enhanced_detections = []
            population = [(0, rng_symptom(), false, false)]
            old_population = []
            p = popfirst!(population)
            time = p[2]
            push!(old_population, p)
            n_cases += 1
            
            fail = false
            bad_detection = false # Note that this isn't used
            while time <= T_sim 
                if (!isnothing(first_detection_time))
                    if (time > (first_detection_time - time_to_detection + t_measure)) 
                       break
                    end
                end
                
                n_infect = rand(Poisson(R0))
                for j = 1:n_infect
                    n_cases += 1
                    t_infect = p[2] + rng_infect()
                    t_symptom = t_infect + rng_symptom()
		    r = rand()
                    detected = r < p_detect
		    enhanced_detected = r < p_detect_enhanced
                    push!(population, (t_infect, t_symptom, detected, enhanced_detected))
		    t_detect = t_symptom + time_to_detection

		    if detected
                        n_detected += 1
                        push!(detections, t_detect)
                        if isnothing(first_detection_time) 
                            first_detection_time = t_detect
			    # Add enhanced detections
			    append!(detections, [ e_t for e_t in enhanced_detections if e_t >= first_detection_time + detection_increase_delay ])
			    filter!(u-> u>first_detection_time + detection_increase_delay, enhanced_detections)
                        else
                            if (t_detect <= first_detection_time)
                                first_detection_time = t_detect
			    	append!(detections, [ e_t for e_t in enhanced_detections if e_t >= first_detection_time + detection_increase_delay ])
				filter!(u-> u>first_detection_time + detection_increase_delay, enhanced_detections)
				
			    end
                        end
                    elseif enhanced_detected
		        if !isnothing(first_detection_time) && t_detect>=first_detection_time + detection_increase_delay
			   push!(detections, t_detect)
			else
			   push!(enhanced_detections, t_detect)
			end
		    end
                end	
                if isempty(population)
                    
                    break
                end
                sort!(population, by=x->x[2])
                p = popfirst!(population)
                push!(old_population, p)
                time = p[2]
            end
	    sort!(detections)
		       if isnothing(first_detection_time) ||(n_detected>1 && (detections[2]-detections[1])<=t_measure)
		       
		       fail = true
		       else
		       fail= false
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
    default="julia_output_enhanced"
    help="base filename"
    
    "-e"
    arg_type=Int64
    default=1
    help="end time"
    "-n"
    arg_type=Int64
    default=100000
    help="number of samples per R0"
end
return parse_args(s)
end

simulate()
