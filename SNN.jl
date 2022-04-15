push!(LOAD_PATH, "./src")

module SNN


using LinearAlgebra
using Random
using Models
using Testing
using Training
using Neurons



export train, test

function mutation_rate(epoch)
	R₀ = 0.05
	tₕ = 25 # half life is 25 epochs
	return R₀*exp(-epoch/tₕ)
end

function train(M::Model, f_in::Vector{Function}, f_out::Vector{Function}, dt::Float64, final_t::Float64, P_mutation = mutation_rate, grad::Bool = false, pop::Int = 30, n_epochs::Int = 50, killrate::Float64 = 0.3)
	n_tasks = size(f_in)[1]
	println("Training network for $n_tasks different tasks.")
	
	# initialize reservoir and voltages
	iRes = M.iRes
	ΔV_L_vecs = Vector{Vector{Float64}}(undef,n_tasks)
	ΔV_R_vecs = Vector{Vector{Float64}}(undef,n_tasks)
	if(iRes == 0)
		println("Training without reservoir, no gradient descent for voltages")
		grad = false
	else
		n_res = size(M.W[iRes,iRes])[1]
		println("Training with $n_res x $n_res reservoir")
		
		# initialize the voltage vectors
		for ΔV_L_vec in ΔV_L_vecs
			ΔV_L_vec = zeros(n_res)
		end
		for ΔV_R_vec in ΔV_R_vecs
			ΔV_R_vec = zeros(n_res)
		end

		if(grad)
			println("Gradient descent loop on, voltages initialized randomly")
			# 0.5 keeps the max ΔV over a synapse under 1 Volt
			Vmax = 0.5
			for ΔV_L_vec in ΔV_L_vecs
				ΔV_L_vec = Vmax*rand(n_res)
			end
			for ΔV_R_vec in ΔV_R_vecs
				ΔV_R_vec = Vmax*rand(n_res)
			end
		else
			println("Gradient descent loop off, voltages set to 0")
		end
	end
	println("Reservoir and voltages initialized")
	println("Creating population of $pop mutated models...")
	Population = Vector{Model}(undef,pop)
	for model in Population
		model = mutate(M,0.5)
	end
	Population[1] = deepcopy(M)
	
	# setting up loop over all epochs, trials, and models for forward prop, error, and mutation
	t = [t for t in 0:dt:maxt]
	println("\n=================== EPOCHS LOOP =================\n")
	for epoch = 1:n_epochs
		err = zeros(pop)
		println("Epoch $epoch: #trials = $n_trials")
		for trial in 1:n_trials
			itask = rand(1:n_tasks)
			f = f_in[itask]
			f_expected = f_out[itask]
			for imodel in eachindex(Population)
				println("Trial $trial for model #$imodel")
				model = Population[imodel]
				output = ForwardProp(model,f,dt,maxt)
				expected = f_expected.(t)
				err[imodel] += Error(output,expected)
			end
		end
		# sort the population based on their performances
		sortError = sortperm(err)
		Population = Population[sortError]
		println("Lowest error: $(minimum(err))")
		# now kill and reproduce
		println("Performing mutations at rate $(P_mutation(epoch)) for pop size $pop and $(killrate*100) death rate")
		for ikill = round(pop - killrate*pop):pop
			iparent = rand(1:killrate*pop)
			Population[ikill] = mutate(Population[iparent],P_mutation(epoch))
		end
	end
	error = minimum(err)
	println("Error = $error")
	return Population[1], error
end

function test(M::Model, f::Function, f_expected::Function, dt::Float64, maxt::Float64, ΔV_L_vec::Vector{Float64}, ΔV_R_vec::Vector{Float64}, visual::Bool=false)
	println("Running model (fix with name) for t = $maxt, dt = $dt")
	println("Voltage vector left and right: ")
	display(ΔV_L_vec)
	display(ΔV_R_vec)
	t = [t for t in 0:dt:maxt]
	Ggen(M,ΔV_L_vec,ΔV_R_vec)
	if(visual)
		output = VisualForwardProp(M, f, dt, maxt)
	else
		output = ForwardProp(M, f, dt, maxt)
	end
	expected = f_expected.(t)
	show(output)
	show(expected)
	err = Error(output,expected)
	return output, expected, t, err
end






end
