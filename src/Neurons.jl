using LinearAlgebra

module SNN

mutable struct Neuron
	threshold::Float64
	refract_rate::Float64
	sat_period::Float64
	output::Float64	
	decay_rate::Float64
	dt::Float64
	# ones that are modified in runtime
	t_since_sat::Float64
	t_since_spike::Float64
	integrated::Float64
	firinglock::Bool
end


function initNeuron(threshold::Float64 = 1, refract_rate::Float64 = -0.1, 


function updateNeuron(n::Neuron) # updates neuron, sends out spike if saturated
	dt = n.dt
	if(n.integrated < n.threshold)
		n.integrated -= dt*n.decay_rate
		return 0
	elseif
		if(firinglock == true)
			if(n.t_since_sat > n.sat_period) #
				n.lock == true
				n.
			else
				return n.output
		else

		end
	end
end


function +(n::Neuron, in::Float64)
	if(n.firinglock == false)
		integrated = max(in,in+integrated)
	end
end
