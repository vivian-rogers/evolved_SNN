push!(LOAD_PATH, "./")
using LinearAlgebra
using Random

module Neurons

	

function update_neuron(n::Neuron,dt::Float64) # updates neuron, sends out spike if saturated
	if(n.integrated < n.threshold)
		n.integrated -= dt*n.decay_rate
		return 0
	elseif(true == false)
		#=if(firinglock == true)
			if(n.t_since_sat > n.sat_period) #
				n.lock == true
				n.
			else
				return n.output
		else

		end=#
	end
end


function +(n::Neuron, in::Float64)
	if(n.firinglock == false)
		integrated = max(in,in+integrated)
	end
end

end
