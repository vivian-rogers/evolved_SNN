push!(LOAD_PATH, "./")
push!(LOAD_PATH, "./src")








module Testing

using Models
using LinearAlgebra
using Random
using Neurons

export VisualForwardProp, ForwardProp

struct Network
	M::Model
	Neurons::Vector{Vector{Neuron}}
end


function VisualForwardProp(M::Model,stimulus::Function,dt::Float64,maxt::Float64,Q::Float64=0.5)
	net = Network(M,initNeurons(M.LayerSizes))
	AllSpikes = Any[]
	OutSpikes = BitArray[]
	Times = [t for t in 0:dt:maxt]
	for t in Times 
		# update neurons
		input = stimulus(t)
		all_spikes_t = updateNeurons(net.Neurons)
		run(`clear`)
		show(all_spikes_t)
		push!(AllSpikes,all_spikes_t)
		OutSpikesᵢ = all_spikes_t[M.n]
		push!(OutSpikes,OutSpikesᵢ)
		# encode to pulse train somehow
		# M.Neurons + inSpikes
		for C in net.M.ConnectionList
			L₁ = C[2]; L₀ = C[1];
			Spikesᵢ = all_spikes_t[C[1]]
			feed(net.Neurons[L₁],Q*net.M.G₀[L₁,L₀]*Spikesᵢ)
		end
	end
	# decode OutSpikes to data 
	# Output = Decode(dt,OutSpikes)
	return OutSpikes 
end




function ForwardProp(M::Model,stimulus::Function,dt::Float64,maxt::Float64,Q::Float64=1)
	net = Network(M,initNeurons())
	AllSpikes = Any[]
	OutSpikes = BitArray[]
	Times = [t for t in 0:dt:maxt]
	for t in Times 
		# update neurons
		input = stimulus(t)
		# Spikesᵢ (Vector of bitarrays) = updateNeurons(net.Neurons)
		push!(AllSpikes,Spikesᵢ)
		OutSpikesᵢ = Spikesᵢ[M.n]
		push!(OutSpikes,OutSpikesᵢ)
		# encode to pulse train somehow
		# M.Neurons + inSpikes
		for C in net.M.ConnectionList
			L₁ = C[2]; L₀ = C[1];
			net.Neurons[Lᵢ] + Q*G₀[L₁,L₀]*Spikesᵢ
		end
	end
	# decode OutSpikes to data 
	# Output = Decode(dt,OutSpikes)
	return Output
end

end	
