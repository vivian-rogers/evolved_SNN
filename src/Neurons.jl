
module Neurons

using LinearAlgebra
using Random
export Neuron, initNeurons, updateNeurons, feed

mutable struct Neuron
	threshold::Float64
	refract_pd::Float64
	decay_pd::Float64
	dt::Float64
	# ones that are modified in runtime
	integrated::Float64
	intlock::Bool
end


function initNeurons(LayerSizes::Vector{Int}, threshold::Float64 = 1.0, decay_pd::Float64 = 10^-8, refract_pd::Float64 = 3*10^-9, dt::Float64 = 10^-9)
	n_layers = size(LayerSizes)[1]
	neuron_layers = Vector{Vector{Neuron}}(undef,n_layers)
	for i = 1:n_layers
		nᵢ = LayerSizes[i]
		neuronLayer = Vector{Neuron}(undef,nᵢ)
		for j = 1:nᵢ
			n = Neuron(threshold,refract_pd,decay_pd,dt,0,false)
			neuronLayer[j] = n
		end
		neuron_layers[i] = neuronLayer
	end
	return neuron_layers
end

function updateNeurons(neurons::Vector{Vector{Neuron}})
	n_layers = size(neurons)[1]
	all_spikes = Vector{Vector{Bool}}(undef,n_layers)
	for i in eachindex(neurons)
		layer = neurons[i]
		nᵢ = size(layer)[1]
		layer_spikes = Vector{Bool}(undef,nᵢ)
		#layer_spikes = updateNeuron.(layer)
		for j in 1:nᵢ
			layer_spikes[j] = updateNeuron(layer[j])
		end
		display(layer_spikes)
		display(all_spikes)
		all_spikes[i] = layer_spikes
	end
	return all_spikes
end

function updateNeuron(n::Neuron) # updates neuron, sends out spike if saturated
	dt = n.dt
	if(n.integrated < n.threshold)
		if(n.integrated >= 0)
			n.intlock = false
		else
			n.integrated = max(n.integrated - dt/n.refract_pd,0)
		end
		return 0
	else
		if(n.intlock == true)
			n.integrated = max(n.integrated - dt/n.refract_pd,0)
			return 0
		else # okay, the neuron is past the point and it has not fired. Needs to fire
			
			#if(t_since_sat <
			n.intlock = true
			return 1
		end
	end
end


function feed(neurons::Vector{Neuron}, in_vec::Vector{Float64})
	for i in eachindex(neurons)
		n = neurons[i]
		in = in_vec[i]
		if(n.intlock == false)
			n.integrated = max(in,in+n.integrated)
		end
	end
end

end

