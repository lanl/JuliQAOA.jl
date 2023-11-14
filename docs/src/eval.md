# Simulation
```@index
Pages = ["eval.md"]
```

```@docs
statevector(angles, mixer, obj_vals)
statevector!(sv, angles, mixer::Mixer{X}, obj_vals)
probabilities(angles, mixer, obj_vals)
probabilities!(sv, angles, mixer, obj_vals)
exp_value(angles::Vector, mixer::Mixer, obj_vals::AbstractVector, measure::AbstractVector=obj_vals)
exp_value!(sv, angles, mixer, obj_vals, measure)
```