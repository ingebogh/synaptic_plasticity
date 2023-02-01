# ``R`` code for inferring synaptic plasticity

This code is used to fit a Bayesian model for inferring synaptic plasticity in spike data. 

## Model 

We assume a binary data model using the Bernoulli likelihood:
  $$s_{i,t} \sim \text{Bern}(\mu_{i,t}) = \mu_{i,t}^{s_{i,t}} (1-\mu_{i,t})^{1-s_{i,t}},$$
where the linear predictor $\eta_{i,t} = \text{logit}(\mu_{i,t})$ is given by
  $$\eta_{i,t} = \sum_{j = 1, j \neq i}^N w_{ij,t} s_{j, t-a} + b_i + \sum_{k = 1}^K \beta_k x_{k,t}.$$
$i$ and $j$ indicate neurons, $a$ is the time lag, $b_i$ is an intercept, and $\beta_k$ is the parameter of the covariate $\mathbf{x}_ k = (x_{k,1}, \dots, x_{k,T})^{\text{T}}$.
$w_{ij,t}$ describes the connection from neuron $j$ to neuron $i$, modelled as a random walk of first order:
  $$w_{ij,t+1} = w_{ij,t} + l(\mathbf{s}_ {i, 1:t}, \mathbf{s}_ {j, 1:t}, \mathbf{\theta}, w_{ij,t}) + \varepsilon_{t+1},$$
where $\varepsilon_{t+1} \sim \mathcal{N}(0, \sigma_w^2)$ is Gaussian noise and $\mathbf{\theta} = (A_+, A_-, \tau_+, \tau_-, w_{\text{min}}, w_{\text{max}})$ are the STDP learning rule parameters.
$l(\cdot)$ is the learning rule, where the user can choose additive:
  $$l(\mathbf{s}_ {i, 1:t}, \mathbf{s}_ {j, 1:t}, \mathbf{\theta}) = A_+ s_{i,t} \sum_{t'=1}^t s_{j,t'} e^{\frac{t'-t}{\tau_+}} - A_- s_{j,t} \sum_{t'=1}^t s_{i,t'} e^{\frac{t'-t}{\tau_-}},$$
or multiplicative:
  $$l(\mathbf{s}_ {i, 1:t}, \mathbf{s}_ {j, 1:t}, \mathbf{\theta}) = (w_{\text{max}} - w_{ij,t}) A_+ s_{i,t} \sum_{t'=1}^t s_{j,t'} e^{\frac{t'-t}{\tau_+}} - (w_{ij,t} - w_{\text{min}}) A_- s_{j,t} \sum_{t'=1}^t s_{i,t'} e^{\frac{t'-t}{\tau_-}}.$$

## Run

To run the model: Download "functions.R" and "run_model.R", open the file "run_model.R" and follow the instructions.

The code includes all necessary functions for running the inference, and some functions for exploring the results.
