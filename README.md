# Robustness against disturbances in power systems under frequency constraints

This repository archives the code used for the paper  "Robustness against disturbances in power systems under frequency constraints." IEEE Transactions on Control of Network Systems 6.3 (2019): 971-979.

## Description of Files

compute_gain.m: compute system gain matrix for the swing equation.

optimize_bound.jl: optimize the perturbation bound. (It requires MATLAB and Julia with IPOPT and JuMP.)

TSA_SGT_2bus.m: creates plots for analyzing 2-bus system.

## Citation

You can find more details in:

```bibtex
  @article{lee2019robustness,
    title={Robustness against disturbances in power systems under frequency constraints},
    author={Lee, Dongchan and Aolaritei, Liviu and Vu, Thanh Long and Turitsyn, Konstantin},
    journal={IEEE Transactions on Control of Network Systems},
    volume={6}, number={3}, pages={971--979}, year={2019}, publisher={IEEE}
  }
```


