# Coverage characterization of STAR-RIS networks: NOMA and OMA
This simulation code is mainly used to reproduce the results of the following paper [1]:

[1] C. Wu, Y. Liu, X. Mu, X. Gu, and O. A. Dobre, "Coverage characterization of STAR-RIS networks: NOMA and OMA," in IEEE Communication Letter, vol. 25, no. 9, pp. 3036-3040, Sep. 2021, doi: [10.1109/LCOMM.2021.3091807](https://ieeexplore.ieee.org/abstract/document/9462949).

***
If you use this simulation code package in any way, please cite the original paper [1] above. The corresponding BiBTeX citation is given below:
```
@article{wu2021coverage,
  title={Coverage characterization of {STAR-RIS} networks: {NOMA} and {OMA}},
  author={Wu, Chenyu and Liu, Yuanwei and Mu, Xidong and Gu, Xuemai and Dobre, Octavia A},
  journal={IEEE Communications Letters},
  volume={25},
  number={9},
  pages={3036--3040},
  year={2021},
  month=sep,
  publisher={IEEE}
}
```

The author in charge of this simulation code package is: Chenyu Wu (email: wuchenyu@hit.edu.cn).

One could run the `plot_figure3.m` directly see the simulation results. 

Some Notes:

1. `rician_fading.m` is used to generate channel coefficient.

2. noma1 and noma2 refer to different decoding order.

3. The default optimization toolbox is 'fmincon', which can also changed to CVX. Please refer to `noma_cvx.m` and `oma_cvx.m`.


Enjoy the reproducible research!

All Rights Reserved. 
