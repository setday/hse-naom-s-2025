# Momentum Strategy for Algorithmic Trading

# How to run locally

- Prepare the dataset. See `lab03/data/README.md`.
  
- Download the NVIDIA nvc++ compiler to utilize the NVIDIA GPU from [here](https://developer.nvidia.com/hpc-sdk-downloads).

- Set up the environment by following the [documentation](https://docs.nvidia.com/hpc-sdk/hpc-sdk-install-guide/index.html#install-linux-end-usr-env-settings).


```bash
$ NVARCH=`uname -s`_`uname -m`; export NVARCH
$ NVCOMPILERS=/opt/nvidia/hpc_sdk; export NVCOMPILERS
$ MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/24.11/compilers/man; export MANPATH
$ PATH=$NVCOMPILERS/$NVARCH/24.11/compilers/bin:$PATH; export PATH

$ export PATH=$NVCOMPILERS/$NVARCH/24.11/comm_libs/mpi/bin:$PATH
$ export MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/24.11/comm_libs/mpi/man
```

- Compile and run the code:

```bash
nvc++ -acc cpp/tune_params.cpp -o train.exe
./train.exe
```


## TODO list:

- ~~Rewrite implementation in C++~~
- [+-] Construct a grid of parameters
- ~~Evaluate the Sortino ratio (utility function "of the first order") to find the best parameters on the training data~~
- Make sure that you integrate over some neighborhood of a set of parameters for better robustness (function "of the second order")
- ~~Take into account transaction costs by introducing a new parameter, Gamma~~
- Tune the final model (using the function of the "second order")
- Test model out of sample (download and prepare test dataset)
- ~~Quick fix for C++ version: parse data before fine-tuning rather than at each iteration!~~
- Implement statistical tests (to distinguish between a real trend and mere chance)
- ~~Explore why Sortino Ratio explodes for small omega (window_size) [say ```window_size=100```]~~