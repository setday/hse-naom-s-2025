# Momentum Strategy for Algorithmic Trading

## TODO list:

- Construct a grid of parameters
- Evaluate the Sortino ratio (utility function "of the first order") to find the best parameters on the training data
- Make sure that you integrate over some neighborhood of a set of parameters for better robustness (function "of the second order")
- Take into account transaction costs by introducing a new parameter, Gamma
- Tune the final model (using the function of the "second order")
