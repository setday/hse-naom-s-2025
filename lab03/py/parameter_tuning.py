import itertools
import polars as pl

from MS import MomentumStrategy


def tune_params(market_data: pl.DataFrame):
    # Define the parameter grids
    a_values = [2, 2.5, 3]
    alpha1_values = [1, 3, 5]
    alpha2_values = [1, 3, 5]
    alpha3_values = [1, 3, 5]
    alpha4_values = [1, 3, 5]
    beta_values = [1, 3, 5]

    # Generate all combinations of parameters
    param_grid = itertools.product(a_values, alpha1_values, alpha2_values, alpha3_values, alpha4_values, beta_values)

    results = []
    for params in param_grid:
        a, alpha1, alpha2, alpha3, alpha4, beta = params

        # Instantiate the strategy with the current parameters
        test = MomentumStrategy(K=144,
                                risk_control=1,
                                window_size=200,
                                data=market_data,
                                a=a,
                                beta=beta,
                                alpha1=alpha1,
                                alpha2=alpha2,
                                alpha3=alpha3,
                                alpha4=alpha4)

        # Get the Sortino Raio for the current set of parameters
        sortino = test.get_sortino_ratio(T=market_data.height)
        print(f"INFO: Sortino Ratio: {round(sortino, 3)}")

        results.append({
            "a": a,
            "alpha1": alpha1,
            "alpha2": alpha2,
            "alpha3": alpha3,
            "alpha4": alpha4,
            "beta": beta,
            "sortino": sortino
        })

    sorted_results = sorted(results, key=lambda x: x["sortino"], reverse=True)
    best_result = sorted_results[0]
    print("Best parameters:", best_result)


if __name__ == "__main__":
    market_data = pl.read_csv("../data/train.csv")
    tune_params(market_data=market_data)
