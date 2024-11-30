import math
import polars as pl


class MomentumStrategy:

    def __init__(self, K: int = 144, risk_control: float = 10_000,
                 window_size: int = 300, data: pl.DataFrame = None,
                 a: float = None, beta: float = None,
                 alpha1: float = None, alpha2: float = None,
                 alpha3: float = None, alpha4: float = None):

        self.params_: dict = {
            "a": a,
            "beta": beta,
            "alpha1": alpha1,
            "alpha2": alpha2,
            "alpha3": alpha3,
            "alpha4": alpha4
        }

        self.hyperparams_: dict = {
            "K": K,
            "risk_control": risk_control,
            "window_size": window_size,  # Measured in seconds
            "tau": 60  # Timestep to compute position again. Measured in seconds
        }
        self.sell_volume_ = data.get_column("SELL_VOLUME").to_list()
        self.buy_volume_ = data.get_column("BUY_VOLUME").to_list()
        self.volume_ = [s + b for s,
                        b in zip(self.sell_volume_, self.buy_volume_)]
        # The prefix sums below are needed for rapid computation of x3_tilda
        self.volume_prefix_sum = [0]
        for v in self.volume_:
            self.volume_prefix_sum.append(v + self.volume_prefix_sum[-1])
        self.buy_volume_prefix_sum = [0]
        for bv in self.buy_volume_:
            self.buy_volume_prefix_sum.append(
                bv + self.buy_volume_prefix_sum[-1])

        self.mid_px_ = data.get_column("MID_PX").to_list()

    def get_position(self, t: int) -> float:
        if t == 0:
            return 0
        x = self.get_x(t)

        if abs(x) < self.params_["a"]:
            return 0

        if x >= self.params_["a"]:
            return 2 * self.hyperparams_["risk_control"] / math.pi * math.atan(math.pow((x - self.params_["a"]), self.params_["beta"]))

        return -2 * self.hyperparams_["risk_control"] / math.pi * math.atan(math.pow((-x - self.params_["a"]), self.params_["beta"]))

    def get_x(self, t: int):
        x1 = self.get_x1(t)
        x2 = self.get_x2(t)
        x3 = self.get_x3(t)
        x4 = self.get_x4(t)
        return self.params_["alpha1"] * x1 + self.params_["alpha2"] * x2 + self.params_["alpha3"] * x3 + self.params_["alpha4"] * x4

    def get_x1(self, t: int):
        return self.get_mu(t) / self.get_mu_std(t)

    def get_x2(self, t: int):
        # The formula below was derived using Taylor expansion for x1(t-w) and x1(t - 2w), where w is a window_size
        return (3 * self.get_x1(t) - 4 * self.get_x1(t - self.hyperparams_["window_size"]) + self.get_x1(t - 2 * self.hyperparams_["window_size"])) / (2 * self.hyperparams_["window_size"])

    # Market Pressure normalized by its standard deviation.
    def get_x3(self, t: int):
        return self.get_x3_tilda(t) / self.get_x3_tilda_std(t)

    def get_x4(self, t: int):
        # The formula below was derived using Taylor expansion for x1(t-w) and x1(t - 2w), where w is a window_size
        return (3 * self.get_x3(t) - 4 * self.get_x3(t - self.hyperparams_["window_size"]) + self.get_x3(t - 2 * self.hyperparams_["window_size"])) / (2 * self.hyperparams_["window_size"] * self.get_x3_tilda_std(t))

    # sigma = volatility
    def get_sigma_squared(self, t: int) -> None:
        sigma_squared = 0
        for i in range(t, t - self.hyperparams_["window_size"] - 1, -1):
            sigma_squared += math.log(self.mid_px_[i] /
                                      self.mid_px_[i - 1]) ** 2

        return sigma_squared / self.hyperparams_["window_size"]

    # mu = trend strength
    def get_mu(self, t: int):
        sigma_squared = self.get_sigma_squared(t)
        return sigma_squared / 2 + math.log(self.mid_px_[t] / self.mid_px_[t - self.hyperparams_["window_size"]]) / self.hyperparams_["window_size"]

    def get_mu_std(self, t: int):
        return math.sqrt(self.get_sigma_squared(t) / self.hyperparams_["window_size"])

    # x3_tilda = Market Pressure, but is not normalized by its standard deviation (x3 is).
    def get_x3_tilda(self, t: int):
        total_volume = self.volume_prefix_sum[t] - \
            self.volume_prefix_sum[t - self.hyperparams_["window_size"]]
        if total_volume == 0:
            return 0

        buy_volume = self.buy_volume_prefix_sum[t] - \
            self.buy_volume_prefix_sum[t - self.hyperparams_["window_size"]]
        return 2 * buy_volume / total_volume - 1

    def get_x3_tilda_mean(self, t: int):
        x3_tilda_mean = 0
        for k in range(0, self.hyperparams_["K"]):
            x3_tilda_mean += self.get_x3_tilda(t -
                                               k * self.hyperparams_["window_size"])

        return x3_tilda_mean / (self.hyperparams_["K"])  # sample mean

    def get_x3_tilda_std(self, t: int):
        x3_tilda_mean = self.get_x3_tilda_mean(t)
        x3_tilda_std = 0
        for k in range(0, self.hyperparams_["K"]):
            assert t - k * self.hyperparams_["window_size"] >= 0, print(t, k * self.hyperparams_["window_size"])
            x3_tilda_std += (self.get_x3_tilda(t - k *
                             self.hyperparams_["window_size"]) - x3_tilda_mean) ** 2
        return x3_tilda_std / (self.hyperparams_["K"] - 1)  # sample variance

    def get_delta_PnL(self, t1: int, t2: int):
        pos = self.get_position(t1)
        return pos * (self.mid_px_[t2] - self.mid_px_[t1])

    def get_PnL(self, T: int):
        total_PnL = 0
        for i in range(self.hyperparams_["window_size"] * (self.hyperparams_["K"] + 1), T - self.hyperparams_["tau"], self.hyperparams_["tau"]):
            total_PnL += self.get_delta_PnL(t1=i,
                                            t2=i + self.hyperparams_["tau"])
        return total_PnL

    def get_sortino_ratio(self, T: int):
        total_PnL = 0
        max_drawdown = 0
        for i in range(self.hyperparams_["window_size"] * (self.hyperparams_["K"] + 1), T - self.hyperparams_["tau"], self.hyperparams_["tau"]):
            delta_pnl = self.get_delta_PnL(t1=i, t2=i + self.hyperparams_["tau"])
            total_PnL += delta_pnl
            if total_PnL <= 0:
                max_drawdown = max(max_drawdown, -total_PnL)

        return total_PnL / max_drawdown
