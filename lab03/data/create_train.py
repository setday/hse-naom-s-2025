import datetime
import polars as pl
import sys


def process_market_data(input_file: str, output_file: str):
    market_data: pl.DataFrame = pl.read_csv(input_file).drop("ID")

    # Binarize Operation type
    market_data = market_data.with_columns(
        SELL=pl.col("OPER") == "S",
        BUY=pl.col("OPER") == "B",
    ).drop(["OPER"])

    # For the further use
    min_time = market_data.get_column("TIME").min()
    max_time = market_data.get_column("TIME").max()

    # Combine DATE and TIME into a DATETIME column
    market_data = market_data.with_columns(
        # Convert DATE to a Date type (YYYY-MM-DD)
        DATE=pl.col("DATE").cast(pl.String).str.to_date(format="%Y%m%d").cast(pl.Date),
        # Convert TIME to a Time type (HH:MM:SS)
        TIME=pl.col("TIME").str.strptime(pl.Time, format="%H:%M:%S")
    )

    market_data = market_data.with_columns(
        # Combine the two into a DATETIME
        DATETIME=pl.col("DATE").dt.combine(pl.col("TIME")).alias("d1"),
    ).drop(["DATE"])

    market_data = market_data.group_by(pl.col("DATETIME"), maintain_order=True
                                       ).agg(
        (pl.col("VOL") * pl.col("SELL")).sum().alias("SELL_VOLUME"),
        (pl.col("VOL") * pl.col("BUY")).sum().alias("BUY_VOLUME"),
        (pl.col("SELL") * pl.col("LAST")).filter(pl.col("SELL") == True).max().alias("BEST_BID"),
        (pl.col("BUY") * pl.col("LAST")).filter(pl.col("BUY") == True).min().alias("BEST_ASK")
    )

    market_data = market_data.with_columns(
        pl.when((~pl.col("BEST_ASK").is_nan()) & (~pl.col("BEST_BID").is_nan()))
        .then
        ((pl.col("BEST_ASK") + pl.col("BEST_BID")) / 2)
        .otherwise
        (pl.coalesce(pl.col("BEST_ASK"), pl.col("BEST_BID"))).alias("MID_PX")
    )

    # Interpolate missing values
    # Get all timeframes between the min and max time
    time_range = pl.DataFrame({
        "DATETIME": pl.datetime_range(
            start=market_data["DATETIME"].min(), end=market_data["DATETIME"].max(), interval="1s", eager=True
        )
    })

    market_data = time_range.join(market_data, on="DATETIME", how="left")

    market_data = market_data.with_columns([
        pl.col("SELL_VOLUME").fill_null(0),
        pl.col("BUY_VOLUME").fill_null(0),
        pl.col("MID_PX").forward_fill()
    ])

    market_data = market_data.with_columns(
        pl.col("DATETIME").cast(pl.Time).alias("TIME")
    )

    max_time_py = datetime.time(int(max_time[0:2]), int(max_time[3:5]), int(max_time[6:8]))
    min_time_py = datetime.time(int(min_time[0:2]), int(min_time[3:5]), int(min_time[6:8]))

    market_data = market_data.filter(
        (pl.col("TIME") <= max_time_py) & ((pl.col("TIME") >= min_time_py))
    )

    # Save processed data to the output file
    market_data.drop(["DATETIME", "BEST_BID", "BEST_ASK", "TIME"]).write_csv(output_file)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 create_train.py data/<filename>.csv")
    else:
        input_file = sys.argv[1]
        process_market_data(input_file, output_file=f"{input_file[:-4]}_train.csv")
        print(f"file {input_file[:-4]}_train.csv was created!")
