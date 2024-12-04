## Download Market Data

You can download market data from one of the following sources:

1. **FINAM** (no registration required):  
   [https://www.finam.ru/quote/moex/gazp/export/](https://www.finam.ru/quote/moex/gazp/export/)

2. **BI Angelite** (requires quick registration, but recommended for better data):  
   [https://tickwize.angelitrade.com/](https://tickwize.angelitrade.com/)

### Data Format

The downloaded data MUST have the following format:

```csv
DATE,TIME,LAST,VOL,ID,OPER
20241101,09:59:44,123.8,10,11370547379,S
```

### Next Steps

1. After downloading the file, place it inside the directory:  
   `lab03/data/<filename>.csv`
   
2. Run the following command to generate training data:

   ```bash
   python3 create_train.py <filename>.csv