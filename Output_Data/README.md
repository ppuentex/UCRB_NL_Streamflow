### This folder contains all the analysis performed on streamflow time series. 

`ts_info.csv` contains the name of each gauge, the year range, number of years, mean flow, flow range, time lag, dimension embedding, rolling window, and global average LLE. 


All `"gauge name"_ts.csv` files are in the same format 

 | year | flow | signal | avgLLE | rolling average | rolling variance |
 | ---- | ---- | ------ | ------ | --------------- | ---------------- |

Where the flow is the reconstructed and observed streamflow combined. The rolling average and rolling variance are calculated based on the rolling window in the `ts_info.csv` file. 