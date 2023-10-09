### This folder contains all the analysis performed on streamflow time series. 

`ts_info.csv` contains the name of each gauge, the year range, number of years, mean flow, flow range, time lag, dimension embedding, rolling window, and global average LLE. 


All `"gauge name"_ts.csv` files are in the same format 

 | year | flow | signal | avgLLE | rolling average | rolling variance |
 | ---- | ---- | ------ | ------ | --------------- | ---------------- |

Where the flow is the reconstructed and observed streamflow combined. The rolling average and rolling variance are calculated based on the rolling window in the `ts_info.csv` file. 

`high-low_pred_years.csv` contains each gauge year and average LLE values where the average LLE value falls below -0.05 (high predictability threshold) and above 0.05 (low predictability threshold). The year can then be used to create windows around the average LLE values for interpretation. 
 - Note: The NA values are present by design. In the algorithm written only 15 values were stored based on the thresholds, if less than 15 values were found, then value filled with NA. 