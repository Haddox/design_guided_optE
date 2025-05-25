* `outlier_designs_data_processing.ipynb`: notebook for processing data on miniprotein designs from Rocklin et al. before fitting random-forest model
* `outlier_designs_rfr_training.ipynb`: notebook for training random-forest models on processed data

* `data/experimental_stability_scores/`: experimentally measured stability scores from Rocklin et al.
	* `Rocklin.v6.experimental_stability_scores.csv`: experimental stability scores from Rocklin et al. The column `stabilityscore` gives the original stability scores measured by Rocklin et al. The column `stabilityscore_cnn_calibrated` gives stability scores recomputed using the updated unfolded state model from Singer et al.

* `data/processed_data/`: the JSON files in this directory list the full set of biophysical features used to train the random-forest model for a given miniprotein topology

* See [zenodo link] for a CSV file with values of all biophysical features for each miniprotein design.

* `data/rfrdisparity/`: directory with CSV files reporting the predicted stability score and experimentally measured stability score of each design. There is one file per topology (e.g., `Rocklin.HHH.rfrdisparity.v1.1527.csv`) with data for all designs from that topology. There is also one file that only shows data for the 21 outlier designs characterized by DMS (`outlier_designs.csv`). In each CSV, each row corresponds to a design, and columns indicate:
	* `name`: the name of the design
	* `stability_pred`: the stability score predicted by the model
	* `stabilityscore_cnn_calibrated`: the experimentally measured stability score used to train the model
