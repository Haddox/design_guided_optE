* Each CSV file gives stability scores predicted by a random-forest model trained on a given miniprotein topology. File names indicate topology as: `Rocklin.{topology}.rfr.20.sqrt.False.113.csv`.
* In each CSV, each row corresponds to a design, and columns indicate:
	* `name`: the name of the design
	* `stability_pred`: the stability score predicted by the model
