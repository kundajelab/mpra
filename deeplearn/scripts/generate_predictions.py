import keras
from keras.models import model_from_json, model_from_yaml
import json

modelID = "record_2_model_Yjv2n_"

json_path = "../model_files/regressionJun24Positives/" + modelID + "modelJson.json"
with open(json_path) as json_file:
    json_string = json.dumps(json.load(json_file))
    model = model_from_json(json_string)  

model.load_weights("../model_files/regressionJun24Positives/" + modelID + "modelWeights.h5")

import h5py
import numpy as np

with h5py.File('../hdf5files/regressionJun24Positives/train_data.hdf5') as f:
    print f.keys()
    x_data = f['X']
    y_data = f['Y']
    x_test = np.array(x_data['sequence'])
    y_test = np.array(y_data['output'])

y_pred = model.predict(x_test)

np.savetxt('../predictions/regressionJun24Positives_trainPredictions.txt', y_pred, delimiter='\t')

with h5py.File('../hdf5files/regressionJun24Positives/valid_data.hdf5') as f:
    print f.keys()
    x_data = f['X']
    y_data = f['Y']  
    x_test = np.array(x_data['sequence'])
    y_test = np.array(y_data['output'])

y_pred = model.predict(x_test)

np.savetxt('../predictions/regressionJun24Positives_validPredictions.txt', y_pred, delimiter='\t')
