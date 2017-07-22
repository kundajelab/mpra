import avutils.util as util
import avutils.file_processing as fp
import numpy as np

def get_other_data_loaders(train_file_prefix):
    return {
        "train": {
            "class": "hdf5_data_loader.MultimodalAtOnceDataLoader",
            "kwargs": {
                "batch_size": 500,
                "path_to_hdf5": "../hdf5files/sharpr_regression_znormed_jul18/train_data.hdf5",
                "num_to_load_for_eval": 100000,
                "bundle_x_and_y_in_generator": False,
                "strip_enclosing_dictionary": True
            }
        }
    }


def get_model_trainer(seed):
    return {
        "class": "keras_model_trainer.KerasFitGeneratorModelTrainer",
        "kwargs": {
            "seed": seed,
            "samples_per_epoch": 200000,
            "stopping_criterion_config": {
                "class": "EarlyStopping" ,
                "kwargs": {
                   "max_epochs": 1000, 
                   "epochs_to_wait": 14,
                   "running_mean_over_epochs": 1
                } 
            },
            # "class_weight": {"0":1, "1":pos_weight}
        }
    }


def get_model_creator(nconvlayers, # number of conv layers
                      filter_lengths, # list of length len(nconvlayers) with filter widths
                      p_conv_dropout, # dropout probability for conv layers
                      pool_length, # max pooling width
                      stride, # max pooling stride
                      nb_filters, # number of filters per layer
                      separableFC, # boolean indicating whether to use separable FC layers
                      symmetric,
                      output_dim, # number of tasks - outputs in the final dense layer
                      lr # learning rate
                     ):

    conv_block = []
    assert nconvlayers <= 4
    #  if (nlayers == 3):
        #  filter_lengths = [15, 14, 14]     
    #  elif (nlayers == 2):
        #  filter_lengths = [21, 21]
    #  elif (nlayers == 1):
        #  filter_lengths = [41]

    filter_lengths = filter_lengths

    for (i, filter_length) in enumerate(filter_lengths): 
        conv_block.extend([
            {
                "class": "keras.layers.convolutional.Convolution1D", 
                "kwargs": {
                    "input_shape": [145,4],
                    "nb_filter": nb_filters[i],
                    "filter_length": filter_length,
                }
            },
            {
                "class": "keras.layers.normalization.BatchNormalization",
                "kwargs": {}
            },
            {
                "class": "keras.layers.core.Activation", 
                "kwargs": {"activation": "relu"}
            },
            {
               "class": "keras.layers.core.Dropout",
               "kwargs": {
                   "p": p_conv_dropout
               }
            },
            {
               "class": "keras.layers.convolutional.MaxPooling1D",
               "kwargs": {
                   "pool_length": pool_length,
                   "stride": stride
               }
            }
            ])

    return {
        "class": "flexible_keras.FlexibleKerasSequential",
        "kwargs": {
            "layers_config": conv_block+[
                #  {
                    #  "class": "keras.layers.convolutional.MaxPooling1D", 
                    #  "kwargs": {"pool_length": pool_length, "stride": stride}
                #  },
                ({
                    "class": "keras.layers.convolutional.SeparableFC",
                    "kwargs": { "output_dim": separableFC,
                                "symmetric": symmetric
                               }
                 }
                 if separableFC == True else 
                 {
                     "class": "keras.layers.core.Flatten", 
                     "kwargs": {}
                 }),
                {
                    "class": "keras.layers.core.Dense", 
                    "kwargs": {"output_dim": output_dim}
                },
                {
                    "class": "keras.layers.core.Activation", 
                    "kwargs": {"activation": "linear"}
                }
            ],
            "optimizer_config": {
                "class": "keras.optimizers.Adam",
                "kwargs": {"lr": lr}
            },
            "loss": "mean_squared_error" 
        } 
    }


def get_hyperparameter_configs(prefix, nconvlayers, filter_lengths,
                               p_conv_dropout, pool_length, stride,
                               nb_filters, separableFC,
                               symmetric, output_dim,
                               lr, seed, train_file_prefix=""):
    message = (prefix
               +("-"+train_file_prefix if (len(train_file_prefix) > 0) else "")
               +"-nconvlay_"+str(nconvlayers)
               +"-filterlens_" + ','.join([str(fl) for fl in filter_lengths])
               +"-pconvdropout_" + str('%.3f' % p_conv_dropout)
               +"-poolw_" + str(pool_length)
               +"-poolstride_" + str(stride)
               +"-nbf_" + ','.join([str(nf) for nf in nb_filters])
               +"-sep_"+('t' if separableFC else 'f')
               +"-dim_"+(str(separableFC) if separableFC else "null")
               +"-symm_"+('t' if symmetric else 'f')
               +"-lr_" + str(lr)
               +"-seed_" + str(seed))
    return {
        "message": message,
        "other_data_loaders": get_other_data_loaders(train_file_prefix=train_file_prefix),
        "model_creator": get_model_creator(nconvlayers=nconvlayers, 
                                           filter_lengths=filter_lengths,
                                           p_conv_dropout=p_conv_dropout,
                                           pool_length=pool_length,
                                           stride=stride,
                                           nb_filters=nb_filters,
                                           separableFC=separableFC, 
                                           symmetric=symmetric,
                                           output_dim=output_dim,
                                           lr=lr),
        "model_trainer": get_model_trainer(seed=seed)
    }


def get_random_hyperparams():
    nconvlayers = np.random.randint(low = 3, high = 4)
    filter_lengths = np.random.randint(low = 4, high = 10, size = nconvlayers)
    p_conv_dropout = np.random.random() / 6 # float from 0 to 0.5
    pool_length = np.random.randint(low = 1, high = 2)
    stride = pool_length
    nb_filters = np.random.randint(low = 75, high = 200, size = nconvlayers)
    separableFC = np.random.random() > 1 # never use sepFC layer / random boolean with p = 1/2
    if separableFC == True: # if True, then define output_dim
        separableFC = np.random.randint(low = 25, high = 100)
    symmetric = np.random.random() > 1 # never use sep FC layer / p = 1/2
    lr = 10 ** -np.random.uniform(low=3, high=4)
    good_params = check_params(filter_lengths = filter_lengths, pool_length = pool_length, stride = stride)
    if good_params:
        return {"nconvlayers": nconvlayers,
                "filter_lengths": filter_lengths,
                "p_conv_dropout": p_conv_dropout,
                "pool_length": pool_length,
                "stride": stride,
                "nb_filters": nb_filters,
                "separableFC": separableFC,
                "symmetric": symmetric,
                "lr": lr}
    else:
        return {}

def check_params(filter_lengths, pool_length, stride, input_length=145):
    length = input_length
    for conv in filter_lengths:
        length = length - conv + 1
        length = (length - pool_length) / stride
    if length < 1:
        return False
    return True

def main(args):
    seed = 723
    np.random.seed(seed)
    
    possible_settings = [get_random_hyperparams() for i in range(25)]
    # possible_settings[0]["filter_lengths"] = np.array([4, 4, 4])
    possible_settings[1]["filter_lengths"] = np.array([4, 5, 6])
    possible_settings[2]["filter_lengths"] = np.array([6, 6, 6])
    possible_settings[3]["filter_lengths"] = np.array([4, 6, 8])
    possible_settings[4]["filter_lengths"] = np.array([4, 7, 10])
    possible_settings[5]["filter_lengths"] = np.array([8, 8, 8])
    possible_settings[1]["filter_lengths"] = np.array([5, 8, 5])
    possible_settings[1]["filter_lengths"] = np.array([5, 10, 4])
    possible_settings = [setting for setting in possible_settings if setting != {}]

    hyperparameter_configs = []
    for seed in range(1):
        for settings in possible_settings: 
            hyperparameter_configs.append(
                get_hyperparameter_configs(prefix=args.prefix, seed=seed, output_dim = 12, **settings)
               )
            
    fp.write_to_file("../yamls/" + args.prefix + "_config/hyperparameter_configs_list_paramsearch.yaml",
                     util.format_as_json(hyperparameter_configs))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("prefix")
    main(parser.parse_args())
