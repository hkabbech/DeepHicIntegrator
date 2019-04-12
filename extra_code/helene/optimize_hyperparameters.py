#! /usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import print_function

import os
import random as rd
import numpy as np
import cooler
from keras.models import load_model
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap

from src.hic import Hic 
from src.predict_hic import PredictHic

from contextlib import redirect_stdout
from sklearn.model_selection import train_test_split
from keras.layers import Conv2D, MaxPooling2D, UpSampling2D
from keras.models import Model
from keras.optimizers import RMSprop
from keras.layers import Input
from keras.models import Sequential

from hyperopt import Trials, STATUS_OK, tpe
from hyperas import optim
from hyperas.distributions import choice, uniform


def data():
    '''
    Data providing function:

    Make sure to have every relevant import statement included here and return data as
    used in model function below. This function is separated from model() so that hyperopt
    won't reload data for each evaluation run.
    '''
    N = 60
    CHR_TRAIN = 10
    RESOLUTION = 25000
    HIC_FILE = 'data/1_binaries/hic/GSE63525_HUVEC_combined_30.hic'
    COOL_FILE = HIC_FILE.split('.hic')[0]+'_'+str(RESOLUTION)+'.cool'
    if not os.path.exists(COOL_FILE):
        hic2cool_convert(HIC_FILE, COOL_FILE, RESOLUTION)
    TRAIN = Hic(cooler.Cooler(COOL_FILE), CHR_TRAIN, N)
    TRAIN.set_matrix()
    TRAIN.set_sub_matrices()
    return train_test_split(TRAIN.sub_matrices, TRAIN.sub_matrices, test_size=0.2, random_state=13)


def model(X_train, Y_train, X_test, Y_test):
    '''
    Model providing function:

    Create Keras model with double curly brackets dropped-in as needed.
    Return value has to be a valid python dictionary with two customary keys:
        - loss: Specify a numeric evaluation metric to be minimized
        - status: Just use STATUS_OK and see hyperopt documentation if not feasible
    The last one is optional, though recommended, namely:
        - model: specify the model just created so that we can later use it again.
    '''
    N = 60
    model = Sequential()
    model.add(Conv2D(32, (3, 3), activation='relu', padding='same', input_shape=Input(shape=(N, N, 1))))
    model.add(MaxPooling2D(pool_size=(2, 2)))
    model.add(Conv2D(64, (3, 3), activation='relu', padding='same'))
    model.add(MaxPooling2D(pool_size=(2, 2)))
    model.add(Conv2D(128, (3, 3), activation='relu', padding='same'))
    model.add(Conv2D(128, (3, 3), activation='relu', padding='same'))
    model.add(UpSampling2D((2, 2)))
    model.add(Conv2D(64, (3, 3), activation='relu', padding='same'))
    model.add(UpSampling2D((2, 2)))
    model.add(Conv2D(1, (3, 3), activation='sigmoid', padding='same'))

    model.compile(loss='mean_squared_error', optimizer=RMSprop(),
                  metrics=['accuracy'])

    model.fit(X_train, Y_train,
              batch_size={{choice([64, 128])}},
              nb_epoch=10,
              show_accuracy=True,
              verbose=2,
              validation_data=(X_test, Y_test))

    score, acc = model.evaluate(X_test, Y_test, show_accuracy=True, verbose=0)
    print('Test accuracy:', acc)
    return {'loss': -acc, 'status': STATUS_OK, 'model': model}

if __name__ == '__main__':

    best_run, best_model = optim.minimize(model=model,
                                          data=data,
                                          algo=tpe.suggest,
                                          max_evals=5,
                                          trials=Trials())
    X_train, Y_train, X_test, Y_test = data()
    print("Evalutation of best performing model:")
    print(best_model.evaluate(X_test, Y_test))