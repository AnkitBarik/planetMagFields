#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from .libgauss import get_data
from .magField import getBr
from .plotlib import plotB

class planet:

    def __init__(self,name='earth',datDir='data/'):
    
        self.name   = name.lower()
        self.datDir = datDir
        self.glm, self.hlm, self.lmax, self.idx = \
                get_data(self.datDir,planet=self.name)

        self.p2D, self.th2D, self.Br, self.dipTheta, self.dipPhi = \
                getBr(datDir=self.datDir,planet=self.name,r=1,info=True)

    def plot(self):
        plt.figure(figsize=(16,9))
        plotB(self.p2D,self.th2D,self.Br,planet=self.name)

