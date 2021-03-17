#!/usr/bin/env python
# coding: utf-8

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# These lines that are commented out for the purpose of testing in jupyter notebook.
'''
UInf = 69.4
nuInf = 1.388e-05
model = "kOmegaSST"
yp = 1
'''
UInf = float(os.environ["UInf"])
nuInf = float(os.environ["nuInf"])
model = os.environ["model"]
yp = os.environ["yp"]
if (yp.find('.') == -1):
    yp = int(os.environ["yp"])
else:
    yp = float(os.environ["yp"])

# This line that is commented out for the purpose of testing in jupyter notebook.
#tauw_path = '../run/turbulentFlatPlate/tauw_{}_{}.csv'.format(model, str(yp))
tauw_path = 'tauw_{}_{}.csv'.format(model, str(yp))
tauw_data = pd.read_csv(tauw_path, delim_whitespace=True)
tauw_data.head()

x0 = 0
Rex = (tauw_data['ccx'] - x0)*UInf/nuInf
Cf = np.sqrt(tauw_data['tau_xx']*tauw_data['tau_xx'] + tauw_data['tau_yy']*tauw_data['tau_yy'] + tauw_data['tau_zz']*tauw_data['tau_zz'])/(0.5*UInf**2)


wieghardt = 0.288*(np.log10(Rex))**(-2.45)

plt.figure()
plt.ylim([0, 0.006])
plt.plot(Rex, Cf, color = 'b', label = "{} y+ = {}".format(model, str(yp)))
plt.plot(Rex, wieghardt, color = 'r', label = "Wieghardt")
plt.xlabel("Re_x")
plt.ylabel("Cf")
plt.grid()
plt.legend()
plt.savefig("tauw_{}_{}.pdf".format(model, str(yp)))


