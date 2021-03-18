#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# In[2]:


# These lines that are commented out for the purpose of testing in jupyter notebook.
'''
model = "kOmegaSST"
yp = 1
'''
model = os.environ["model"]
yp = os.environ["yp"]
if (yp.find('.') == -1):
    yp = int(os.environ["yp"])
else:
    yp = float(os.environ["yp"])


# In[7]:


# This line that is commented out for the purpose of testing in jupyter notebook.
#solverInfo_path = '../run/turbulentFlatPlate1/solverInfo_{}_{}.csv'.format(model, str(yp))
solverInfo_path = 'solverInfo_{}_{}.csv'.format(model, str(yp))
solverInfo_data = pd.read_csv(solverInfo_path, delim_whitespace=True, skiprows = 1)
solverInfo_data.head()


# In[8]:


Ux_init_res = solverInfo_data['Ux_initial']
Uy_init_res = solverInfo_data['Uy_initial']
p_init_res = solverInfo_data['p_initial']
t = solverInfo_data['Time']


# In[12]:


plt.figure()
plt.yscale("log")
plt.plot(t, Ux_init_res, color = 'b', label = "Ux")
plt.plot(t, Uy_init_res, color = 'r', label = "Uy")
plt.plot(t, p_init_res, color = 'g', label = "p")
plt.xlabel("Time")
plt.ylabel("Initial Residual (log scale)")
plt.grid()
plt.legend()
plt.savefig("solverInfo_{}_{}.pdf".format(model, str(yp)))

