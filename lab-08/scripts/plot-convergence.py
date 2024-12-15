#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import sys

convergence_data = pd.read_csv(sys.argv[1], sep = ",")

plt.rcParams.update({"font.size": 14})

plt.plot(convergence_data.dt,
         convergence_data.eL2,
         marker = 'o',
         label = 'L2')
plt.plot(convergence_data.dt,
         convergence_data.eH1,
         marker = 'o',
         label = 'H1')
plt.plot(convergence_data.dt,
         convergence_data.dt,
         '--',
         label = 'dt')
plt.plot(convergence_data.dt,
         convergence_data.dt**2,
         '--',
         label = 'dt^2')

plt.xscale("log")
plt.yscale("log")
plt.xlabel("dt")
plt.ylabel("error")
plt.legend()

plt.savefig("convergence.pdf")