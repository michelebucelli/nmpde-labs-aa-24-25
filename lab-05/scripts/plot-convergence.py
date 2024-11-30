#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import sys

convergence_data = pd.read_csv(sys.argv[1], sep = ",")

plt.rcParams.update({"font.size": 14})

plt.plot(convergence_data.h,
         convergence_data.eL2,
         marker = 'o',
         label = 'L2')
plt.plot(convergence_data.h,
         convergence_data.eH1,
         marker = 'o',
         label = 'H1')
plt.plot(convergence_data.h,
         convergence_data.h,
         '--',
         label = 'h')
plt.plot(convergence_data.h,
         convergence_data.h**2,
         '--',
         label = 'h^2')
plt.plot(convergence_data.h,
         convergence_data.h**3,
         '--',
         label = 'h^3')

plt.xscale("log")
plt.yscale("log")
plt.xlabel("h")
plt.ylabel("error")
plt.legend()

plt.savefig("convergence.pdf")