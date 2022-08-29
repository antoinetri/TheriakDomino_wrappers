# TheriakDomino_wrappers
Some tools to wrap Theriak-Domino routines and plot vectorized phase diagrams

Here are some modest wrapper of the Theriak-Domino software package written by Christian de Capitani (de Capitani & Brown, 1987; de Capitani & Petrakakis, 2010) used to calculate phase equilibria and built phase diagrams.

The python code TheriakDomino_wrapper_plot_pseudosection.py will basically read the .plt file produced by the guzzler app (see an example for a mafic rock in the repo: domplt.plt), read phase diagram parameters and assemblages and produce a plot using matplotlib library for which the color marks the system variance of each space of the diagram. It will also add a label in each space displaying the list of stable phases. Then this plot can be exported as a vectorized file (.eps) and be treated for publication purpose (example below).

![image](https://user-images.githubusercontent.com/14851413/187172630-2010c869-b7eb-4257-8eff-41f9018ed4bd.png)

After a few illustrator editing: 

![image](https://user-images.githubusercontent.com/14851413/187173452-5efa6aee-0b5e-496d-98e5-c7a1348adf9a.png)

The raw result is far from perfection at the moment but saves a lots of time in drawing phases diagrams. This code will be improved soon. All contributions are welcome. Enjoy !
