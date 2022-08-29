# TheriakDomino_wrappers
Some tools to wrap Theriak-Domino routines and plot vectorized phase diagrams

Here are some modest wrapper of the Theriak-Domino software package written by Christian de Capitani (de Capitani & Brown, 1987; de Capitani & Petrakakis, 2010) used to calculate phase equilibria and built phase diagrams.

The python code TheriakDomino_wrapper_plot_pseudosection.py will basically read the .plt file produced by the guzzler app, read phase diagram parameters and assemblages and produce a plot using matplotlib library for which the color marks the system variance of each space of the diagram. It will also add a label in each space displaying the list of stable phases. Then this plot can be exported as a vectorized file (.eps) and be treated for publication purpose (example below).

![image](https://user-images.githubusercontent.com/14851413/187167882-9bdc29ae-de51-4c13-a38b-c5f0bb054370.png)

After a few illustrator editing: 

The raw result is far from perfection at the moment but saves a lots of time in drawing phases diagrams. This code will be improved soon. All contributions are welcome. Enjoy !
