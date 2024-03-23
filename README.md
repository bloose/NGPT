# NGPT
Matlab implementation of the modified Noble Gas Paleothermometer
# Sea Ice Formation, Glacial Melt and the Solubility Pump Boundary Conditions in the Ross Sea
============

Authors
--------
[Brice Loose](https://bloose.github.io)<sup>1</sup>, Sharon Stammerjohn<sup>2</sup> , Peter Sedwick<sup>3</sup> , and Stephen Ackley<sup>4</sup>

1: [URI Graduate School of Oceanography](https://web.uri.edu/gso/), Narragansett, RI, USA.
2: [Institute of Arctic and Alpine Research, CU Boulder](https://www.colorado.edu/instaar/),  Boulder, CO, USA.
3: [Department of Ocean and Earth Sciences, Old Dominion University](https://www.odu.edu/oes) Norfolk, VA, USA.
4: [Center for Advanced Measurements in Extreme Environments, University of Texas at San Antonio](https://www.utsa.edu/NASA-CAMEE/team.html), San Antonio, TX, USA.

Abstract
--------
In-situ sensors for environmental chemistry promise more thorough observations, which are necessary for high-confidence predictions in earth systems science. However, these can be a challenge to interpret, because the sensors are strongly influenced by temperature, humidity, pressure, or other secondary environmental conditions that are not of direct interest. We present a comparison of two statistical learning methods - a Generalized Additive Model, and a Long-Short Term Memory (LSTM) Neural Network model for bias correction of in-situ sensor data. We discuss their performance and tradeoffs, when the two bias correction methods are applied to data from submersible and shipboard mass spectrometers. Both instruments measure the most abundant gases dissolved in water, and can be used to reconstruct biochemical metabolisms, including those that regulate atmospheric carbon dioxide. Both models demonstrate a high degree of skill at correcting for instrument bias using correlated environmental measurements; the difference in their respective performance is less than 1% in terms of root mean squared error. Overall the LSTM bias correction produced an error of 5% for O2 and 8.5% for CO2, when compared against independent membrane DO and laser spectrometer instruments. This represents a predictive accuracy of 92-95% for both gases. It is apparent that the most important factor in a skillful bias correction is the measurement of the secondary environmental conditions that are likely to correlate with the instrument bias. These statistical learning methods are extremely flexible and permit the inclusion of nearly an infinite number of correlates in finding the best bias correction solution.

Status
----------
The paper is in press. Comments, questions, and suggestions are welcome and warmly appreciated. Please email me at bloose@uri.edu.

Code
----
Coding performed in Python.  The GAM backfit algorithm was used for the iterative fit.  The LSTM RNN model was implemented using the Keras interface to Tensorflow.

Data
------


Support
-------
This work was supported by a grant from the National Science Foundation, Award # 1429940.

Acknowledgments
----------------
This research was supported by an award from the National Science Foundation Chemical and Biological Oceanography Program #1429940. We thank two anonymous reviewers for the comments and suggestions that have improved this manuscript. The GAM backfit algorithm is available at https://github.com/bloose/Python_GAM_Backfit. The supplemental contains annotated Python scripts and SWIMS example data to demonstrate application of the GAM and LSTM to bias correction.

**Thanks** to Cesar Rocha (crocha@ucsd.edu) for providing this template and example to follow.
