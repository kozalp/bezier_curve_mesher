# bezier_curve_mesher

Hello! The ultimate goal of this project is to write a code that meshes a Bezier curve.
To do that, I used the same approach as in the paper below by Dr. Michael D. Senter.

https://iopscience.iop.org/article/10.1088/1748-3190/ababb0/pdf?casa_token=9cPbFTHHpkoAAAAA:WwSZzRXZmBxv-6IVI3gMkWqwbTxL88IrEMDuAE6hvp8gQVWhVkhitQ9hukRQIEheCZSwH_klCg

The workflow is as follows:

1. Run curve.m and select weight points for the Bézier curve.
2. Once done, click enter to plot the Bézier curve.
3. Press escape and terminate curve.m
4. Run mesher.m. By default, mesher.m separates the curve into 20 roughly equal pieces. 

Currently, the points are not roughly equal. I am suspecting that the problem is in the calculation of the cost function.
