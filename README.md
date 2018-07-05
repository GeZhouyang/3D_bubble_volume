## Brief description of the code

Given a 3D image of a periodic bubble dispersion, how to count the actual number of bubbles inside and get the individual volume of each one? Here, a periodic dispersion means the organization/shape of the bubbles has periodicity in the spatial directions, and our 3D image is colorless, i.e. it is an intensity field (or local volume fraction).

To solve this problem, we divide the task into two parts.

1. We count and label the bubbles from the raw input picture, using a n-dimensional image processing library `scipy/ndimage`. The function `ndimage.label` takes two inputs, counting condition and neighborhood, and returns the connected regions tagged with a unique label. See the related pages for details.

2. We account for periodicity by merging bubbles that are actually connected. This is the main part of the present code. Essentially, we apply the same labeling algorithm to each 2D plane, compare its opposite directions using labels obtained in the 3D field, then merge if multiple bubbles intersect the same plane. The detailed algorithm is ommitted since the implementation is specific to one type of periodic boundary conditions. It may look a bit tedious; but the methodology is general, and the solution is so-far verified to be complete.

This script is coded during the 2018 CTR summer program at Stanford University. I thank Michael Dodd for mentioning the `ndimage` library. I also thank Marco Rosti for providing the data set to verify the algorithm.