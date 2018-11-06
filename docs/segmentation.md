# Segmentation

Various functions for segmenting cells. All use connected areas to create distinct objects

`thresholdSegment`

Defines a pixel threshold using various different methods.

`edgeSegment`

Detects edges in the image and dilates to fill in regions in between

`watershedSegment`

Requires input segmentation. Applies a watershed transform to the image, nucleating from the seedMask.
