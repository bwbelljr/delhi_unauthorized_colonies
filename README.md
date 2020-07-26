# Delhi Unauthorized Colonies Data Cleaning

This repository includes code to de-duplicate unauthorized colonies of Delhi extracted from the Delhi Development Authority [website](https://dda.org.in/ddaweb/planningmap.aspx).

Specifically, for each unauthorized colony there is a map with that colony and its neighboring colonies present. When these polygons are extracted from the georeferenced PDF, multiple polygons refers to the same map. This code deduplicates the dataset by selecting the polygon whose centroid is closest to the centroid of the bounding box for the PDF.
